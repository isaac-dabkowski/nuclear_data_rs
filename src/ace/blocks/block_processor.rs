use std::error::Error;
use std::fs::File;
use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};

use strum::IntoEnumIterator;

use crate::ace::blocks::{
    DataBlockType,
    DataBlock,
    ESZ,
    MTR,
    LSIG,
    SIG,
};
use crate::ace::arrays::{JxsArray, NxsArray};
use crate::async_task_dag::{AsyncTaskDag, Task, TaskResults, GetTaskResult};

#[derive(Clone, Debug, Default)]
pub struct DataBlocks {
    pub ESZ: Option<ESZ>,
    pub MTR: Option<MTR>,
    pub LSIG: Option<LSIG>,
    pub SIG: Option<SIG>
}

impl DataBlocks {
    // Create a new BlockProcessor from a XXS array, the NXS and JXS array are used to
    // determine the start and end locations of each block
    pub async fn from_ascii_file(reader: &mut BufReader<File>, nxs_array: &NxsArray, jxs_array: &JxsArray) -> Result<Self, Box<dyn Error>> {
        // Read the entire XXS array into a vector, which we will then partition into the blocks
        let mut xxs_array: Vec<&str> = Vec::with_capacity(nxs_array.xxs_len + 1);
        let time = std::time::SystemTime::now();

        let mut buffer = String::new();
        reader.read_to_string(&mut buffer)?;

        xxs_array.push("INDEX PLACEHOLDER"); // 1-based index

        for line in buffer.lines() {
            for i in 0..4 {
                let start = i * 20;
                if start >= line.len() {
                    break;
                }
                let end = usize::min(start + 20, line.len());
                xxs_array.push(line[start..end].trim_ascii_start());
            }
        }
        let xxs_array = xxs_array.to_vec();

        // // Add a dummy entry to make XXS 1-indexable to match the ACE spec better.
        // xxs_array.push("INDEX PLACEHOLDER".to_string());

        // for line in reader.lines() {
        //     let line = line?;
        //     // Split the line into chunks of 20 characters
        //     for chunk in line.as_bytes().chunks(20) {
        //         let block = String::from_utf8_lossy(chunk).trim_start().to_string();
        //         xxs_array.push(block);
        //     }
        // }

        println!(
            "⚛️  Time to read in XXS array ⚛️ : {} ms",
            std::time::SystemTime::now().duration_since(time).unwrap().as_millis()
        );
        let time = std::time::SystemTime::now();
        
        // Split XXS array into raw text correspoding to each block
        let block_map = DataBlocks::split_ascii_xxs_into_blocks(nxs_array, jxs_array, &xxs_array);
        // println!("{:?}", block_map);
        println!(
            "⚛️  Time to split XXS into blocks ⚛️ : {} ms",
            std::time::SystemTime::now().duration_since(time).unwrap().as_millis()
        );
        let time = std::time::SystemTime::now();

        // Build an AsyncTaskDag to process all of our blocks
        let dag = DataBlocks::construct_dag(block_map, nxs_array);
        println!(
            "⚛️  Time to construct DAG ⚛️ : {} ms",
            std::time::SystemTime::now().duration_since(time).unwrap().as_millis()
        );
        let time = std::time::SystemTime::now();

        // Execute the DAG
        dag.execute().await.unwrap();
        println!(
            "⚛️  Time to execute DAG ⚛️ : {} ms",
            std::time::SystemTime::now().duration_since(time).unwrap().as_millis()
        );
        let time = std::time::SystemTime::now();

        // Pass the DAG results back onto our DataBlocks object
        let data_blocks = DataBlocks::from_dag_results(dag);
        Ok( data_blocks )
    }

    fn split_ascii_xxs_into_blocks<'a>(nxs_array: &NxsArray, jxs_array: &JxsArray, xxs_array: &'a Vec<&'a str>) -> HashMap<DataBlockType, &'a [&'a str]> {
        let mut block_map: HashMap<DataBlockType, &'a [&'a str]> = HashMap::default();
        // Loop over all possible DataBlockTypes
        for block_type in DataBlockType::iter() {
            // If the block type's start index is non-zero, the block is present in the XXS array
            let start_index = jxs_array.get(&block_type);
            if start_index != 0 {
                // Pull the block from the XXS array (if procedure to do so has been implemented)
                if let Some(block_text) = DataBlocks::pull_block_from_ascii_xxs_array(&block_type, nxs_array, jxs_array, xxs_array) {
                    block_map.insert(block_type, block_text);
                }
            }
        }
        block_map
    }

    fn pull_block_from_ascii_xxs_array<'a>(block_type: &DataBlockType, nxs_array: &NxsArray, jxs_array: &JxsArray, xxs_array: &'a [&'a str]) -> Option<&'a [&'a str]> {
        match block_type {
            DataBlockType::ESZ => Some(ESZ::pull_from_ascii_xxs_array(nxs_array, jxs_array, xxs_array)),
            DataBlockType::MTR => Some(MTR::pull_from_ascii_xxs_array(nxs_array, jxs_array, xxs_array)),
            DataBlockType::LSIG => Some(LSIG::pull_from_ascii_xxs_array(nxs_array, jxs_array, xxs_array)),
            DataBlockType::SIG => Some(SIG::pull_from_ascii_xxs_array(nxs_array, jxs_array, xxs_array)),
            _ => {
                // println!("DataBlockType {} was found in XXS array, but its parsing has not been implemented yet!", block_type);
                None
            }
        }
    }

    // Build a DAG for block processing based on what blocks are present
    fn construct_dag(block_map: HashMap<DataBlockType, &[&str]>, nxs_array: &NxsArray) -> AsyncTaskDag<DataBlockType, DataBlock> {
        let mut dag: AsyncTaskDag<DataBlockType, DataBlock> = AsyncTaskDag::new();
        let nxs_array = nxs_array.clone();

        // Energy grid
        let esz_text = block_map.get(&DataBlockType::ESZ).unwrap().iter().map(|s| s.to_string()).collect::<Vec<String>>().clone();
        let esz_closure = {
            let nxs_array = nxs_array.clone();
            move |_| async move {
                Ok(DataBlock::ESZ(ESZ::process(esz_text.clone(), &nxs_array)))
            }
        };
        let esz_task = Task::new(DataBlockType::ESZ, esz_closure);
        let esz_task_id = dag.add_task(esz_task);

        // Reaction MT values
        let mtr_text = block_map.get(&DataBlockType::MTR).unwrap().iter().map(|s| s.to_string()).collect::<Vec<String>>().clone();
        let mtr_closure = {
            move |_| async move {
                Ok(DataBlock::MTR(MTR::process(mtr_text.clone())))
            }
        };
        let mtr_task = Task::new(DataBlockType::MTR, mtr_closure);
        let mtr_task_id = dag.add_task(mtr_task);

        // Cross section locations
        let lsig_text = block_map.get(&DataBlockType::LSIG).unwrap().iter().map(|s| s.to_string()).collect::<Vec<String>>().clone();
        let lsig_closure = {
            move |_| async move {
                Ok(DataBlock::LSIG(LSIG::process(lsig_text.clone())))
            }
        };
        let lsig_task = Task::new(DataBlockType::LSIG, lsig_closure);
        let lsig_task_id = dag.add_task(lsig_task);

        // Cross section values
        let sig_text = block_map.get(&DataBlockType::SIG).unwrap().iter().map(|s| s.to_string()).collect::<Vec<String>>().clone();
        let sig_closure = {
            move |results: TaskResults<DataBlockType, DataBlock>| async move {
                let esz = match results.get_result(&DataBlockType::ESZ)? {
                    DataBlock::ESZ(val) => val,
                    _ => panic!("ESZ block was likely improperly parsed!")
                };
                let mtr = match results.get_result(&DataBlockType::MTR)? {
                    DataBlock::MTR(val) => val,
                    _ => panic!("MTR block was likely improperly parsed!")
                };
                let lsig = match results.get_result(&DataBlockType::LSIG)? {
                    DataBlock::LSIG(val) => val,
                    _ => panic!("LSIG block was likely improperly parsed!")
                };
                Ok(DataBlock::SIG(SIG::process(sig_text.clone(), mtr, lsig, esz)))
            }
        };
        let sig_task = Task::new(DataBlockType::SIG, sig_closure);
        let sig_task_id = dag.add_task(sig_task);
        dag.add_task_dependency(esz_task_id, sig_task_id).unwrap();
        dag.add_task_dependency(mtr_task_id, sig_task_id).unwrap();
        dag.add_task_dependency(lsig_task_id, sig_task_id).unwrap();

        dag
    }

    // Construct DataBlocks from results of a DAG
    fn from_dag_results(dag: AsyncTaskDag<DataBlockType, DataBlock>) -> Self {
        let mut data_blocks = DataBlocks::default();
        for result in dag.get_all_results().iter() {
            let (block_type, block_value) = result.pair();
            match (block_type, block_value) {
                (DataBlockType::ESZ, DataBlock::ESZ(esz)) => data_blocks.ESZ = Some(esz.clone()),
                (DataBlockType::MTR, DataBlock::MTR(mtr)) => data_blocks.MTR = Some(mtr.clone()),
                (DataBlockType::LSIG, DataBlock::LSIG(lsig)) => data_blocks.LSIG = Some(lsig.clone()),
                (DataBlockType::SIG, DataBlock::SIG(sig)) => data_blocks.SIG = Some(sig.clone()),
                _ => println!("Block type {} has been processed but is not passed back onto DataBlocks!", block_type),
            }
        }
        data_blocks
    }
}

impl std::fmt::Display for DataBlocks {
    fn fmt(&self, _: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}
