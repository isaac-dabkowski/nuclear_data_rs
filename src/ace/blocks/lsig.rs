// Represents the LSIG data block - contains locations of incident neutron cross section values.
use crate::ace::arrays::{NxsArray, JxsArray};
use crate::ace::blocks::DataBlockType;

// See page 16 of the ACE format spec for a description of the LSIG block
#[derive(Debug, Clone, PartialEq)]
pub struct LSIG {
    pub xs_locs: Vec<usize>
}

impl LSIG {
    pub fn process(data: &[f64]) -> Self {
        let xs_locs: Vec<usize> = data
            .iter()
            .map(|val| val.to_bits() as usize)
            .collect();

        Self { xs_locs }
    }

    pub fn pull_from_xxs_array<'a>(nxs_array: &NxsArray, jxs_array: &JxsArray, xxs_array: &'a [f64]) -> &'a [f64] {
        // Block start index (binary XXS is zero indexed for speed)
        let block_start = jxs_array.get(&DataBlockType::LSIG) - 1;
        // Calculate the block end index, see the LSIG description in the ACE spec
        let num_reactions = nxs_array.ntr;
        let mut block_end = block_start + num_reactions;
        // Avoid issues if this is the last block in the file
        if block_end == xxs_array.len() + 1 {
            block_end -= 1;
        }
        // Return the block
        &xxs_array[block_start..block_end]
    }
}

impl std::fmt::Display for LSIG {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LSIG({} xs)", self.xs_locs.len())
    }
}

#[cfg(test)]
mod tests {
    use crate::ace::utils::get_parsed_test_file;

    #[tokio::test]
    async fn test_lsig_parsing() {
        let parsed_ace = get_parsed_test_file().await;

        // Check contents
        let lsig = parsed_ace.data_blocks.LSIG.unwrap();
        assert_eq!(lsig.xs_locs, vec![1]);
    }
}