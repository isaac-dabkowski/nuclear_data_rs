use std::collections::HashMap;
use std::ops::Deref;

use rayon::prelude::*;

use crate::arrays::Arrays;
use crate::blocks::block_traits::{Process, PullFromXXS, block_range_to_slice, get_block_start};
use crate::blocks::{BlockType, ESZ, LSIG, MTR};
use crate::helpers::reaction_type_from_MT;

//=====================================================================
// SIG data block
//
// Contains incident neutron cross section data for the ACE file. See
// the ACE format spec for a description of the SIG block.
//=====================================================================
#[derive(Debug, Clone)]
pub struct SIG(pub CrossSectionMap);

impl Deref for SIG {
    type Target = CrossSectionMap;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> PullFromXXS<'a> for SIG {
    fn pull_from_xxs_array(arrays: &'a Arrays) -> Option<&'a [f64]> {
        // We expect SIG if NXS(4) (NTR) != 0
        let has_xs_other_than_elastic = arrays.nxs.ntr != 0;

        // Get the starting index of the block in the XXS array
        let block_start = get_block_start(
            &BlockType::SIG,
            arrays,
            has_xs_other_than_elastic,
            "SIG is expected if NXS(4) (NTR) != 0, but SIG was not found.".to_string(),
        )?;

        // Calculate the block length, see the SIG description in the ACE spec
        // Loop over the number of cross sections
        let mut block_length: usize = 1;
        for _ in 0..arrays.nxs.ntr {
            // Get the number of energy points in the cross section
            let num_entries = arrays.xxs[block_start + block_length].to_bits() as usize;
            // Jump forward to the next cross section
            block_length += num_entries + 2;
        }

        // Return the block's raw data as a slice
        Some(block_range_to_slice(block_start, block_length, arrays))
    }
}

impl<'a> Process<'a> for SIG {
    type Dependencies = (&'a Option<MTR>, &'a Option<LSIG>, &'a Option<ESZ>);

    fn process(
        data: &[f64],
        _arrays: &Arrays,
        dependencies: (&Option<MTR>, &Option<LSIG>, &Option<ESZ>),
    ) -> Self {
        let mtr = dependencies
            .0
            .as_ref()
            .expect("MTR is required to construct SIG");
        let lsig = dependencies
            .1
            .as_ref()
            .expect("LSIG is required to construct SIG");
        let esz = dependencies
            .2
            .as_ref()
            .expect("ESZ is required to construct SIG");

        // Build one CrossSection per (MT, LSIG) pair in parallel, then collect.
        let entries: Vec<(usize, CrossSection)> = mtr
            .par_iter()
            .zip(lsig.par_iter())
            .map(|(mt, start_pos)| {
                let sig_start = *start_pos;

                // Index of the first energy point in ESZ.energy for this XS
                let energy_start_index: usize = data[sig_start - 1].to_bits() as usize;
                // Number of XS values for this reaction
                let num_xs_values: usize = data[sig_start].to_bits() as usize;

                // Cross-section values
                let xs_slice = &data[sig_start + 1..sig_start + 1 + num_xs_values];
                let xs_val = xs_slice.to_vec();

                // Corresponding energy values
                let energy_start = energy_start_index - 1;
                let energy_slice =
                    &esz.energy[energy_start..energy_start + num_xs_values];
                let energy = energy_slice.to_vec();

                (
                    *mt,
                    CrossSection {
                        mt: *mt,
                        energy,
                        xs_val,
                    },
                )
            })
            .collect();

        let xs: CrossSectionMap = entries.into_iter().collect();
        Self(xs)
    }
}

impl std::fmt::Display for SIG {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut sorted_xs: Vec<CrossSection> = self.values().cloned().collect();
        sorted_xs.sort_by_key(|xs| xs.mt);
        let xs_string = sorted_xs
            .iter()
            .map(|xs| format!("{}", xs))
            .collect::<Vec<String>>()
            .join(", ");
        write!(f, "SIG({})", xs_string)
    }
}

//=====================================================================
// Helper struct to represent a cross section.
//=====================================================================
type CrossSectionMap = HashMap<usize, CrossSection>;

#[derive(Debug, Clone)]
pub struct CrossSection {
    pub mt: usize,
    pub energy: Vec<f64>,
    pub xs_val: Vec<f64>,
}

impl std::fmt::Display for CrossSection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "CrossSection(MT={} {})",
            self.mt,
            reaction_type_from_MT(self.mt)
        )
    }
}

#[cfg(not(feature = "local"))]
#[cfg(test)]
mod tests {
    use crate::utils::get_parsed_test_file;

    #[tokio::test]
    async fn test_sig_parsing() {
        let parsed_ace = get_parsed_test_file().await;

        // Check contents
        let sig = parsed_ace.data_blocks.SIG.unwrap();
        assert!(sig.contains_key(&18));

        let fission_xs = sig.get(&18).unwrap();
        assert_eq!(fission_xs.energy.len(), 3);
        assert_eq!(fission_xs.xs_val.len(), fission_xs.energy.len());
        assert_eq!(fission_xs.energy, vec![1.0, 2.0, 3.0]);
        assert_eq!(fission_xs.xs_val, vec![17.0, 38.0, 100.0]);
    }
}
