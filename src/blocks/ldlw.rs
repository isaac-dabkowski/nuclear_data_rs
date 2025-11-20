use std::ops::Deref;

use crate::arrays::Arrays;
use crate::blocks::BlockType;
use crate::blocks::block_traits::{Process, PullFromXXS, block_range_to_slice, get_block_start};

//=====================================================================
// LDLW data block
//
// Contains locations of energy distributions of secondary neutrons.
//=====================================================================
#[derive(Debug, Clone, PartialEq)]
pub struct LDLW(pub Vec<usize>);

impl Deref for LDLW {
    type Target = Vec<usize>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<'a> PullFromXXS<'a> for LDLW {
    fn pull_from_xxs_array(arrays: &'a Arrays) -> Option<&'a [f64]> {
        // We expect LDLW if NXS(5) (NR) != 0.
        let has_secondary_neutron_reactions = arrays.nxs.nr != 0;

        // Get the starting index of the block in the XXS array
        let block_start = get_block_start(
            &BlockType::LDLW,
            arrays,
            has_secondary_neutron_reactions,
            "LDLW is expected if NXS(5) (NR) != 0, but LDLW was not found.".to_string(),
        )?;

        // Calculate the block length, see the LDLW description in the ACE spec
        let num_reactions_with_secondary_neutrons = arrays.nxs.nr;
        let block_length = num_reactions_with_secondary_neutrons;

        // Return the block's raw data as a slice
        Some(block_range_to_slice(block_start, block_length, arrays))
    }
}

impl<'a> Process<'a> for LDLW {
    type Dependencies = ();

    fn process(data: &[f64], _arrays: &Arrays, _dependencies: ()) -> Self {
        Self(data.iter().map(|val| val.to_bits() as usize).collect())
    }
}

impl std::fmt::Display for LDLW {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "LDLW({} reactions)", self.len())
    }
}

#[cfg(not(feature = "local"))]
#[cfg(test)]
mod tests {
    use crate::utils::get_parsed_test_file;

    #[tokio::test]
    async fn test_ldlw_parsing() {
        let parsed_ace = get_parsed_test_file().await;

        // Check contents
        let ldlw = parsed_ace.data_blocks.LDLW.unwrap();
        assert_eq!(ldlw.len(), 1);
        assert_eq!(ldlw[0], 1);
    }
}
