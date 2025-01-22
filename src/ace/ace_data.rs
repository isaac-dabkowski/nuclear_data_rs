use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;

use crate::ace::header::AceHeader;
use crate::ace::arrays::{IzawPair, IzawArray, JxsArray, NxsArray};
use crate::ace::data_blocks::DataBlocks;
use crate::ace::utils::is_ascii_file;

use super::data_blocks;

#[derive(Clone)]
pub struct AceIsotopeData {
    pub header: AceHeader,
    pub izaw_array: IzawArray,
    pub nxs_array: NxsArray,
    pub jxs_array: JxsArray,
    pub data_blocks: DataBlocks
}

impl AceIsotopeData {
    pub fn from_file<P: AsRef<Path>>(file_path: P) -> Result<Self, Box<dyn Error>> {
        let path = file_path.as_ref();

        // Invoke ASCII or binary parsing based on file type
        if is_ascii_file(path)? {
            // Parse ASCII file
            let ace_data = AceIsotopeData::from_ascii_file(path)?;
            Ok(ace_data)
        } else {
            // Parse binary file
            todo!()
        }
    }

    // Create an AceIsotopeData object from an ASCII file
    pub fn from_ascii_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let file = File::open(path).map_err(|e| format!("Error opening ACE ASCII file: {}", e))?;
        let mut reader = BufReader::new(file);

        // Process the header
        let header = AceHeader::from_ascii_file(&mut reader)?;

        // Process the IZAW array
        let izaw_array = IzawArray::from_ascii_file(&mut reader)?;

        // Process the NXS array
        let nxs_array = NxsArray::from_ascii_file(&mut reader)?;

        // Process the JXS array
        let jxs_array = JxsArray::from_ascii_file(&mut reader, &nxs_array)?;

        // Process the XXS array into each block's raw text
        let data_blocks = DataBlocks::from_ascii_file(&mut reader, &nxs_array, &jxs_array)?;

        Ok(Self { header, izaw_array, nxs_array, jxs_array, data_blocks })
    }

    // ZAID of the isotope
    #[inline]
    pub fn zaid(&self) -> String {
        self.header.zaid.clone()
    }

    // SZAID of the isotope (version 2.0.0 and later)
    #[inline]
    pub fn szaid(&self) -> Option<String> {
        self.header.szaid.clone()
    }

    // Atomic mass fraction
    #[inline]
    pub fn atomic_mass_fraction(&self) -> f64 {
        self.header.atomic_mass_fraction
    }

    // kT
    #[inline]
    pub fn kT(&self) -> f64 {
        self.header.kT
    }

    // Temperature in Kelvin
    #[inline]
    pub fn temperature(&self) -> f64 {
        self.header.temperature
    }

    // S alpha beta pairs of ZAIDs and atomic weight ratios
    #[inline]
    pub fn s_a_b_pairs(&self) -> Vec<IzawPair> {
        self.izaw_array.pairs.clone()
    }

    // Number of entries in the main data array
    #[inline]
    pub fn num_entries(&self) -> usize {
        self.nxs_array.xxs_len
    }

    // Number of energies
    #[inline]
    pub fn num_energies(&self) -> usize {
        self.nxs_array.nes
    }

    // ZA of the isotope
    #[inline]
    pub fn za(&self) -> usize {
        self.nxs_array.za
    }

    // Atomic number
    #[inline]
    pub fn z(&self) -> usize {
        self.nxs_array.z
    }

    // Mass number
    #[inline]
    pub fn a(&self) -> usize {
        self.nxs_array.a
    }

    // Energy grid from ESZ block
    #[inline]
    pub fn energies(&self) -> Vec<f64> {
        self.data_blocks.ESZ.as_ref().unwrap().energy.clone()
    }
}

#[cfg(test)]
mod ascii_tests {
    use crate::ace::utils::get_parsed_ascii_for_testing;

    #[test]
    fn test_szaid_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.szaid(), Some(String::from("1001.800nc")));
    }

    #[test]
    fn test_zaid_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.zaid(), String::from("1001.00c"));
    }

    #[test]
    fn test_atomic_mass_fraction_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.atomic_mass_fraction(), 0.999167);
    }

    #[test]
    fn test_kT_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.kT(), 2.5301e-08);
    }

    #[test]
    fn test_temperature_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.temperature(), 293.6059129982851);
    }

    #[test]
    fn test_izaw_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        for za_iz_pair in parsed_ace.s_a_b_pairs() {
            assert_eq!(za_iz_pair.za, 0);
            assert_eq!(za_iz_pair.iz, 0.0);
        }
        assert_eq!(parsed_ace.s_a_b_pairs().len(), 16)
    }

    #[test]
    fn test_num_entries_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.num_entries(), 10257);
    }

    #[test]
    fn test_num_energies_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.num_energies(), 631);
    }

    #[test]
    fn test_za_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.za(), 1001);
    }

    #[test]
    fn test_z_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.z(), 1);
    }

    #[test]
    fn test_a_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.a(), 1);
    }

    #[test]
    fn test_esz_energies_parsing() {
        let parsed_ace = get_parsed_ascii_for_testing();
        assert_eq!(parsed_ace.energies().len(), parsed_ace.num_energies());
    }
}

#[cfg(test)]
mod binary_tests {
}