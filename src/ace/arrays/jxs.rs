#![allow(dead_code)]

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use crate::ace::arrays::NxsArray;
use crate::ace::utils;

// Represents an entry within the JXS array containing location and length of a data block
#[derive(Debug, Clone, PartialEq)]
pub struct JxsEntry {
    pub loc: usize,  // Starting location of the data block
    pub len: usize,  // Length of the data block
}

impl JxsEntry {
    /// Creates a new JxsEntry from location and length values
    fn new(loc: usize, len: usize) -> Self {
        Self { loc, len }
    }

    /// Creates an Option<JxsEntry> from a pair of values, returning None if loc is 0
    fn from_pair(loc: usize, len: usize) -> Option<Self> {
        if loc == 0 {
            None
        } else {
            Some(Self::new(loc, len))
        }
    }
}

// Indices for different values within JXS array, corresponding to the ACE format specification.
#[repr(usize)]
#[derive(Debug, Clone, Copy)]
pub enum JxsIndex {
    Esz = 0,    // Energy table
    Nu = 1,     // Fission nu data
    Mtr = 2,    // MT array
    Lqr = 3,    // Q-value array
    Tyr = 4,    // Reaction type array
    Lsig = 5,   // Table of cross section locators
    Sig = 6,    // Cross sections
    Land = 7,   // Table of angular distribution locators
    And = 8,    // Angular distributions
    Ldlw = 9,   // Table of energy distribution locators
    Dlw = 10,   // Energy distributions
    Gpd = 11,   // Photon production data
    Mtrp = 12,  // Photon production MT array
    Lsigp = 13, // Table of photon production cross section locators
    Sigp = 14,  // Photon production cross sections
    Landp = 15, // Table of photon production angular distribution locators
    Andp = 16,  // Photon production angular distributions
    Ldlwp = 17, // Table of photon production energy distribution locators
    Dlwp = 18,  // Photon production energy distributions
    Yp = 19,    // Table of yield multipliers
    Fis = 20,   // Total fission cross section
    End = 21,   // Last word of the conventional table
    Lund = 22,  // Probability tables
    Dnu = 23,   // Delayed nu-bar data
    Bdd = 24,   // Basic delayed neutron precursor data
    Dnedl = 25, // Table of delayed neutron energy distribution locators
    Dned = 26,  // Delayed neutron energy distributions
    Ptype = 29, // Particle type array
    Ntro = 30,  // Array containing number of particle production reactions
    Next = 31,  // Table of particle production locators
}

// Represents the complete JXS array from an ACE file
#[derive(Clone, Debug, PartialEq)]
pub struct JxsArray {
    esz: Option<JxsEntry>,     // Energy table
    nu: Option<JxsEntry>,      // Fission nu data
    mtr: Option<JxsEntry>,     // MT array
    lqr: Option<JxsEntry>,     // Q-value array
    tyr: Option<JxsEntry>,     // Reaction type array
    lsig: Option<JxsEntry>,    // Table of cross section locators
    sig: Option<JxsEntry>,     // Cross sections
    land: Option<JxsEntry>,    // Table of angular distribution locators
    and: Option<JxsEntry>,     // Angular distributions
    ldlw: Option<JxsEntry>,   // Table of energy distribution locators
    dlw: Option<JxsEntry>,    // Energy distributions
    gpd: Option<JxsEntry>,    // Photon production data
    mtrp: Option<JxsEntry>,   // Photon production MT array
    lsigp: Option<JxsEntry>,  // Table of photon production cross section locators
    sigp: Option<JxsEntry>,   // Photon production cross sections
    landp: Option<JxsEntry>,  // Table of photon production angular distribution locators
    andp: Option<JxsEntry>,   // Photon production angular distributions
    ldlwp: Option<JxsEntry>,  // Table of photon production energy distribution locators
    dlwp: Option<JxsEntry>,   // Photon production energy distributions
    yp: Option<JxsEntry>,     // Table of yield multipliers
    fis: Option<JxsEntry>,    // Total fission cross section
    end: Option<JxsEntry>,    // Last word of the conventional table
    lund: Option<JxsEntry>,   // Probability tables
    dnu: Option<JxsEntry>,    // Delayed nu-bar data
    bdd: Option<JxsEntry>,    // Basic delayed neutron precursor data
    dnedl: Option<JxsEntry>,  // Table of delayed neutron energy distribution locators
    dned: Option<JxsEntry>,   // Delayed neutron energy distributions
    ptype: Option<JxsEntry>,  // Particle type array
    ntro: Option<JxsEntry>,   // Array containing number of particle production reactions
    next: Option<JxsEntry>,   // Table of particle production locators
}

impl JxsArray {
    // Creates a new JxsArray from an ASCII file reader and NXS array information.
    pub fn from_ascii_file(reader: &mut BufReader<File>, nxs_array: &NxsArray) -> Result<Self, Box<dyn Error>> {
        // A JXS array consists of 4 lines, each with eight integers.
        let jxs_array_text = utils::read_lines(reader, 4)?;

        // Parse to integers
        let parsed_jxs_array: Vec<usize> = jxs_array_text
            .iter()
            .flat_map(
                |s| {
                    s.split_whitespace() 
                        .map(|num| num.parse::<usize>())
                        .filter_map(Result::ok)
                    }
                )
            .collect();

        // Builder scheme to ease intialization
        let mut jxs_builder = JxsArrayBuilder::default();
        for i in 0..32 {
            let loc = parsed_jxs_array[i];
            match loc == 0 {
                // Block does not exist
                true => {
                    jxs_builder.set_field(i, None)
                },
                // Block exists
                false => {
                    // Loop forward to find the length of the block
                    let mut next_i = i + 1;
                    while next_i < 32 && parsed_jxs_array[next_i] == 0 {
                        next_i += 1;
                    }
                    let len = if next_i != 32 {
                        parsed_jxs_array[next_i] - loc
                    } else {
                        nxs_array.xxs_len - loc
                    };
                    jxs_builder.set_field(i, Some(JxsEntry::new(loc, len)));
                },
            }
        }
        Ok(jxs_builder.build())
    }
}

#[derive(Default)]
struct JxsArrayBuilder {
    fields: [Option<JxsEntry>; 32],
}

impl JxsArrayBuilder {
    fn set_field(&mut self, index: usize, value: Option<JxsEntry>) {
        self.fields[index] = value;
    }

    fn build(self) -> JxsArray {
        JxsArray {
            esz: self.fields[JxsIndex::Esz as usize].clone(),       // Energy table
            nu: self.fields[JxsIndex::Nu as usize].clone(),         // Fission nu data
            mtr: self.fields[JxsIndex::Mtr as usize].clone(),       // MT array
            lqr: self.fields[JxsIndex::Lqr as usize].clone(),       // Q-value array
            tyr: self.fields[JxsIndex::Tyr as usize].clone(),       // Reaction type array
            lsig: self.fields[JxsIndex::Lsig as usize].clone(),     // Table of cross section locators
            sig: self.fields[JxsIndex::Sig as usize].clone(),       // Cross sections
            land: self.fields[JxsIndex::Land as usize].clone(),     // Table of angular distribution locators
            and: self.fields[JxsIndex::And as usize].clone(),       // Angular distributions
            ldlw: self.fields[JxsIndex::Ldlw as usize].clone(),     // Table of energy distribution locators
            dlw: self.fields[JxsIndex::Dlw as usize].clone(),       // Energy distributions
            gpd: self.fields[JxsIndex::Gpd as usize].clone(),       // Photon production data
            mtrp: self.fields[JxsIndex::Mtrp as usize].clone(),     // Photon production MT array
            lsigp: self.fields[JxsIndex::Lsigp as usize].clone(),   // Table of photon production cross section locators
            sigp: self.fields[JxsIndex::Sigp as usize].clone(),     // Photon production cross sections
            landp: self.fields[JxsIndex::Landp as usize].clone(),   // Table of photon production angular distribution locators
            andp: self.fields[JxsIndex::Andp as usize].clone(),     // Photon production angular distributions
            ldlwp: self.fields[JxsIndex::Ldlwp as usize].clone(),   // Table of photon production energy distribution locators
            dlwp: self.fields[JxsIndex::Dlwp as usize].clone(),     // Photon production energy distributions
            yp: self.fields[JxsIndex::Yp as usize].clone(),         // Table of yield multipliers
            fis: self.fields[JxsIndex::Fis as usize].clone(),       // Total fission cross section
            end: self.fields[JxsIndex::End as usize].clone(),       // Last word of the conventional table
            lund: self.fields[JxsIndex::Lund as usize].clone(),     // Probability tables
            dnu: self.fields[JxsIndex::Dnu as usize].clone(),       // Delayed nu-bar data
            bdd: self.fields[JxsIndex::Bdd as usize].clone(),       // Basic delayed neutron precursor data
            dnedl: self.fields[JxsIndex::Dnedl as usize].clone(),   // Table of delayed neutron energy distribution locators
            dned: self.fields[JxsIndex::Dned as usize].clone(),     // Delayed neutron energy distributions
            ptype: self.fields[JxsIndex::Ptype as usize].clone(),   // Particle type array
            ntro: self.fields[JxsIndex::Ntro as usize].clone(),     // Array containing number of particle production reactions
            next: self.fields[JxsIndex::Next as usize].clone(),     // Table of particle production locators
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jxs_entry_new() {
        let entry = JxsEntry::new(10, 20);
        assert_eq!(entry.loc, 10);
        assert_eq!(entry.len, 20);
    }

    #[test]
    fn test_jxs_entry_from_pair() {
        let entry = JxsEntry::from_pair(10, 20);
        assert_eq!(entry, Some(JxsEntry::new(10, 20)));

        let none_entry = JxsEntry::from_pair(0, 20);
        assert_eq!(none_entry, None);
    }

    #[test]
    fn test_jxs_parsing() {
        // Simulate ACE JXS array
        let jxs_text = concat!(
            "    1    0    3    4    5    6    7    8\n",
            "    9   10    0    0    0   14   15   16\n",
            "   17   18   19   20   21   22   23   24\n",
            "   25   26   27   28   29   30   31   32\n"
        );
        let mut reader = utils::create_reader_from_string(jxs_text);

        // Simulate NXS array
        let nxs = NxsArray {
            xxs_len: 100,
            za: 5010,
            nes: 941,
            ntr: 55,
            nr: 35,
            ntrp: 38,
            ntype: 2,
            npcr: 0,
            s: 0,
            z: 5,
            a: 10,
        };

        // Parse the array
        let jxs = JxsArray::from_ascii_file(&mut reader, &nxs).expect("Failed to parse JXS array");

        // Check fields
        assert_eq!(jxs.esz, Some(JxsEntry::new(1, 2)));
        assert_eq!(jxs.nu, None);
        assert_eq!(jxs.mtr, Some(JxsEntry::new(3, 1)));
        assert_eq!(jxs.lqr, Some(JxsEntry::new(4, 1)));
        assert_eq!(jxs.tyr, Some(JxsEntry::new(5, 1)));
        assert_eq!(jxs.lsig, Some(JxsEntry::new(6, 1)));
        assert_eq!(jxs.sig, Some(JxsEntry::new(7, 1)));
        assert_eq!(jxs.land, Some(JxsEntry::new(8, 1)));
        assert_eq!(jxs.and, Some(JxsEntry::new(9, 1)));
        assert_eq!(jxs.ldlw, Some(JxsEntry::new(10, 4)));
        assert_eq!(jxs.dlw, None);
        assert_eq!(jxs.gpd, None);
        assert_eq!(jxs.mtrp, None);
        assert_eq!(jxs.lsigp, Some(JxsEntry::new(14, 1)));
        assert_eq!(jxs.sigp, Some(JxsEntry::new(15, 1)));
        assert_eq!(jxs.landp, Some(JxsEntry::new(16, 1)));
        assert_eq!(jxs.andp, Some(JxsEntry::new(17, 1)));
        assert_eq!(jxs.ldlwp, Some(JxsEntry::new(18, 1)));
        assert_eq!(jxs.dlwp, Some(JxsEntry::new(19, 1)));
        assert_eq!(jxs.yp, Some(JxsEntry::new(20, 1)));
        assert_eq!(jxs.fis, Some(JxsEntry::new(21, 1)));
        assert_eq!(jxs.end, Some(JxsEntry::new(22, 1)));
        assert_eq!(jxs.lund, Some(JxsEntry::new(23, 1)));
        assert_eq!(jxs.dnu, Some(JxsEntry::new(24, 1)));
        assert_eq!(jxs.bdd, Some(JxsEntry::new(25, 1)));
        assert_eq!(jxs.dnedl, Some(JxsEntry::new(26, 1)));
        assert_eq!(jxs.dned, Some(JxsEntry::new(27, 1)));
        assert_eq!(jxs.ptype, Some(JxsEntry::new(30, 1)));
        assert_eq!(jxs.ntro, Some(JxsEntry::new(31, 1)));
        assert_eq!(jxs.next, Some(JxsEntry::new(32, 68)));
    }
}