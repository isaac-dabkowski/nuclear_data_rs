use std::error::Error;

#[cfg(test)]
use std::time::Instant;

use crate::arrays::{Arrays, JxsArray, NxsArray, XxsArray};
use crate::blocks::block_traits::Parse;
use crate::blocks::{AND, BDD, DNU, ESZ, LAND, LQR, LSIG, MTR, NU, SIG, TYR};
use crate::time_it;
use crate::utils::PaceMmap;

#[derive(Clone, Debug, Default)]
pub struct DataBlocks {
    pub ESZ: Option<ESZ>,
    pub MTR: Option<MTR>,
    pub LSIG: Option<LSIG>,
    pub SIG: Option<SIG>,
    pub LQR: Option<LQR>,
    pub NU: Option<NU>,
    pub DNU: Option<DNU>,
    pub BDD: Option<BDD>,
    pub TYR: Option<TYR>,
    pub LAND: Option<LAND>,
    pub AND: Option<AND>,
}

impl DataBlocks {
    pub fn from_PACE(
        mmap: &PaceMmap,
        nxs_array: &NxsArray,
        jxs_array: &JxsArray,
    ) -> Result<Self, Box<dyn Error>> {
        // Recall that this array is returned as f64's, we will parse these values back to
        // integers where appropriate later
        let xxs_array: &XxsArray = mmap.xxs_array();

        // Construct the Arrays struct
        let arrays = Arrays {
            nxs: nxs_array,
            jxs: jxs_array,
            xxs: xxs_array,
        };

        // Process the data blocks from the binary ACE file
        // -------------------------------
        // Blocks which are always present
        // -------------------------------
        // Energy grid
        let esz = time_it!("ESZ", ESZ::parse(&arrays, ()));

        // -------------------------------------------
        // Blocks present if isotope has reactions
        // other than elastic scattering (NXS(4) != 0)
        // -------------------------------------------
        // Reaction MT values
        let mtr = time_it!("MTR", MTR::parse(&arrays, ()));
        // Q values
        let lqr = time_it!("LQR", LQR::parse(&arrays, &mtr));
        // Cross section locations
        let lsig = time_it!("LSIG", LSIG::parse(&arrays, ()));
        // Cross section values
        let sig = time_it!("SIG", SIG::parse(&arrays, (&mtr, &lsig, &esz)));
        // Secondary neutron information
        let tyr = time_it!("TYR", TYR::parse(&arrays, &mtr));

        // -------------------------------------------
        // Blocks present if fission nu data is
        // available (JXS(2) != 0)
        // -------------------------------------------
        // Fission nu values
        let nu = time_it!("NU", NU::parse(&arrays, ()));
        // Fission dnu values
        let dnu = time_it!("DNU", DNU::parse(&arrays, ()));
        // Fission precursor data values
        let bdd = time_it!("BDD", BDD::parse(&arrays, ()));

        // --------------------------------------------------------------------------------
        // Blocks which are always present, but where having MTR makes them easier to parse
        // --------------------------------------------------------------------------------
        // Secondary neutron angular distribution locations
        let land = time_it!("LAND", LAND::parse(&arrays, &mtr));
        // Secondary neutron angular distributions
        let and = time_it!("AND", AND::parse(&arrays, (&tyr, &land)));

        Ok(Self {
            ESZ: esz,
            MTR: mtr,
            LSIG: lsig,
            SIG: sig,
            LQR: lqr,
            DNU: dnu,
            NU: nu,
            BDD: bdd,
            TYR: tyr,
            LAND: land,
            AND: and,
        })
    }
}

impl std::fmt::Display for DataBlocks {
    fn fmt(&self, _: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        todo!()
    }
}
