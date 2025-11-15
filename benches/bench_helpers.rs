#![allow(dead_code)]

use std::sync::OnceLock;
use std::time::Duration;

use criterion::Criterion;
use tokio::runtime::Runtime;

use pace_rs::arrays::{JxsArray, NxsArray};
use pace_rs::blocks::DataBlocks;
use pace_rs::utils::PaceMmap;
use pace_rs::{Arrays, PaceData};

pub fn default_benchmark_config() -> Criterion {
    Criterion::default()
        .measurement_time(Duration::from_secs_f64(0.1))
        .warm_up_time(Duration::from_secs_f64(0.1))
        .configure_from_args()
}

//=====================================================================
// Bench-only helpers that keep values parsed from disk so other
// benches can reuse it. This shares a single Tokio runtime and lazily
// parses the data when first requested.
//=====================================================================

#[cfg(feature = "local")]
const TEST_PACE_PATH: &str = "local_test_files/92235.800nc.pace";

#[cfg(not(feature = "local"))]
const TEST_PACE_PATH: &str = "test_nuclear_data_files/1100.800nc.pace";

// Returns the cached `PaceData` parsed from `TEST_PACE_PATH`.
pub fn cached_pace_data() -> &'static PaceData {
    static PACE_DATA: OnceLock<PaceData> = OnceLock::new();
    static RUNTIME: OnceLock<Runtime> = OnceLock::new();

    let runtime = RUNTIME
        .get_or_init(|| Runtime::new().expect("Failed to create tokio runtime for benchmarks."));

    PACE_DATA.get_or_init(|| {
        runtime
            .block_on(PaceData::from_file(TEST_PACE_PATH))
            .expect("failed to parse pace data for benchmarks")
    })
}

// Returns the cached `PaceMmap` parsed from `TEST_PACE_PATH`.
pub fn cached_mmap() -> &'static PaceMmap {
    static MMAP: OnceLock<PaceMmap> = OnceLock::new();
    static RUNTIME: OnceLock<Runtime> = OnceLock::new();

    let runtime = RUNTIME
        .get_or_init(|| Runtime::new().expect("Failed to create tokio runtime for benchmarks."));

    MMAP.get_or_init(|| {
        runtime
            .block_on(async { PaceMmap::from_file(TEST_PACE_PATH) })
            .expect("Failed to get mem map.")
    })
}

// Returns the cached `DataBlocks` parsed from `TEST_PACE_PATH`.
pub fn cached_blocks() -> &'static DataBlocks {
    static BLOCKS: OnceLock<DataBlocks> = OnceLock::new();
    static RUNTIME: OnceLock<Runtime> = OnceLock::new();

    let runtime = RUNTIME
        .get_or_init(|| Runtime::new().expect("Failed to create tokio runtime for benchmarks."));

    let mmap = cached_mmap();
    let nxs_array = NxsArray::from_PACE(&mmap).expect("Failed to parse NXS.");
    let jxs_array = JxsArray::from_PACE(&mmap).expect("Failed to parse JXS.");

    BLOCKS.get_or_init(|| {
        runtime
            .block_on(async { DataBlocks::from_PACE(&mmap, &nxs_array, &jxs_array) })
            .expect("Failed to get DataBlocks.")
    })
}

// Returns the cached `Arrays` parsed from `TEST_PACE_PATH`.
pub fn cached_arrays() -> &'static Arrays<'static> {
    // Cache the parsed NXS and JXS arrays once, then build a single
    // `Arrays<'static>` that borrows them and the XXS slice.
    static NXS_ARRAY: OnceLock<NxsArray> = OnceLock::new();
    static JXS_ARRAY: OnceLock<JxsArray> = OnceLock::new();
    static ARRAYS: OnceLock<Arrays<'static>> = OnceLock::new();

    let mmap: &'static PaceMmap = cached_mmap();

    let nxs: &'static NxsArray =
        NXS_ARRAY.get_or_init(|| NxsArray::from_PACE(mmap).expect("Failed to parse NXS."));

    let jxs: &'static JxsArray =
        JXS_ARRAY.get_or_init(|| JxsArray::from_PACE(mmap).expect("Failed to parse JXS."));

    ARRAYS.get_or_init(|| Arrays {
        nxs,
        jxs,
        xxs: mmap.xxs_array(),
    })
}
