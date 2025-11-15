use std::sync::OnceLock;

use tokio::runtime::Runtime;

use pace_rs::PaceData;
use pace_rs::utils::PaceMmap;

//=====================================================================
// Bench-only helper that keeps values parsed from disk so other
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

    let runtime = RUNTIME.get_or_init(|| {
        Runtime::new().expect("Failed to create tokio runtime for benchmarks.")
    });

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

    let runtime = RUNTIME.get_or_init(|| {
        Runtime::new().expect("Failed to create tokio runtime for benchmarks.")
    });

    MMAP.get_or_init(|| {
        runtime
            .block_on(async{ PaceMmap::from_file(TEST_PACE_PATH) })
            .expect("Failed to get mem map.")
    })
}

// WRITE FUNCTION THAT CAN CACHE DATA BLOCKS AND EXPOSE THEM ALL, MIGHT GET TRICKY
// pub fn cached_blocks() -> &'static PaceMmap {
//     static MMAP: OnceLock<PaceMmap> = OnceLock::new();
//     static RUNTIME: OnceLock<Runtime> = OnceLock::new();

//     let runtime = RUNTIME.get_or_init(|| {
//         Runtime::new().expect("Failed to create tokio runtime for benchmarks.")
//     });

//     MMAP.get_or_init(|| {
//         runtime
//             .block_on(async{ PaceMmap::from_file(TEST_PACE_PATH) })
//             .expect("Failed to get DataBlocks.")
//     })
// }