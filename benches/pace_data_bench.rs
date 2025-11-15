#[path = "bench_helpers.rs"]
mod bench_helpers;

use criterion::{Criterion, black_box};
use tokio::runtime::Runtime;

use pace_rs::PaceData;

const TEST_PACE_PATH: &str = "test_nuclear_data_files/1100.800nc.pace";

pub fn bench_pace_data(_: &mut Criterion) {
    let mut c = bench_helpers::default_benchmark_config();

    let runtime = Runtime::new().expect("Failed to construct Tokio runtime");

    c.bench_function("pace_parse::test_data", |b| {
        b.iter(|| {
            let parsed = runtime
                .block_on(PaceData::from_file(TEST_PACE_PATH))
                .expect("Parse should succeed");
            black_box(parsed);
        })
    });
}
