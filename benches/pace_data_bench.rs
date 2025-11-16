#[path = "bench_helpers.rs"]
mod bench_helpers;

use std::time::Duration;

use criterion::{Criterion, black_box};
use tokio::runtime::Runtime;

use pace_rs::PaceData;

pub fn bench_pace_data(_: &mut Criterion) {
    let mut c = Criterion::default()
        .measurement_time(Duration::from_secs_f64(1.0))
        .warm_up_time(Duration::from_secs_f64(0.1))
        .configure_from_args();

    let runtime = Runtime::new().expect("Failed to construct Tokio runtime");

    c.bench_function("pace_parse::test_data", |b| {
        b.iter(|| {
            let parsed = runtime
                .block_on(PaceData::from_file(bench_helpers::TEST_PACE_PATH))
                .expect("Parse should succeed");
            black_box(parsed);
        })
    });
}
