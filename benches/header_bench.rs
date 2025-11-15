#[path = "bench_helpers.rs"]
mod bench_helpers;

use criterion::{Criterion, black_box, criterion_group, criterion_main};

use pace_rs::header::Header;

//=====================================================================
// Benchmark the header parsing.
//=====================================================================

pub fn bench_header(_: &mut Criterion) {
    let mut c = bench_helpers::default_benchmark_config();

    let mmap = bench_helpers::cached_mmap();

    c.bench_function("header_bench:header", |b| {
        b.iter(|| {
            let parsed = Header::from_PACE(&mmap).expect("Failed to parse header.");
            black_box(parsed);
        })
    });
}

criterion_group!(header_bench, bench_header);
criterion_main!(header_bench);
