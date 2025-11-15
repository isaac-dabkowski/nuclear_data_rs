#[path = "bench_helpers.rs"]
mod bench_helpers;

use criterion::{Criterion, black_box};

use pace_rs::arrays::{IzawArray, JxsArray, NxsArray};

//=====================================================================
// Benchmark the array (IZAW/NXS/JXS/XXS) parsing.
//=====================================================================

pub fn bench_izaw(_: &mut Criterion) {
    let mut c = bench_helpers::default_benchmark_config();

    let mmap = bench_helpers::cached_mmap();

    c.bench_function("array_bench:izaw", |b| {
        b.iter(|| {
            let parsed = IzawArray::from_PACE(&mmap).expect("Failed to parse IZAW.");
            black_box(parsed);
        })
    });
}

pub fn bench_nxs(_: &mut Criterion) {
    let mut c = bench_helpers::default_benchmark_config();

    let mmap = bench_helpers::cached_mmap();

    c.bench_function("array_bench:nxs", |b| {
        b.iter(|| {
            let parsed = NxsArray::from_PACE(&mmap).expect("Failed to parse NXS.");
            black_box(parsed);
        })
    });
}

pub fn bench_jxs(_: &mut Criterion) {
    let mut c = bench_helpers::default_benchmark_config();

    let mmap = bench_helpers::cached_mmap();

    c.bench_function("array_bench:jxs", |b| {
        b.iter(|| {
            let parsed = JxsArray::from_PACE(&mmap).expect("Failed to parse JXS.");
            black_box(parsed);
        })
    });
}
