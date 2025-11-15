use criterion::{criterion_group, criterion_main};

mod array_bench;
mod bench_helpers;
mod block_bench;
mod header_bench;
mod pace_data_bench;

criterion_group!(header_bench, header_bench::bench_header);
criterion_group!(
    array_bench,
    array_bench::bench_izaw,
    array_bench::bench_nxs,
    array_bench::bench_jxs
);
criterion_group!(
    block_bench,
    block_bench::bench_esz,
    block_bench::bench_mtr,
    block_bench::bench_lqr,
    block_bench::bench_lsig,
    block_bench::bench_sig,
    block_bench::bench_tyr,
    block_bench::bench_nu,
    block_bench::bench_dnu,
    block_bench::bench_bdd,
    block_bench::bench_land,
    block_bench::bench_and
);
criterion_group!(pace_data, pace_data_bench::bench_pace_data);
criterion_main!(header_bench, array_bench, block_bench, pace_data);
