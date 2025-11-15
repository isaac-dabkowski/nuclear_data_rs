#[path = "bench_helpers.rs"]
mod bench_helpers;

use criterion::{Criterion, black_box};

use pace_rs::arrays::Arrays;
use pace_rs::blocks::{
    AND, BDD, DNU, DataBlocks, ESZ, LAND, LQR, LSIG, MTR, NU, SIG, TYR, block_traits::Parse,
};

fn run_block_bench<B, F>(label: &'static str, f: F)
where
    F: Fn(&'static Arrays<'static>, &'static DataBlocks) -> Option<B>,
{
    let mut c = bench_helpers::default_benchmark_config();
    let arrays: &'static Arrays<'static> = bench_helpers::cached_arrays();
    let blocks: &'static DataBlocks = bench_helpers::cached_blocks();

    c.bench_function(label, |b| {
        b.iter(|| {
            let parsed = f(arrays, blocks);
            black_box(parsed);
        })
    });
}

pub fn bench_esz(_: &mut Criterion) {
    run_block_bench("block_bench::ESZ", |arrays, _| ESZ::parse(arrays, ()));
}

pub fn bench_mtr(_: &mut Criterion) {
    run_block_bench("block_bench::MTR", |arrays, _| MTR::parse(arrays, ()));
}

pub fn bench_lsig(_: &mut Criterion) {
    run_block_bench("block_bench::LSIG", |arrays, _| LSIG::parse(arrays, ()));
}

pub fn bench_nu(_: &mut Criterion) {
    run_block_bench("block_bench::NU", |arrays, _| NU::parse(arrays, ()));
}

pub fn bench_dnu(_: &mut Criterion) {
    run_block_bench("block_bench::DNU", |arrays, _| DNU::parse(arrays, ()));
}

pub fn bench_bdd(_: &mut Criterion) {
    run_block_bench("block_bench::BDD", |arrays, _| BDD::parse(arrays, ()));
}

pub fn bench_lqr(_: &mut Criterion) {
    run_block_bench("block_bench::LQR", |arrays, blocks| {
        LQR::parse(arrays, &blocks.MTR)
    });
}

pub fn bench_sig(_: &mut Criterion) {
    run_block_bench("block_bench::SIG", |arrays, blocks| {
        SIG::parse(arrays, (&blocks.MTR, &blocks.LSIG, &blocks.ESZ))
    });
}

pub fn bench_tyr(_: &mut Criterion) {
    run_block_bench("block_bench::TYR", |arrays, blocks| {
        TYR::parse(arrays, &blocks.MTR)
    });
}

pub fn bench_land(_: &mut Criterion) {
    run_block_bench("block_bench::LAND", |arrays, blocks| {
        LAND::parse(arrays, &blocks.MTR)
    });
}

pub fn bench_and(_: &mut Criterion) {
    run_block_bench("block_bench::AND", |arrays, blocks| {
        AND::parse(arrays, (&blocks.TYR, &blocks.LAND))
    });
}
