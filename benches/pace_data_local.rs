// #![cfg(feature = "local")]

// use criterion::{Criterion, black_box, criterion_group, criterion_main};
// use tokio::runtime::Runtime;

// use pace_rs::PaceData;
// use pace_rs::utils::binary_format::convert_ACE_to_PACE;

// const LOCAL_TEST_PACE_PATH: &str = "local_test_files/92235.800nc.pace";
// const LOCAL_TEST_ACE_PATH: &str = "local_test_files/uranium_test_file";

// fn bench_parse_local_test_file(c: &mut Criterion) {
//     let runtime = Runtime::new().expect("failed to construct Tokio runtime");

//     c.bench_function("pace_parse::local_data", |b| {
//         b.iter(|| {
//             let parsed = runtime
//                 .block_on(PaceData::from_file(LOCAL_TEST_PACE_PATH))
//                 .expect("parse should succeed");
//             black_box(parsed);
//         })
//     });
// }

// criterion_group!(pace_data_local_benches, bench_parse_local_test_file);
// criterion_main!(pace_data_local_benches);
