#![allow(unused)]

//=====================================================================
// Utility functions to aid in testing
//=====================================================================

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, Seek, Write};
use std::path::Path;
use std::time::Instant;
use tempfile::NamedTempFile;
use tokio::sync::OnceCell;

use crate::pace_data::PaceData;
use crate::utils::binary_format::convert_ACE_to_PACE;

// Helper macro: in tests, time a parse and print; in other builds, just
// evaluate the expression without any timing or println noise.
#[cfg(test)]
#[macro_export]
macro_rules! time_it {
    ($label:expr, $expr:expr) => {{
        let start = Instant::now();
        let result = $expr;
        println!("⚛️  {}  ⚛️ : {} μs", $label, start.elapsed().as_micros());
        result
    }};
}

#[cfg(not(test))]
#[macro_export]
macro_rules! time_it {
    ($label:expr, $expr:expr) => {{ $expr }};
}

// These constants and cells hold test file paths and parsed `PaceData`
// so that they are available to all tests and parsed only once.

static TEST_PACE_DATA: OnceCell<PaceData> = OnceCell::const_new();

#[cfg(not(feature = "local"))]
pub const TEST_ACE: &str = "test_nuclear_data_files/test_ascii_ace";

// For local testing (enabled with the `local` feature)
#[cfg(feature = "local")]
pub const TEST_ACE: &str = "local_test_files/uranium_test_file";

// Checks if a file is ASCII by reading the first 1 kB of the file
pub fn is_ascii_file<P: AsRef<Path>>(path: P) -> Result<bool> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut buffer = vec![0; 1024];

    match reader.read(&mut buffer)? {
        0 => Ok(true),
        n => Ok(!buffer[..n]
            .iter()
            .any(|&byte| byte >= 128 || (byte < 32 && !matches!(byte, 9 | 10 | 13)))),
    }
}

// This function simply removes comments from the specially-constructed ASCII ACE test file
fn uncomment_ace_test_file<P: AsRef<Path>>(path: P) -> Result<NamedTempFile> {
    // Open input
    let commented_file = File::open(path)?;
    let reader = BufReader::new(commented_file);

    // Create temp file in the OS temp dir
    let mut tmp = NamedTempFile::new()?;

    // Rewrite without comment lines
    for line_result in reader.lines() {
        let line = line_result?;
        if !line.starts_with("//") {
            writeln!(tmp, "{}", line)?;
        }
    }

    // Ensure contents are on disk before handing out the handle/path
    tmp.flush()?;
    Ok(tmp)
}

// The following code parses example ACE files and caches the resulting
// `PaceData` so it is globally accessible for testing.

// Default test data based on the canonical hydrogen test file.
pub async fn get_parsed_test_file() -> PaceData {
    // In effect, this acts as a sloppy integration test as it involves
    // the parsing of an actual ACE file.
    let data = TEST_PACE_DATA
        .get_or_init(|| async {
            // Path to ACE test file
            let ace_path = Path::new(TEST_ACE);

            // Deal with comments if we are using a custom commented ACE ASCII file, store in temp file
            let tmp = uncomment_ace_test_file(ace_path).expect("Failed to uncomment ACE test file");

            // Convert the ACE file to a PACE file in the temp directory
            let tmp_pace_path = time_it!(
                format!("Time to convert ACE test file {} to PACE", TEST_ACE),
                convert_ACE_to_PACE(tmp.path()).expect("Failed to convert ACE to PACE.")
            );

            // Also copy the generated PACE file into the same directory as the
            // original ACE test file so it can be inspected later if desired.
            if let Some(dir) = ace_path.parent() {
                if let Some(fname) = tmp_pace_path.file_name() {
                    let saved_pace_path = dir.join(fname);
                    // Ignore copy errors here; the temp file is still usable
                    let _ = std::fs::copy(&tmp_pace_path, &saved_pace_path);
                }
            }

            // Parse into PaceData
            let parsed_ace = time_it!(
                format!("Time to parse test PACE file {}", tmp_pace_path.display()),
                PaceData::from_file(&tmp_pace_path).await.unwrap()
            );
            parsed_ace
        })
        .await;

    data.clone()
}
