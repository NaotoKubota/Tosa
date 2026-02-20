use std::env;
use std::fs;
use std::path::Path;

fn main() {
    // Write the VERSION file from Cargo.toml's version
    let version = env::var("CARGO_PKG_VERSION").unwrap();
    let manifest_dir = env::var("CARGO_MANIFEST_DIR").unwrap();
    let version_path = Path::new(&manifest_dir).join("VERSION");
    fs::write(version_path, format!("{}\n", version)).expect("Failed to write VERSION file");
}
