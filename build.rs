// Build script "inspired" by
// https://github.com/imartayan/CBL/blob/328bcc6932a2169c1d7a03b047c685d18f21643a/build.rs#L9

/// This is a hack to support "dynamic" NB_HASH values.
/// NB_HASH is implemented as a const generic in our code
/// as we expect it to remain constant across executions
/// and benefit from compile-time optimizations.
/// This build script will set the value of NB_HASH at compile-time
/// from an environment variable, so one can easily build
/// the project with the desired NB_HASH value.
/// This will not re-build if the NB_HASH value does not change.
fn build_constants() {
    let out_dir: std::path::PathBuf = std::env::var("OUT_DIR")
        .expect("Failed to obtain OUT_DIR")
        .into();
    let mut code = Vec::new();

    println!("cargo:rerun-if-env-changed=NB_HASH");
    let nb_hash: usize = std::env::var("NB_HASH")
        .unwrap_or_else(|_| "1".into())
        .parse()
        .expect("Failed to parse NB_HASH");
    assert!(nb_hash >= 1, "NB_HASH must be ≥ 1");
    assert!(nb_hash <= 16, "NB_HASH must be ≤ 16");
    code.push(format!("pub const NB_HASH: usize = {nb_hash};"));

    std::fs::write(out_dir.join("constants.rs"), code.join("\n"))
        .expect("Failed to write const file");
}

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    build_constants();
}
