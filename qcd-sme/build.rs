use std::env;
use std::io::Write;
use std::path;

#[cfg(debug_assertions)]
const HEADERS_DIR: &str = "debug/include";
#[cfg(not(debug_assertions))]
const HEADERS_DIR: &str = "release/include";
const C_HEADER_FILE: &str = "qcd_sme.h";
const CPP_HEADER_FILE: &str = "qcd_sme.hpp";
const CYTHON_HEADER_FILE: &str = "qcd_sme.pxd";

const INCLUDE_GUARD: &str = "__QCD_SME_H__";

fn main() {
    println!("cargo::warning=There's a known bug in the *_w_field_config functions: you still have to call qcd_sme::consts::set_number_of_colors(nc) for nc != 3 and nf != 0. This bug also makes the functions not thread-safe with respect to a change of nc. I will eventually fix this but I'm too busy right now.");

    let crate_dir = path::PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let mut target_dir = crate_dir.join("target");
    // This crate depends on other crates and it will be built last, by which
    // time hopefully we'll already have a 'target' directory somewhere.
    // Workaround for a proper CARGO_TARGET_DIR not existing.
    if !target_dir.is_dir() {
        target_dir = crate_dir.parent().unwrap().join("target");
    }
    let c_header_path = target_dir.join(HEADERS_DIR).join(C_HEADER_FILE);
    let cpp_header_path = target_dir.join(HEADERS_DIR).join(CPP_HEADER_FILE);
    let cython_header_path = target_dir.join(HEADERS_DIR).join(CYTHON_HEADER_FILE);

    // C header
    cbindgen::Builder::new()
        .with_language(cbindgen::Language::C)
        .with_crate(&crate_dir)
        .with_include_guard(INCLUDE_GUARD)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(&c_header_path);

    let mut c_header = std::fs::read_to_string(&c_header_path).expect("Unable to read bindings");
    let rx = regex::Regex::new(r"typedef\s+Complex<R>\s+C\s*;").unwrap();
    if rx.find(&c_header).is_some() {
        c_header = rx
            .replace(
                &c_header,
                "struct Complex {R re; R im;};\n\ntypedef struct Complex C;",
            )
            .to_string();
        std::fs::File::create(&c_header_path)
            .expect("Unable to modify bindings")
            .write_all(c_header.as_bytes())
            .expect("Unable to modify bindings");
    };
    drop(c_header);

    // C++ header
    cbindgen::Builder::new()
        .with_crate(&crate_dir)
        .with_include_guard(INCLUDE_GUARD)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(&cpp_header_path);

    let mut cpp_header =
        std::fs::read_to_string(&cpp_header_path).expect("Unable to read bindings");
    let rx = regex::Regex::new(r"using\s+C\s+=\s+Complex<R>\s*;").unwrap();
    if rx.find(&cpp_header).is_some() {
        cpp_header = rx
            .replace(
                &cpp_header,
                "struct Complex {R re; R im;};\n\nusing C = Complex;",
            )
            .to_string();
        std::fs::File::create(&cpp_header_path)
            .expect("Unable to modify bindings")
            .write_all(cpp_header.as_bytes())
            .expect("Unable to modify bindings");
    };
    drop(cpp_header);

    // Cython header
    let mut cython_config = cbindgen::Config::default();
    cython_config.language = cbindgen::Language::Cython;
    cython_config.cython.header = Some(String::from("\"") + C_HEADER_FILE + "\"");
    cbindgen::Builder::new()
        .with_config(cython_config)
        .with_crate(&crate_dir)
        .generate()
        .expect("Unable to generate bindings")
        .write_to_file(&cython_header_path);

    let mut cython_header =
        std::fs::read_to_string(&cython_header_path).expect("Unable to read bindings");
    let rx = regex::Regex::new(r"ctypedef\s+Complex<R>\s+C\s*;").unwrap();
    if rx.find(&cython_header).is_some() {
        cython_header = rx
            .replace(
                &cython_header,
                "cdef struct Complex:\n    double re;\n    double im;\n\n  ctypedef Complex C;",
            )
            .to_string();
        std::fs::File::create(&cython_header_path)
            .expect("Unable to modify bindings")
            .write_all(cython_header.as_bytes())
            .expect("Unable to modify bindings");
    };
    drop(cython_header);
}
