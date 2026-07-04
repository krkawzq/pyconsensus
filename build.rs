//! build.rs — statically link htslib (which bundles htscodecs).
//!
//! bcftools is NOT compiled or linked — it is only an algorithm reference.
//!
//! Strategy:
//!   * If `libs/htslib/libhts.a` already exists, reuse it (fast path).
//!   * Otherwise ensure a htslib checkout exists under `libs/htslib/`. If it
//!     is missing, clone htslib, check out `<HTSLIB_REF>`, and init the
//!     submodules. If present, reuse it; re-init the submodule if needed and,
//!     when the working tree is clean, sync HEAD to `<HTSLIB_REF>`.
//!   * Then run `autoreconf -i && ./configure ... && make -j libhts.a` once.
//!
//! The default ref is a tested commit (htslib 1.23.1-64-g9d53dcaa); override
//! with the `HTSLIB_REF` env var (any commit/tag/branch). The remote can be
//! overridden with `HTSLIB_URL`.
//!
//! We disable libcurl/s3/gcs/plugins: this tool only reads local files, so the
//! remote-URL machinery and its curl dependency are unnecessary. That yields a
//! minimal static libhts.a linked against system zlib/bz2/lzma/zstd.
//!
//! Network: the first build clones htslib from GitHub. On restricted networks
//! (e.g. the PJLab dev box), enable a proxy first — `labpon` or export
//! `http_proxy`/`https_proxy` — before `cargo build` / `maturin build`.

use std::env;
use std::path::{Path, PathBuf};
use std::process::Command;

/// Default htslib remote. Override with the `HTSLIB_URL` env var.
const HTSLIB_URL: &str = "https://github.com/samtools/htslib.git";
/// Pinned, tested htslib commit (1.23.1-64-g9d53dcaa). Override with
/// `HTSLIB_REF` — any ref: commit hash, tag, or branch name.
const HTSLIB_REF: &str = "9d53dcaa";

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").unwrap());
    let libs_dir = manifest_dir.join("libs");
    let htslib_dir = libs_dir.join("htslib");
    let libhts = htslib_dir.join("libhts.a");

    if !libhts.exists() {
        ensure_htslib(&htslib_dir);
        build_htslib(&htslib_dir);
    }

    // Sanity: make sure the archive is really there now.
    if !libhts.exists() {
        panic!(
            "libhts.a not found at {} after build. See the messages above.",
            libhts.display()
        );
    }

    // Compile the tiny C shim (struct-field accessors) into its own static
    // archive. OUT_DIR is provided by cargo.
    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    build_shim(&manifest_dir, &htslib_dir, &out_dir);

    // Link order matters for static archives: shim + libhts.a first, then the
    // system libs that resolve their undefined symbols.
    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=hts_shim");
    println!("cargo:rustc-link-search=native={}", htslib_dir.display());
    println!("cargo:rustc-link-lib=static=hts");

    // Statically link the compression libs that htslib pulls in (z, libdeflate,
    // bz2, lzma, zstd). Doing so keeps the resulting cdylib's NEEDED list down
    // to manylinux-allowlisted system libs (libm/libgcc_s/libc/ld-linux plus
    // pthread/dl below), so the wheel is self-contained and auditwheel has no
    // external library to repair. libdeflate.so.0 in particular is NOT on the
    // manylinux allow-list and breaks `auditwheel repair` on this host.
    //
    // These `.a` archives ship with the system's -dev packages; verify each is
    // present so a missing one fails loudly here instead of as a link error.
    for lib in ["z", "deflate", "bz2", "lzma", "zstd"] {
        ensure_static_lib(lib);
        println!("cargo:rustc-link-lib=static={}", lib);
    }
    // pthread/m/dl ARE on the manylinux allow-list — keep them dynamic.
    for lib in ["pthread", "m", "dl"] {
        println!("cargo:rustc-link-lib=dylib={}", lib);
    }

    println!("cargo:rerun-if-changed={}", libhts.display());
    println!(
        "cargo:rerun-if-changed={}",
        manifest_dir.join("src/hts_shim.c").display()
    );
    println!("cargo:rerun-if-changed=build.rs");
}

/// Make sure `dir` is a usable htslib checkout at `HTSLIB_REF` with the
/// bundled htscodecs submodule populated. Clones if absent; otherwise reuses
/// the existing checkout (re-inits the submodule, syncs the ref only when the
/// working tree is clean so local edits are never clobbered).
fn ensure_htslib(dir: &Path) {
    let url = env::var("HTSLIB_URL").unwrap_or_else(|_| HTSLIB_URL.to_string());
    let target_ref = env::var("HTSLIB_REF").unwrap_or_else(|_| HTSLIB_REF.to_string());

    if dir.exists() {
        if !is_git_repo(dir) {
            panic!(
                "`{}` exists but is not a git clone. Remove it so build.rs can \
                 clone htslib fresh:\n  rm -rf {}",
                dir.display(),
                dir.display()
            );
        }
        ensure_submodules(dir);
        checkout_if_needed(dir, &target_ref);
        return;
    }

    if let Some(parent) = dir.parent() {
        std::fs::create_dir_all(parent).unwrap_or_else(|e| {
            panic!("failed to create {}: {}", parent.display(), e);
        });
    }

    eprintln!("[build.rs] cloning htslib from {}", url);
    let dest = dir.to_str().expect("non-UTF-8 libs/htslib path");
    let display = format!("git clone {} {}", url, dest);
    let status = Command::new("git")
        .args(["clone", url.as_str(), dest])
        .status();
    check_status(status, &display);

    checkout_ref(dir, &target_ref);
    ensure_submodules(dir);
}

fn is_git_repo(dir: &Path) -> bool {
    dir.join(".git").exists()
}

/// `git submodule update --init --recursive` when the bundled htscodecs is
/// missing. Idempotent.
fn ensure_submodules(dir: &Path) {
    let bundled_htscodecs = dir.join("htscodecs").join("htscodecs");
    if bundled_htscodecs.exists() {
        return;
    }
    eprintln!("[build.rs] initializing htslib submodules (htscodecs)");
    run(
        dir,
        "git",
        &["submodule", "update", "--init", "--recursive"],
    );
}

/// Force-checkout `ref_` (used right after a fresh clone).
fn checkout_ref(dir: &Path, ref_: &str) {
    run(dir, "git", &["checkout", ref_]);
}

/// If HEAD != target_ref and the working tree is clean, check out target_ref.
/// A dirty tree is left untouched so we never clobber local edits.
fn checkout_if_needed(dir: &Path, ref_: &str) {
    let head = match git_output_opt(dir, &["rev-parse", "HEAD"]) {
        Some(h) => h,
        None => return,
    };
    let target = match git_output_opt(dir, &["rev-parse", ref_]) {
        Some(t) => t,
        None => {
            eprintln!(
                "[build.rs] could not resolve ref `{}` in libs/htslib; keeping \
                 HEAD ({})",
                ref_,
                head.trim()
            );
            return;
        }
    };
    if head == target {
        return;
    }
    let dirty = git_output_opt(dir, &["status", "--porcelain"])
        .map(|s| !s.trim().is_empty())
        .unwrap_or(true);
    if dirty {
        eprintln!(
            "[build.rs] libs/htslib HEAD ({}) != target ref `{}` but the working \
             tree is dirty; keeping current checkout. To pin `{}`, commit/stash \
             your changes and rerun.",
            head.trim(),
            ref_,
            ref_
        );
        return;
    }
    eprintln!(
        "[build.rs] checking out htslib `{}` (HEAD was {})",
        ref_,
        head.trim()
    );
    run(dir, "git", &["checkout", ref_]);
}

fn build_htslib(dir: &Path) {
    // htscodecs submodule must be populated (ensure_htslib did this; re-check
    // so a half-initialized checkout fails loudly here instead of in make).
    let bundled_htscodecs = dir.join("htscodecs").join("htscodecs");
    if !bundled_htscodecs.exists() {
        panic!(
            "htslib bundled htscodecs missing at {}. Run inside libs/htslib:\n  \
             git submodule update --init --recursive",
            bundled_htscodecs.display()
        );
    }

    let nproc = env::var("CARGO_BUILD_JOBS")
        .ok()
        .or_else(num_cpus_hint)
        .unwrap_or_else(|| "4".to_string());

    run(dir, "autoreconf", &["-i"]);
    // -fPIC is mandatory: the archive links into a cdylib (Python extension).
    let configure_status = Command::new("./configure")
        .args([
            "--disable-libcurl",
            "--disable-s3",
            "--disable-gcs",
            "--disable-plugins",
        ])
        .env("CFLAGS", "-O2 -fPIC")
        .current_dir(dir)
        .status();
    check_status(
        configure_status,
        "./configure --disable-libcurl --disable-s3 --disable-gcs --disable-plugins",
    );
    run(dir, "make", &["-j", &nproc, "libhts.a"]);
}

fn num_cpus_hint() -> Option<String> {
    // /proc/cpuinfo fallback; avoids pulling in the num_cpus crate.
    std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .map(|s| s.matches("processor\t:").count())
        .map(|c| c.to_string())
}

/// Panic with a actionable message if `lib<name>.a` is not on the system.
/// Statically linking htslib's compression deps requires the matching -dev
/// packages; a missing archive would otherwise surface as an opaque link error.
fn ensure_static_lib(name: &str) {
    let archive = format!("lib{}.a", name);
    let candidates = [
        "/usr/lib/x86_64-linux-gnu",
        "/lib/x86_64-linux-gnu",
        "/usr/local/lib",
        "/usr/lib",
        "/lib",
    ];
    let found = candidates
        .iter()
        .any(|dir| Path::new(dir).join(&archive).exists());
    if !found {
        let pkg = match name {
            "z" => "zlib1g-dev",
            "deflate" => "libdeflate-dev",
            "bz2" => "libbz2-dev",
            "lzma" => "liblzma-dev",
            "zstd" => "libzstd-dev",
            other => other,
        };
        panic!(
            "static library `{}` not found in standard search paths. Install \
             the matching -dev package, e.g.:\n  sudo apt install {}",
            archive, pkg
        );
    }
}

fn build_shim(manifest_dir: &Path, htslib_dir: &Path, out_dir: &Path) {
    let src = manifest_dir.join("src").join("hts_shim.c");
    let obj = out_dir.join("hts_shim.o");
    let archive = out_dir.join("libhts_shim.a");

    let display = format!(
        "gcc -O2 -fPIC -I{} -c {} -o {}",
        htslib_dir.display(),
        src.display(),
        obj.display()
    );
    eprintln!("[build.rs] running: {}", display);
    let status = Command::new("gcc")
        .args(["-O2", "-fPIC", "-I"])
        .arg(htslib_dir)
        .arg("-c")
        .arg(&src)
        .arg("-o")
        .arg(&obj)
        .status();
    check_status(status, &display);

    let display = format!("ar rcs {} {}", archive.display(), obj.display());
    eprintln!("[build.rs] running: {}", display);
    let status = Command::new("ar")
        .arg("rcs")
        .arg(&archive)
        .arg(&obj)
        .status();
    check_status(status, &display);
}

fn run(dir: &Path, program: &str, args: &[&str]) {
    let display = format!("{} {}", program, args.join(" "));
    eprintln!("[build.rs] running: {} (in {})", display, dir.display());
    let status = Command::new(program).args(args).current_dir(dir).status();
    check_status(status, &display);
}

/// Run a git subcommand and return its stdout, or `None` on failure.
fn git_output_opt(dir: &Path, args: &[&str]) -> Option<String> {
    let out = Command::new("git")
        .args(args)
        .current_dir(dir)
        .output()
        .ok()?;
    if !out.status.success() {
        return None;
    }
    Some(String::from_utf8_lossy(&out.stdout).to_string())
}

fn check_status(status: std::io::Result<std::process::ExitStatus>, display: &str) {
    let status = status.unwrap_or_else(|e| panic!("failed to execute `{}`: {}", display, e));
    if !status.success() {
        panic!("`{}` failed with status {}", display, status);
    }
}
