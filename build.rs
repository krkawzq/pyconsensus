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
//! remote-URL machinery and its curl dependency are unnecessary.
//!
//! htslib's compression deps (zlib, bzip2, liblzma, libdeflate, zstd) are also
//! built from source with `-fPIC` and statically linked. The system -dev `.a`
//! archives are non-PIC and cannot link into a cdylib; building our own keeps
//! the wheel self-contained and manylinux-clean (no external `.so` to repair).
//!
//! Network: the first build clones htslib and downloads the compression-lib
//! tarballs. On restricted networks (e.g. the PJLab dev box), enable a proxy
//! first — `labpon` or export `http_proxy`/`https_proxy` — before `cargo
//! build` / `maturin build`.

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

    // PIC static archives for htslib's compression dependencies, built from
    // source with -fPIC (see ensure_* below). The system -dev `.a` archives are
    // non-PIC and fail with R_X86_64_PC32 relocation errors when linked into a
    // shared object, so we compile our own under libs/. Each ensure_* is a
    // no-op if the archive already exists (fast incremental rebuilds).
    let zlib_a = ensure_zlib(&libs_dir);
    let bz2_a = ensure_bzip2(&libs_dir);
    let lzma_a = ensure_lzma(&libs_dir);
    let deflate_a = ensure_libdeflate(&libs_dir);
    let zstd_a = ensure_zstd(&libs_dir);

    // Link order: shim + libhts.a first, then the PIC compression archives,
    // then the manylinux-allowlisted system libs (kept dynamic).
    println!("cargo:rustc-link-search=native={}", out_dir.display());
    println!("cargo:rustc-link-lib=static=hts_shim");
    println!("cargo:rustc-link-search=native={}", htslib_dir.display());
    println!("cargo:rustc-link-lib=static=hts");
    for (archive, name) in [
        (&zlib_a, "z"),
        (&bz2_a, "bz2"),
        (&lzma_a, "lzma"),
        (&deflate_a, "deflate"),
        (&zstd_a, "zstd"),
    ] {
        let dir = archive.parent().expect("archive has no parent dir");
        println!("cargo:rustc-link-search=native={}", dir.display());
        println!("cargo:rustc-link-lib=static={}", name);
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

// --- PIC static compression libs (compiled from source with -fPIC) ---------
// The system -dev `.a` archives (zlib1g-dev, libbz2-dev, ...) are built without
// -fPIC and cannot be linked into a cdylib: they trigger R_X86_64_PC32
// relocation errors against symbols like `stderr`. We download release
// tarballs and build PIC archives under libs/<name>/. Override the source URL
// with the `<NAME>_URL` env var. A proxy (`labpon`) is required on the PJLab
// dev box for the first build — curl inherits http_proxy/https_proxy from the
// environment, so run `labpon` before `cargo build` / `maturin build`.

const ZLIB_URL: &str = "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz";
const BZIP2_URL: &str = "https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz";
const LZMA_URL: &str =
    "https://github.com/tukaani-project/xz/releases/download/v5.6.2/xz-5.6.2.tar.gz";
const LIBDEFLATE_URL: &str =
    "https://github.com/ebiggers/libdeflate/releases/download/v1.22/libdeflate-1.22.tar.gz";
const ZSTD_URL: &str =
    "https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz";

/// Parallelism hint for `make -j`. Prefers cargo's job count, falls back to
/// /proc/cpuinfo, then 4.
fn nproc() -> String {
    env::var("CARGO_BUILD_JOBS")
        .ok()
        .or_else(num_cpus_hint)
        .unwrap_or_else(|| "4".to_string())
}

/// Download `url` to `dest` with curl, verifying the tarball is intact before
/// returning. No-op if `dest` already exists. curl honors http_proxy/
/// https_proxy from the env. Up to 3 attempts: a truncated download (curl
/// exits 0 on a proxy hiccup but produces a short file) is caught by the
/// `tar -tf` integrity probe and retried after deleting the partial file.
fn download(url: &str, dest: &Path) {
    if dest.exists() {
        return;
    }
    if let Some(parent) = dest.parent() {
        std::fs::create_dir_all(parent).unwrap_or_else(|e| {
            panic!("failed to create {}: {}", parent.display(), e);
        });
    }
    for attempt in 1..=3 {
        eprintln!("[build.rs] downloading {} (attempt {})", url, attempt);
        let display = format!(
            "curl -fL --retry 3 --max-time 300 -o {} {}",
            dest.display(),
            url
        );
        let status = Command::new("curl")
            .args(["-fL", "--retry", "3", "--max-time", "300", "-o"])
            .arg(dest)
            .arg(url)
            .status();
        check_status(status, &display);
        // Integrity check: list the archive; a truncated gzip/xz fails here.
        let probe = Command::new("tar")
            .args(["-tf", dest.to_str().expect("non-UTF-8 tarball path")])
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status();
        if probe.map(|s| s.success()).unwrap_or(false) {
            return;
        }
        eprintln!("[build.rs] tarball failed integrity probe, retrying");
        let _ = std::fs::remove_file(dest);
    }
    panic!(
        "failed to download a valid tarball from {} after 3 attempts",
        url
    );
}

/// Extract a tarball into `dest_dir`, stripping the top-level `name-1.2.3/`
/// component so source files land directly in dest_dir. GNU tar auto-detects
/// gzip/bzip2/xz compression.
fn extract_tarball(tarball: &Path, dest_dir: &Path) {
    std::fs::create_dir_all(dest_dir).unwrap_or_else(|e| {
        panic!("failed to create {}: {}", dest_dir.display(), e);
    });
    run(
        dest_dir,
        "tar",
        &[
            "-xf",
            tarball.to_str().expect("non-UTF-8 tarball path"),
            "--strip-components=1",
            "--auto-compress",
        ],
    );
}

/// Run `program args...` in `dir` with extra environment variables.
fn run_env(dir: &Path, program: &str, args: &[&str], envs: &[(&str, &str)]) {
    let display = format!("{} {}", program, args.join(" "));
    eprintln!("[build.rs] running: {} (in {})", display, dir.display());
    let mut cmd = Command::new(program);
    cmd.args(args).current_dir(dir);
    for (k, v) in envs {
        cmd.env(k, v);
    }
    let status = cmd.status();
    check_status(status, &display);
}

/// `libz.a` built with -fPIC from zlib 1.3.1. configure writes CFLAGS into the
/// generated Makefile, so `make` picks it up without a further override.
fn ensure_zlib(libs_dir: &Path) -> PathBuf {
    let dir = libs_dir.join("zlib");
    let archive = dir.join("libz.a");
    if archive.exists() {
        return archive;
    }
    let url = env::var("ZLIB_URL").unwrap_or_else(|_| ZLIB_URL.to_string());
    let tarball = libs_dir.join("zlib.tar.gz");
    download(&url, &tarball);
    extract_tarball(&tarball, &dir);
    run_env(
        &dir,
        "./configure",
        &["--static"],
        &[("CFLAGS", "-O2 -fPIC")],
    );
    run(&dir, "make", &["-j", &nproc()]);
    assert!(
        archive.exists(),
        "libz.a not produced at {}",
        archive.display()
    );
    archive
}

/// `libbz2.a` built with -fPIC from bzip2 1.0.8. bzip2 has no configure step;
/// CC/CFLAGS are passed as make command-line variables so they override the
/// Makefile's own assignments.
fn ensure_bzip2(libs_dir: &Path) -> PathBuf {
    let dir = libs_dir.join("bzip2");
    let archive = dir.join("libbz2.a");
    if archive.exists() {
        return archive;
    }
    let url = env::var("BZIP2_URL").unwrap_or_else(|_| BZIP2_URL.to_string());
    let tarball = libs_dir.join("bzip2.tar.gz");
    download(&url, &tarball);
    extract_tarball(&tarball, &dir);
    run(
        &dir,
        "make",
        &["-j", &nproc(), "CC=gcc", "CFLAGS=-O2 -fPIC", "libbz2.a"],
    );
    assert!(
        archive.exists(),
        "libbz2.a not produced at {}",
        archive.display()
    );
    archive
}

/// `liblzma.a` built with -fPIC from xz 5.6.2. The release tarball ships a
/// pre-generated configure, so no autoreconf is needed.
fn ensure_lzma(libs_dir: &Path) -> PathBuf {
    let dir = libs_dir.join("xz");
    let archive = dir.join("src/liblzma/.libs/liblzma.a");
    if archive.exists() {
        return archive;
    }
    let url = env::var("LZMA_URL").unwrap_or_else(|_| LZMA_URL.to_string());
    let tarball = libs_dir.join("xz.tar.gz");
    download(&url, &tarball);
    extract_tarball(&tarball, &dir);
    run_env(
        &dir,
        "./configure",
        &[
            "--disable-shared",
            "--enable-static",
            "--disable-xz",
            "--disable-xzdec",
            "--disable-lzmadec",
            "--disable-lzmainfo",
            "--disable-scripts",
            "--disable-doc",
        ],
        &[("CFLAGS", "-O2 -fPIC")],
    );
    run(&dir, "make", &["-j", &nproc()]);
    assert!(
        archive.exists(),
        "liblzma.a not produced at {}",
        archive.display()
    );
    archive
}

/// `libdeflate.a` built with -fPIC from libdeflate 1.22. Since 1.22 the
/// upstream Makefile was replaced by CMake, so we drive cmake directly: static
/// lib only, no shared lib / gzip program / tests, with CMAKE_C_FLAGS=-fPIC.
fn ensure_libdeflate(libs_dir: &Path) -> PathBuf {
    let dir = libs_dir.join("libdeflate");
    let build_dir = dir.join("build");
    let archive = build_dir.join("libdeflate.a");
    if archive.exists() {
        return archive;
    }
    let url = env::var("LIBDEFLATE_URL").unwrap_or_else(|_| LIBDEFLATE_URL.to_string());
    let tarball = libs_dir.join("libdeflate.tar.gz");
    download(&url, &tarball);
    extract_tarball(&tarball, &dir);
    std::fs::create_dir_all(&build_dir).unwrap_or_else(|e| {
        panic!("failed to create {}: {}", build_dir.display(), e);
    });
    run_env(
        &build_dir,
        "cmake",
        &[
            dir.to_str().expect("non-UTF-8 libdeflate path"),
            "-DCMAKE_BUILD_TYPE=Release",
            "-DLIBDEFLATE_BUILD_STATIC_LIB=ON",
            "-DLIBDEFLATE_BUILD_SHARED_LIB=OFF",
            "-DLIBDEFLATE_BUILD_GZIP=OFF",
            "-DLIBDEFLATE_BUILD_TESTS=OFF",
            "-DCMAKE_C_FLAGS=-O2 -fPIC",
        ],
        &[],
    );
    run(&build_dir, "cmake", &["--build", ".", "-j", &nproc()]);
    assert!(
        archive.exists(),
        "libdeflate.a not produced at {}",
        archive.display()
    );
    archive
}

/// `libzstd.a` built with -fPIC from zstd 1.5.6 (built in the lib/ subdir).
fn ensure_zstd(libs_dir: &Path) -> PathBuf {
    let dir = libs_dir.join("zstd");
    let archive = dir.join("lib/libzstd.a");
    if archive.exists() {
        return archive;
    }
    let url = env::var("ZSTD_URL").unwrap_or_else(|_| ZSTD_URL.to_string());
    let tarball = libs_dir.join("zstd.tar.gz");
    download(&url, &tarball);
    extract_tarball(&tarball, &dir);
    let lib_dir = dir.join("lib");
    run(
        &lib_dir,
        "make",
        &["-j", &nproc(), "CFLAGS=-O2 -fPIC", "libzstd.a"],
    );
    assert!(
        archive.exists(),
        "libzstd.a not produced at {}",
        archive.display()
    );
    archive
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
