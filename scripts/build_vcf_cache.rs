use engine::logging::{set_htslib_log_level, HtsLogLevel, LogControl, LogLevel};
use engine::vcf_store::VcfStore;
use rayon::prelude::*;
use std::collections::HashSet;
use std::env;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::time::Instant;

struct Args {
    threads: Option<usize>,
    log_level: LogLevel,
    htslib_log_level: HtsLogLevel,
    paths: Vec<PathBuf>,
}

fn cache_path(path: &Path) -> PathBuf {
    VcfStore::default_cache_path_for(path)
}

fn absolute_path(path: &Path) -> PathBuf {
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        env::current_dir()
            .unwrap_or_else(|_| PathBuf::from("."))
            .join(path)
    }
}

fn normalize_path_for_dedup(path: &Path) -> PathBuf {
    let absolute = absolute_path(path);
    if let Ok(canonical) = absolute.canonicalize() {
        return canonical;
    }
    match (absolute.parent(), absolute.file_name()) {
        (Some(parent), Some(name)) => parent
            .canonicalize()
            .unwrap_or_else(|_| parent.to_path_buf())
            .join(name),
        _ => absolute,
    }
}

fn file_mb(path: &Path) -> Option<f64> {
    let bytes = std::fs::metadata(path).ok()?.len();
    Some(bytes as f64 / 1_000_000.0)
}

fn parse_args() -> Result<Args, String> {
    let mut threads = None;
    let mut log_level = LogLevel::Info;
    let mut htslib_log_level = HtsLogLevel::Info;
    let mut paths = Vec::new();
    let mut it = env::args_os().skip(1);
    while let Some(arg) = it.next() {
        if arg == "--threads" {
            let Some(value) = it.next() else {
                return Err("--threads requires a value".to_string());
            };
            let value = value
                .to_string_lossy()
                .parse::<usize>()
                .map_err(|e| format!("invalid --threads value: {e}"))?;
            threads = Some(value);
        } else if arg == "--log-level" {
            let Some(value) = it.next() else {
                return Err("--log-level requires a value".to_string());
            };
            log_level = LogLevel::parse(&value.to_string_lossy())
                .ok_or_else(|| "log level must be one of off,error,warn,info,debug".to_string())?;
        } else if arg == "--htslib-log-level" {
            let Some(value) = it.next() else {
                return Err("--htslib-log-level requires a value".to_string());
            };
            htslib_log_level = HtsLogLevel::parse(&value.to_string_lossy()).ok_or_else(|| {
                "htslib log level must be one of off,error,warn,info,debug,trace".to_string()
            })?;
        } else if arg == "-h" || arg == "--help" {
            return Err(String::new());
        } else {
            paths.push(PathBuf::from(arg));
        }
    }
    Ok(Args {
        threads,
        log_level,
        htslib_log_level,
        paths,
    })
}

fn thread_count(requested: Option<usize>, n: usize) -> usize {
    if n == 0 {
        return 1;
    }
    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    match requested {
        Some(v) if v > 0 => v.min(n).max(1),
        _ => available.min(n).max(1),
    }
}

fn usage() {
    eprintln!(
        "usage: build_vcf_cache [--threads N] [--log-level off|error|warn|info|debug] [--htslib-log-level off|error|warn|info|debug|trace] <vcf-or-bcf> [<vcf-or-bcf> ...]"
    );
}

fn main() {
    let args = match parse_args() {
        Ok(args) => args,
        Err(err) => {
            usage();
            if !err.is_empty() {
                eprintln!("{err}");
                std::process::exit(2);
            }
            return;
        }
    };
    if args.paths.is_empty() {
        usage();
        std::process::exit(2);
    }

    set_htslib_log_level(args.htslib_log_level);
    let log = Arc::new(LogControl::new(args.log_level));
    let mut seen = HashSet::new();
    let mut paths = Vec::new();
    for path in args.paths {
        let key = normalize_path_for_dedup(&cache_path(&path));
        if seen.insert(key) {
            paths.push(path);
        } else if log.enabled(LogLevel::Warn) {
            eprintln!("CACHE_SKIP_DUPLICATE path={}", path.display());
        }
    }

    let nthr = thread_count(args.threads, paths.len());
    if log.enabled(LogLevel::Info) {
        eprintln!("CACHE_BUILD_START n={} threads={}", paths.len(), nthr);
    }

    let failed = AtomicUsize::new(0);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(nthr)
        .build()
        .expect("failed to build rayon pool");
    pool.install(|| {
        paths.par_iter().for_each(|path| {
            let cache = cache_path(path);
            let before = cache.exists();
            if log.enabled(LogLevel::Info) {
                eprintln!(
                    "CACHE_START path={} cache={} existing={}",
                    path.display(),
                    cache.display(),
                    before
                );
            }

            let t0 = Instant::now();
            match VcfStore::load_with_log(path, None, log.as_ref()) {
                Ok(store) => {
                    let elapsed = t0.elapsed().as_secs_f64();
                    let cache_mb = file_mb(&cache).unwrap_or(0.0);
                    println!(
                        "CACHE_DONE path={} records={} samples={} cache={} cache_mb={:.1} elapsed_sec={:.3}",
                        path.display(),
                        store.n_records(),
                        store.sample_names().len(),
                        cache.display(),
                        cache_mb,
                        elapsed
                    );
                }
                Err(err) => {
                    failed.fetch_add(1, Ordering::Relaxed);
                    if log.enabled(LogLevel::Error) {
                        eprintln!("CACHE_ERROR path={} error={}", path.display(), err);
                    }
                }
            }
        });
    });

    let failed = failed.load(Ordering::Relaxed);
    if failed != 0 {
        eprintln!("CACHE_FAILED n={failed}");
        std::process::exit(1);
    }
}
