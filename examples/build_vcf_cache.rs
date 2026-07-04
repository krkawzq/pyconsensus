use engine::vcf_store::VcfStore;
use std::env;
use std::ffi::OsString;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn cache_path(path: &Path) -> PathBuf {
    let mut p: OsString = path.as_os_str().to_os_string();
    p.push(".cvcf");
    PathBuf::from(p)
}

fn file_mb(path: &Path) -> Option<f64> {
    let bytes = std::fs::metadata(path).ok()?.len();
    Some(bytes as f64 / 1_000_000.0)
}

fn main() {
    let paths: Vec<PathBuf> = env::args_os().skip(1).map(PathBuf::from).collect();
    if paths.is_empty() {
        eprintln!("usage: build_vcf_cache <vcf-or-bcf> [<vcf-or-bcf> ...]");
        std::process::exit(2);
    }

    let mut failed = 0usize;
    for path in paths {
        let cache = cache_path(&path);
        let before = cache.exists();
        eprintln!(
            "CACHE_START path={} cache={} existing={}",
            path.display(),
            cache.display(),
            before
        );

        let t0 = Instant::now();
        match VcfStore::load(&path) {
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
                failed += 1;
                eprintln!("CACHE_ERROR path={} error={}", path.display(), err);
            }
        }
    }

    if failed != 0 {
        eprintln!("CACHE_FAILED n={failed}");
        std::process::exit(1);
    }
}
