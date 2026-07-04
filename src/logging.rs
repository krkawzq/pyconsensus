use crate::htslib_ffi as ffi;
use std::sync::atomic::{AtomicBool, AtomicU8, Ordering};

#[repr(u8)]
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub enum LogLevel {
    Off = 0,
    Error = 1,
    Warn = 2,
    #[default]
    Info = 3,
    Debug = 4,
}

impl LogLevel {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "off" | "none" | "quiet" => Some(LogLevel::Off),
            "error" | "err" => Some(LogLevel::Error),
            "warn" | "warning" => Some(LogLevel::Warn),
            "info" => Some(LogLevel::Info),
            "debug" => Some(LogLevel::Debug),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            LogLevel::Off => "off",
            LogLevel::Error => "error",
            LogLevel::Warn => "warn",
            LogLevel::Info => "info",
            LogLevel::Debug => "debug",
        }
    }

    #[inline]
    fn from_u8(v: u8) -> Self {
        match v {
            0 => LogLevel::Off,
            1 => LogLevel::Error,
            2 => LogLevel::Warn,
            3 => LogLevel::Info,
            _ => LogLevel::Debug,
        }
    }
}

#[derive(Debug)]
pub struct LogControl {
    level: AtomicU8,
}

impl LogControl {
    pub fn new(level: LogLevel) -> Self {
        LogControl {
            level: AtomicU8::new(level as u8),
        }
    }

    #[inline]
    pub fn enabled(&self, level: LogLevel) -> bool {
        let current = self.level.load(Ordering::Relaxed);
        current != LogLevel::Off as u8 && current >= level as u8
    }

    #[inline]
    pub fn level(&self) -> LogLevel {
        LogLevel::from_u8(self.level.load(Ordering::Relaxed))
    }

    #[inline]
    pub fn set_level(&self, level: LogLevel) {
        self.level.store(level as u8, Ordering::Relaxed);
    }
}

impl Default for LogControl {
    fn default() -> Self {
        Self::new(LogLevel::default())
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum HtsLogLevel {
    Off,
    Error,
    Warn,
    #[default]
    Info,
    Debug,
    Trace,
}

impl HtsLogLevel {
    pub fn parse(s: &str) -> Option<Self> {
        match s.to_ascii_lowercase().as_str() {
            "off" | "none" | "quiet" => Some(HtsLogLevel::Off),
            "error" | "err" => Some(HtsLogLevel::Error),
            "warn" | "warning" => Some(HtsLogLevel::Warn),
            "info" => Some(HtsLogLevel::Info),
            "debug" => Some(HtsLogLevel::Debug),
            "trace" => Some(HtsLogLevel::Trace),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            HtsLogLevel::Off => "off",
            HtsLogLevel::Error => "error",
            HtsLogLevel::Warn => "warn",
            HtsLogLevel::Info => "info",
            HtsLogLevel::Debug => "debug",
            HtsLogLevel::Trace => "trace",
        }
    }

    #[inline]
    fn to_htslib(self) -> std::os::raw::c_int {
        match self {
            HtsLogLevel::Off => ffi::HTS_LOG_OFF,
            HtsLogLevel::Error => ffi::HTS_LOG_ERROR,
            HtsLogLevel::Warn => ffi::HTS_LOG_WARNING,
            HtsLogLevel::Info => ffi::HTS_LOG_INFO,
            HtsLogLevel::Debug => ffi::HTS_LOG_DEBUG,
            HtsLogLevel::Trace => ffi::HTS_LOG_TRACE,
        }
    }

    #[inline]
    fn from_htslib(level: std::os::raw::c_int) -> Self {
        match level {
            ffi::HTS_LOG_OFF => HtsLogLevel::Off,
            ffi::HTS_LOG_ERROR => HtsLogLevel::Error,
            ffi::HTS_LOG_WARNING => HtsLogLevel::Warn,
            ffi::HTS_LOG_INFO => HtsLogLevel::Info,
            ffi::HTS_LOG_DEBUG => HtsLogLevel::Debug,
            ffi::HTS_LOG_TRACE => HtsLogLevel::Trace,
            v if v <= ffi::HTS_LOG_OFF => HtsLogLevel::Off,
            v if v <= ffi::HTS_LOG_ERROR => HtsLogLevel::Error,
            v if v <= ffi::HTS_LOG_WARNING => HtsLogLevel::Warn,
            v if v <= ffi::HTS_LOG_INFO => HtsLogLevel::Info,
            v if v <= ffi::HTS_LOG_DEBUG => HtsLogLevel::Debug,
            _ => HtsLogLevel::Trace,
        }
    }
}

static HTSLIB_LOG_INITIALIZED: AtomicBool = AtomicBool::new(false);

/// htslib defaults to warning. pyconsensus makes the process-global default
/// info unless the user explicitly configured it first.
#[inline]
pub fn ensure_default_htslib_log_level() {
    if HTSLIB_LOG_INITIALIZED
        .compare_exchange(false, true, Ordering::Relaxed, Ordering::Relaxed)
        .is_ok()
    {
        unsafe {
            ffi::hts_set_log_level(HtsLogLevel::Info.to_htslib());
        }
    }
}

#[inline]
pub fn set_htslib_log_level(level: HtsLogLevel) {
    unsafe {
        ffi::hts_set_log_level(level.to_htslib());
    }
    HTSLIB_LOG_INITIALIZED.store(true, Ordering::Relaxed);
}

#[inline]
pub fn htslib_log_level() -> HtsLogLevel {
    ensure_default_htslib_log_level();
    let level = unsafe { ffi::hts_get_log_level() };
    HtsLogLevel::from_htslib(level)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_updates_log_level() {
        assert_eq!(LogLevel::parse("WARN"), Some(LogLevel::Warn));
        assert_eq!(LogLevel::parse("quiet"), Some(LogLevel::Off));
        assert_eq!(LogLevel::parse("verbose"), None);

        let log = LogControl::new(LogLevel::Info);
        assert!(log.enabled(LogLevel::Warn));
        assert!(log.enabled(LogLevel::Info));
        assert!(!log.enabled(LogLevel::Debug));
        log.set_level(LogLevel::Off);
        assert_eq!(log.level(), LogLevel::Off);
        assert!(!log.enabled(LogLevel::Error));
    }

    #[test]
    fn parses_htslib_log_level() {
        assert_eq!(HtsLogLevel::parse("warning"), Some(HtsLogLevel::Warn));
        assert_eq!(HtsLogLevel::parse("TRACE"), Some(HtsLogLevel::Trace));
        assert_eq!(HtsLogLevel::parse("verbose"), None);
    }
}
