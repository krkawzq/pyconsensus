//! stats — lightweight counters for fastpath-first execution.
//!
//! Counters are intentionally plain Rust data structures. They can be carried
//! by tests, benchmark harnesses, or future Python diagnostics without pulling
//! atomics into the hot path unless cross-thread aggregation is explicitly
//! needed.

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FastPathLane {
    EmptyRegion = 0,
    SameLenOnly = 1,
    SameLenIupac = 2,
    BiallelicPhasedBatch = 3,
    NormalizedEditScript = 4,
    MixedSimpleEdits = 5,
    FallbackStateMachine = 6,
}

impl FastPathLane {
    pub const COUNT: usize = 7;

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            FastPathLane::EmptyRegion => "EmptyRegion",
            FastPathLane::SameLenOnly => "SameLenOnly",
            FastPathLane::SameLenIupac => "SameLenIupac",
            FastPathLane::BiallelicPhasedBatch => "BiallelicPhasedBatch",
            FastPathLane::NormalizedEditScript => "NormalizedEditScript",
            FastPathLane::MixedSimpleEdits => "MixedSimpleEdits",
            FastPathLane::FallbackStateMachine => "FallbackStateMachine",
        }
    }
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FallbackReason {
    ChainEnabled = 0,
    MarkEnabled = 1,
    MaskOverlap = 2,
    VariantOverlap = 3,
    SymbolicAllele = 4,
    RefMismatch = 5,
    ComplexAllele = 6,
    MissingGt = 7,
    LengthChangingEdit = 8,
    UnsupportedMode = 9,
}

impl FallbackReason {
    pub const COUNT: usize = 10;

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            FallbackReason::ChainEnabled => "ChainEnabled",
            FallbackReason::MarkEnabled => "MarkEnabled",
            FallbackReason::MaskOverlap => "MaskOverlap",
            FallbackReason::VariantOverlap => "VariantOverlap",
            FallbackReason::SymbolicAllele => "SymbolicAllele",
            FallbackReason::RefMismatch => "RefMismatch",
            FallbackReason::ComplexAllele => "ComplexAllele",
            FallbackReason::MissingGt => "MissingGt",
            FallbackReason::LengthChangingEdit => "LengthChangingEdit",
            FallbackReason::UnsupportedMode => "UnsupportedMode",
        }
    }
}

#[repr(u8)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ApplyFailureKind {
    RefMismatch = 0,
    BrokenVcf = 1,
    InvalidOverlapTrim = 2,
    UnsupportedSymbolicAllele = 3,
}

impl ApplyFailureKind {
    pub const COUNT: usize = 4;

    #[inline]
    pub fn as_usize(self) -> usize {
        self as usize
    }

    pub fn name(self) -> &'static str {
        match self {
            ApplyFailureKind::RefMismatch => "RefMismatch",
            ApplyFailureKind::BrokenVcf => "BrokenVcf",
            ApplyFailureKind::InvalidOverlapTrim => "InvalidOverlapTrim",
            ApplyFailureKind::UnsupportedSymbolicAllele => "UnsupportedSymbolicAllele",
        }
    }
}

#[derive(Clone, Debug, Default)]
pub struct RuntimeStats {
    pub regions_total: u64,
    pub tasks_total: u64,
    pub records_seen: u64,
    pub lane_counts: [u64; FastPathLane::COUNT],
    pub fallback_reasons: [u64; FallbackReason::COUNT],
    pub same_len_fastpath_records: u64,
    pub edit_script_fastpath_records: u64,
    pub fallback_records: u64,
    pub apply_failures: [u64; ApplyFailureKind::COUNT],
    pub alloc_bytes: u64,
}

impl RuntimeStats {
    pub fn observe_region(&mut self) {
        self.regions_total += 1;
    }

    pub fn observe_task(&mut self) {
        self.tasks_total += 1;
    }

    pub fn observe_tasks(&mut self, n: u64) {
        self.tasks_total += n;
    }

    pub fn observe_lane(&mut self, lane: FastPathLane) {
        self.lane_counts[lane.as_usize()] += 1;
    }

    pub fn observe_record(&mut self) {
        self.records_seen += 1;
    }

    pub fn observe_records(&mut self, n: u64) {
        self.records_seen += n;
    }

    pub fn observe_same_len_fastpath(&mut self) {
        self.same_len_fastpath_records += 1;
    }

    pub fn observe_same_len_fastpath_records(&mut self, n: u64) {
        self.same_len_fastpath_records += n;
    }

    pub fn observe_edit_script_fastpath(&mut self) {
        self.edit_script_fastpath_records += 1;
    }

    pub fn observe_edit_script_fastpath_records(&mut self, n: u64) {
        self.edit_script_fastpath_records += n;
    }

    pub fn observe_fallback_records(&mut self, n: u64) {
        self.fallback_records += n;
    }

    pub fn observe_alloc_bytes(&mut self, n: u64) {
        self.alloc_bytes += n;
    }

    pub fn observe_fallback_reason(&mut self, reason: FallbackReason) {
        self.fallback_reasons[reason.as_usize()] += 1;
    }

    pub fn observe_apply_failure(&mut self, kind: ApplyFailureKind) {
        self.apply_failures[kind.as_usize()] += 1;
    }

    pub fn apply_failure_count(&self, kind: ApplyFailureKind) -> u64 {
        self.apply_failures[kind.as_usize()]
    }

    pub fn observe_fallback(&mut self, reason: FallbackReason) {
        self.observe_fallback_records(1);
        self.observe_fallback_reason(reason);
        self.observe_lane(FastPathLane::FallbackStateMachine);
    }

    pub fn lane_count(&self, lane: FastPathLane) -> u64 {
        self.lane_counts[lane.as_usize()]
    }

    pub fn fallback_reason_count(&self, reason: FallbackReason) -> u64 {
        self.fallback_reasons[reason.as_usize()]
    }

    pub fn merge(&mut self, other: RuntimeStats) {
        self.regions_total += other.regions_total;
        self.tasks_total += other.tasks_total;
        self.records_seen += other.records_seen;
        for (dst, src) in self.lane_counts.iter_mut().zip(other.lane_counts) {
            *dst += src;
        }
        for (dst, src) in self.fallback_reasons.iter_mut().zip(other.fallback_reasons) {
            *dst += src;
        }
        for (dst, src) in self.apply_failures.iter_mut().zip(other.apply_failures) {
            *dst += src;
        }
        self.same_len_fastpath_records += other.same_len_fastpath_records;
        self.edit_script_fastpath_records += other.edit_script_fastpath_records;
        self.fallback_records += other.fallback_records;
        self.alloc_bytes += other.alloc_bytes;
    }

    pub fn summary_lines(&self) -> Vec<String> {
        let mut lines = vec![
            format!("regions_total={}", self.regions_total),
            format!("tasks_total={}", self.tasks_total),
            format!("records_seen={}", self.records_seen),
            format!(
                "same_len_fastpath_records={}",
                self.same_len_fastpath_records
            ),
            format!(
                "edit_script_fastpath_records={}",
                self.edit_script_fastpath_records
            ),
            format!("fallback_records={}", self.fallback_records),
            format!("alloc_bytes={}", self.alloc_bytes),
        ];
        for lane in [
            FastPathLane::EmptyRegion,
            FastPathLane::SameLenOnly,
            FastPathLane::SameLenIupac,
            FastPathLane::BiallelicPhasedBatch,
            FastPathLane::NormalizedEditScript,
            FastPathLane::MixedSimpleEdits,
            FastPathLane::FallbackStateMachine,
        ] {
            lines.push(format!("lane.{}={}", lane.name(), self.lane_count(lane)));
        }
        for reason in [
            FallbackReason::ChainEnabled,
            FallbackReason::MarkEnabled,
            FallbackReason::MaskOverlap,
            FallbackReason::VariantOverlap,
            FallbackReason::SymbolicAllele,
            FallbackReason::RefMismatch,
            FallbackReason::ComplexAllele,
            FallbackReason::MissingGt,
            FallbackReason::LengthChangingEdit,
            FallbackReason::UnsupportedMode,
        ] {
            lines.push(format!(
                "fallback.{}={}",
                reason.name(),
                self.fallback_reason_count(reason)
            ));
        }
        for kind in [
            ApplyFailureKind::RefMismatch,
            ApplyFailureKind::BrokenVcf,
            ApplyFailureKind::InvalidOverlapTrim,
            ApplyFailureKind::UnsupportedSymbolicAllele,
        ] {
            lines.push(format!(
                "apply_failure.{}={}",
                kind.name(),
                self.apply_failure_count(kind)
            ));
        }
        lines
    }
}
