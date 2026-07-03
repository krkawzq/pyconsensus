//! chain — `-c FILE` UCSC chain output.
//!
//! Ports `init_chain` / `push_chain_gap` / `print_chain`
//! (consensus.c:136/196/150). A chain records the ungapped blocks and gaps
//! between ref and alt (consensus) coordinates, in UCSC chain format.

/// One ungapped block + the gaps that follow it on ref/alt.
pub struct Chain {
    /// block length (ungapped run before this gap)
    block_lengths: Vec<i64>,
    /// ref gap length following the block
    ref_gaps: Vec<i64>,
    /// alt gap length following the block
    alt_gaps: Vec<i64>,
    pub ori_pos: i64,
    ref_last_block_ori: i64,
    alt_last_block_ori: i64,
    /// auto-increment chain id
    pub num_id: i64,
    /// total ref length (fa_length) for the final block / chain header.
    pub fa_length: i64,
    pub chr: String,
}

impl Chain {
    pub fn new(chr: String, ref_ori_pos: i64, fa_length: i64) -> Self {
        Chain {
            block_lengths: Vec::new(),
            ref_gaps: Vec::new(),
            alt_gaps: Vec::new(),
            ori_pos: ref_ori_pos,
            ref_last_block_ori: ref_ori_pos,
            alt_last_block_ori: ref_ori_pos,
            num_id: 0,
            fa_length,
            chr,
        }
    }

    /// Number of recorded gaps (blocks = gaps + 1 implicitly via the tail).
    pub fn n(&self) -> usize {
        self.block_lengths.len()
    }

    /// `push_chain_gap` (consensus.c:196). Records a ref↔alt gap. Back-to-back
    /// gaps (ref_start <= ref_last_block_ori) extend the previous block.
    pub fn push_gap(&mut self, ref_start: i64, ref_len: i64, alt_start: i64, alt_len: i64) {
        let num = self.block_lengths.len();
        if num > 0 && ref_start <= self.ref_last_block_ori {
            self.ref_last_block_ori = ref_start + ref_len;
            self.alt_last_block_ori = alt_start + alt_len;
            self.ref_gaps[num - 1] += ref_len;
            self.alt_gaps[num - 1] += alt_len;
        } else {
            self.block_lengths.push(ref_start - self.ref_last_block_ori);
            self.ref_gaps.push(ref_len);
            self.alt_gaps.push(alt_len);
            self.ref_last_block_ori = ref_start + ref_len;
            self.alt_last_block_ori = alt_start + alt_len;
        }
    }

    /// `print_chain` (consensus.c:150). Render to UCSC chain text. Increments
    /// the chain id (caller passes the running id; we keep it simple: use
    /// self.num_id and leave it incremented).
    pub fn render(&mut self) -> String {
        let ref_end_pos = self.fa_length + self.ori_pos;
        let last_block_size = ref_end_pos - self.ref_last_block_ori;
        let alt_end_pos = self.alt_last_block_ori + last_block_size;
        let mut score: i64 = 0;
        for &b in &self.block_lengths {
            score += b;
        }
        score += last_block_size;
        self.num_id += 1;

        let mut out = String::new();
        out.push_str(&format!(
            "chain {} {} {} + {} {} {} {} + {} {} {}\n",
            score,
            self.chr,
            ref_end_pos,
            self.ori_pos,
            ref_end_pos,
            self.chr,
            alt_end_pos,
            self.ori_pos,
            alt_end_pos,
            self.num_id
        ));
        for i in 0..self.block_lengths.len() {
            out.push_str(&format!(
                "{} {} {}\n",
                self.block_lengths[i], self.ref_gaps[i], self.alt_gaps[i]
            ));
        }
        out.push_str(&format!("{}\n\n", last_block_size));
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_chain_renders_header_and_tail() {
        let mut c = Chain::new("chr1".into(), 0, 100);
        let s = c.render();
        // score=100, ref_end=100, alt_end=100, one tail block of 100
        assert!(s.starts_with("chain 100 chr1 100 + 0 100 chr1 100 + 0 100 1\n"));
        assert!(s.ends_with("100\n\n"));
    }
}
