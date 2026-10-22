//! Seeded random subsampling to a target number of bases.
//!
//! Reads are drawn in a random order until the target is reached, then kept in
//! their input order. The generator is implemented here (SplitMix64 feeding a
//! Fisher-Yates shuffle) rather than taken from a crate, so a given seed picks
//! the same reads in every build and on every platform.

pub struct SplitMix64(u64);

impl SplitMix64 {
    pub fn new(seed: u64) -> Self {
        SplitMix64(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `0..n` (Lemire's multiply-shift; the bias is negligible for read counts)
    pub fn below(&mut self, n: u64) -> u64 {
        ((self.next_u64() as u128 * n as u128) >> 64) as u64
    }
}

/// Indices (ascending) of the reads to keep so that their lengths sum to at
/// least `target` bases, or every index if the input is already that small.
pub fn choose(lengths: &[usize], target: u64, seed: u64) -> Vec<usize> {
    let total: u64 = lengths.iter().map(|&l| l as u64).sum();
    if target == 0 || total <= target {
        return (0..lengths.len()).collect();
    }
    let mut order: Vec<usize> = (0..lengths.len()).collect();
    let mut rng = SplitMix64::new(seed);
    for i in (1..order.len()).rev() {
        let j = rng.below(i as u64 + 1) as usize;
        order.swap(i, j);
    }
    let mut kept = Vec::new();
    let mut bases = 0u64;
    for i in order {
        if bases >= target {
            break;
        }
        kept.push(i);
        bases += lengths[i] as u64;
    }
    kept.sort_unstable();
    kept
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_everything_under_target() {
        assert_eq!(choose(&[10, 20, 30], 100, 1), vec![0, 1, 2]);
        assert_eq!(choose(&[10, 20, 30], 0, 1), vec![0, 1, 2]);
    }

    #[test]
    fn reaches_target_and_is_seeded() {
        let lengths: Vec<usize> = (0..1000).map(|i| 100 + i % 50).collect();
        let a = choose(&lengths, 20_000, 11);
        let b = choose(&lengths, 20_000, 11);
        let c = choose(&lengths, 20_000, 12);
        assert_eq!(a, b);
        assert_ne!(a, c);
        let bases: usize = a.iter().map(|&i| lengths[i]).sum();
        assert!((20_000..20_000 + 150).contains(&bases));
        assert!(a.windows(2).all(|w| w[0] < w[1]));
    }

    #[test]
    fn generator_is_stable() {
        // pinned so that a change of generator, which would silently change
        // which reads a seed picks, fails here
        let mut rng = SplitMix64::new(11);
        assert_eq!(rng.next_u64(), 0x50f5_647d_2380_309d);
    }
}
