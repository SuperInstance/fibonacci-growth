//! Fibonacci growth: Penrose inflation outward, Mandelbrot roughness inward.
//! CR is scale-invariant by construction.

/// Fibonacci sequence generator with arbitrary starting values
pub fn fibonacci_sequence(a: u64, b: u64, n: usize) -> Vec<u64> {
    if n == 0 { return vec![]; }
    if n == 1 { return vec![a]; }
    let mut seq = vec![a, b];
    for i in 2..n {
        let next = seq[i - 1].checked_add(seq[i - 2]).unwrap_or(u64::MAX);
        seq.push(next);
    }
    seq.truncate(n);
    seq
}

/// Standard Fibonacci starting from (1, 1)
pub fn standard_fibonacci(n: usize) -> Vec<u64> {
    fibonacci_sequence(1, 1, n)
}

/// Golden ratio approximation from Fibonacci numbers: F(n+1)/F(n)
pub fn golden_ratio_approx(n: usize) -> f64 {
    if n < 2 { return 1.0; }
    let seq = standard_fibonacci(n + 1);
    seq[n] as f64 / seq[n - 1] as f64
}

/// Penrose inflation: given an initial patch of tiles (represented as a sequence
/// of edge types A=1, B=2), inflate by substituting according to Penrose rules.
/// A → AB (deflation: short becomes long + short)
/// B → A (deflation: long becomes short)
/// We track ratio A/B which converges to φ.
pub fn penrose_inflate(edges: &[u8], steps: usize) -> Vec<u8> {
    let mut current = edges.to_vec();
    for _ in 0..steps {
        let mut next = Vec::with_capacity(current.len() * 2);
        for &edge in &current {
            match edge {
                1 => { next.push(1); next.push(2); } // A → AB
                _ => { next.push(1); }                // B → A
            }
        }
        current = next;
    }
    current
}

/// Ratio of edge type 1 to edge type 2 in a sequence (converges to φ)
pub fn edge_ratio(edges: &[u8]) -> f64 {
    let a = edges.iter().filter(|&&e| e == 1).count();
    let b = edges.iter().filter(|&&e| e == 2).count();
    if b == 0 { return if a > 0 { f64::INFINITY } else { 0.0 }; }
    a as f64 / b as f64
}

/// Mandelbrot roughness: estimate the fractal dimension via box-counting on
/// the escape boundary. Returns approximate Hausdorff dimension.
pub fn mandelbrot_roughness(max_iter: u32, _samples: usize) -> f64 {
    // Scan a focused region around the main cardioid/period-2 bulb boundary
    let resolution = 200usize;
    let r_min = -1.5; let r_max = 0.5;
    let i_min = -1.0; let i_max = 1.0;
    let dr = (r_max - r_min) / resolution as f64;
    let di = (i_max - i_min) / resolution as f64;

    // Compute escape times on grid
    let mut escape = vec![vec![0u32; resolution]; resolution];
    for ri in 0..resolution {
        for ii in 0..resolution {
            let cr = r_min + ri as f64 * dr;
            let ci = i_min + ii as f64 * di;
            let (_, _, it) = mandelbrot_iterate(cr, ci, max_iter);
            escape[ri][ii] = it;
        }
    }

    // Extract boundary: pixels where a neighbor has different in/out status
    let mut boundary: Vec<(usize, usize)> = Vec::new();
    for ri in 1..resolution-1 {
        for ii in 1..resolution-1 {
            let inside = escape[ri][ii] >= max_iter;
            let has_diff = (escape[ri-1][ii] >= max_iter) != inside
                || (escape[ri+1][ii] >= max_iter) != inside
                || (escape[ri][ii-1] >= max_iter) != inside
                || (escape[ri][ii+1] >= max_iter) != inside;
            if has_diff { boundary.push((ri, ii)); }
        }
    }
    if boundary.is_empty() { return 1.0; }

    // Box-counting at pixel scales
    let scales: [usize; 4] = [1, 2, 4, 8];
    let mut counts = Vec::new();
    for &s in &scales {
        let mut boxes = std::collections::HashSet::new();
        for &(r, c) in &boundary {
            boxes.insert((r / s, c / s));
        }
        counts.push((s, boxes.len() as f64));
    }

    // Regression: log(N) vs log(1/scale)
    let pts: Vec<(f64, f64)> = counts.iter()
        .filter(|(_, n)| *n > 0.0)
        .map(|&(s, n)| (1.0 / s as f64, n))
        .map(|(inv_s, n)| (inv_s.ln(), n.ln()))
        .collect();
    if pts.len() < 2 { return 1.0; }
    let n = pts.len() as f64;
    let sx: f64 = pts.iter().map(|(x, _)| *x).sum();
    let sy: f64 = pts.iter().map(|(_, y)| *y).sum();
    let sxx: f64 = pts.iter().map(|(x, _)| x * x).sum();
    let sxy: f64 = pts.iter().map(|(x, y)| x * y).sum();
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-30 { return 1.0; }
    let slope = (n * sxy - sx * sy) / denom;
    slope.max(0.0).min(2.5)
}

/// Iterate z = z² + c, return final z and iteration count
fn mandelbrot_iterate(cr: f64, ci: f64, max_iter: u32) -> (f64, f64, u32) {
    let mut zr = 0.0_f64;
    let mut zi = 0.0_f64;
    for i in 0..max_iter {
        let zr2 = zr * zr;
        let zi2 = zi * zi;
        if zr2 + zi2 > 4.0 { return (zr, zi, i); }
        zi = 2.0 * zr * zi + ci;
        zr = zr2 - zi2 + cr;
    }
    (zr, zi, max_iter)
}

/// Conservation ratio at Fibonacci scales.
/// For a graph with Fibonacci-structured degrees, compute how CR
/// changes across scales.
pub fn fibonacci_cr_profile(adj: &[Vec<f64>], levels: usize) -> Vec<f64> {
    let n = adj.len();
    if n == 0 { return vec![]; }
    let mut profile = Vec::new();
    let lap = laplacian(adj);
    let eigenvalues = jacobi_eigenvalues(lap);

    // Base CR: λ₂/λ_n (algebraic connectivity / max eigenvalue)
    let lambda_2 = eigenvalues.get(1).copied().unwrap_or(0.0);
    let lambda_n = eigenvalues.last().copied().unwrap_or(1.0);
    let base_cr = if lambda_n > 0.0 { lambda_2 / lambda_n } else { 0.0 };
    profile.push(base_cr);

    // At each level, coarsen by level using Fibonacci-based aggregation
    for level in 1..=levels {
        let fib = standard_fibonacci(level + 2);
        let block_size = fib.last().copied().unwrap_or(1) as usize;
        if block_size >= n { break; }
        let coarse = coarsen_graph(adj, block_size);
        if coarse.is_empty() { break; }
        let lap = laplacian(&coarse);
        let eigs = jacobi_eigenvalues(lap);
        let l2 = eigs.get(1).copied().unwrap_or(0.0);
        let ln = eigs.last().copied().unwrap_or(1.0);
        let cr = if ln > 0.0 { l2 / ln } else { 0.0 };
        profile.push(cr);
    }
    profile
}

/// Build Laplacian from adjacency matrix
fn laplacian(adj: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let n = adj.len();
    let mut lap = vec![vec![0.0_f64; n]; n];
    for i in 0..n {
        let deg: f64 = adj[i].iter().sum();
        lap[i][i] = deg;
        for j in 0..n {
            if i != j { lap[i][j] = -adj[i][j]; }
        }
    }
    lap
}

/// Coarsen a graph by aggregating nodes into blocks of given size
fn coarsen_graph(adj: &[Vec<f64>], block_size: usize) -> Vec<Vec<f64>> {
    let n = adj.len();
    let num_blocks = (n + block_size - 1) / block_size;
    let mut coarse = vec![vec![0.0_f64; num_blocks]; num_blocks];
    for bi in 0..num_blocks {
        for bj in 0..num_blocks {
            let mut total = 0.0;
            let si = bi * block_size;
            let ei = ((bi + 1) * block_size).min(n);
            let sj = bj * block_size;
            let ej = ((bj + 1) * block_size).min(n);
            for i in si..ei {
                for j in sj..ej {
                    if bi != bj || i != j { total += adj[i][j]; }
                }
            }
            coarse[bi][bj] = total;
        }
    }
    coarse
}

/// Jacobi eigenvalue decomposition
fn jacobi_eigenvalues(mut a: Vec<Vec<f64>>) -> Vec<f64> {
    let n = a.len();
    if n == 0 { return vec![]; }
    for _ in 0..100 * n * n {
        let (mut p, mut q) = (0, 1);
        let mut max_val = 0.0_f64;
        for i in 0..n {
            for j in (i + 1)..n {
                if a[i][j].abs() > max_val {
                    max_val = a[i][j].abs();
                    p = i; q = j;
                }
            }
        }
        if max_val < 1e-14 { break; }
        let app = a[p][p]; let aqq = a[q][q]; let apq = a[p][q];
        let theta = if (app - aqq).abs() < 1e-30 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * (2.0 * apq / (app - aqq)).atan()
        };
        let (c, s) = (theta.cos(), theta.sin());
        for i in 0..n {
            if i != p && i != q {
                let aip = a[i][p]; let aiq = a[i][q];
                a[i][p] = c * aip + s * aiq; a[p][i] = a[i][p];
                a[i][q] = -s * aip + c * aiq; a[q][i] = a[i][q];
            }
        }
        a[p][p] = c * c * app + 2.0 * s * c * apq + s * s * aqq;
        a[q][q] = s * s * app - 2.0 * s * c * apq + c * c * aqq;
        a[p][q] = 0.0; a[q][p] = 0.0;
    }
    let mut eigs: Vec<f64> = (0..n).map(|i| a[i][i]).collect();
    eigs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    eigs
}

/// Zeckendorf representation: express n as sum of non-consecutive Fibonacci numbers
pub fn zeckendorf(n: u64) -> Vec<u64> {
    if n == 0 { return vec![]; }
    // Generate Fibonacci numbers up to n
    let mut fibs = vec![1u64, 2];
    while *fibs.last().unwrap() <= n {
        let next = fibs[fibs.len() - 1].checked_add(fibs[fibs.len() - 2]).unwrap_or(u64::MAX);
        if next > n { break; }
        fibs.push(next);
    }
    // Greedy: pick largest Fibonacci ≤ remainder
    let mut remainder = n;
    let mut result = Vec::new();
    for &f in fibs.iter().rev() {
        if f <= remainder {
            result.push(f);
            remainder -= f;
        }
        if remainder == 0 { break; }
    }
    result
}

/// Phyllotaxis spiral: position of the nth element in a Fibonacci/phi spiral arrangement
pub fn phyllotaxis_position(n: usize, scale: f64) -> (f64, f64) {
    let golden_angle = std::f64::consts::TAU / (1.0 + (5.0_f64).sqrt() / 2.0);
    let angle = n as f64 * golden_angle;
    let radius = scale * (n as f64).sqrt();
    (radius * angle.cos(), radius * angle.sin())
}

/// Verify phyllotaxis has uniform density (no large gaps)
pub fn phyllotaxis_uniformity(count: usize, bins: usize) -> f64 {
    let mut angle_bins = vec![0usize; bins];
    let golden_angle = std::f64::consts::TAU / (1.0 + (5.0_f64).sqrt() / 2.0);
    for i in 0..count {
        let angle = (i as f64 * golden_angle) % std::f64::consts::TAU;
        let bin = (angle / std::f64::consts::TAU * bins as f64) as usize % bins;
        angle_bins[bin] += 1;
    }
    let mean = count as f64 / bins as f64;
    let variance: f64 = angle_bins.iter()
        .map(|&c| (c as f64 - mean).powi(2))
        .sum::<f64>() / bins as f64;
    let std_dev = variance.sqrt();
    // Coefficient of variation (lower = more uniform)
    if mean > 0.0 { std_dev / mean } else { f64::NAN }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fibonacci_sequence_standard() {
        let fib = standard_fibonacci(10);
        assert_eq!(fib, vec![1, 1, 2, 3, 5, 8, 13, 21, 34, 55]);
    }

    #[test]
    fn fibonacci_custom_start() {
        let seq = fibonacci_sequence(2, 5, 6);
        assert_eq!(seq, vec![2, 5, 7, 12, 19, 31]);
    }

    #[test]
    fn golden_ratio_converges() {
        let ratio = golden_ratio_approx(20);
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((ratio - phi).abs() < 0.001, "Got {}, expected ~{}", ratio, phi);
    }

    #[test]
    fn penrose_inflate_ratio() {
        // Start with a single A edge, inflate many times
        let inflated = penrose_inflate(&[1], 15);
        let ratio = edge_ratio(&inflated);
        let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
        assert!((ratio - phi).abs() < 0.01, "Penrose ratio {} vs φ={}", ratio, phi);
    }

    #[test]
    fn mandelbrough_roughness_sane() {
        let dim = mandelbrot_roughness(100, 2000);
        // Mandelbrot boundary dimension is ~1.5; grid-based estimate is rough
        assert!(dim > 0.5 && dim < 2.5, "Got dimension {}", dim);
    }

    #[test]
    fn fibonacci_cr_profile_runs() {
        // Build a small Fibonacci-structured graph
        let adj = vec![
            vec![0.0, 1.0, 1.0, 0.0],
            vec![1.0, 0.0, 0.0, 1.0],
            vec![1.0, 0.0, 0.0, 1.0],
            vec![0.0, 1.0, 1.0, 0.0],
        ];
        let profile = fibonacci_cr_profile(&adj, 2);
        assert!(!profile.is_empty(), "Should have at least one CR value");
        // All CR values should be in [0, 1]
        for &cr in &profile {
            assert!(cr >= 0.0 && cr <= 1.0, "CR {} out of range", cr);
        }
    }

    #[test]
    fn zeckendorf_representation() {
        let rep = zeckendorf(100);
        // 100 = 89 + 8 + 3
        assert_eq!(rep, vec![89, 8, 3]);
        // Verify no consecutive Fibonacci numbers
        for w in rep.windows(2) {
            let fibs = standard_fibonacci(30);
            let idx_a = fibs.iter().position(|&f| f == w[0]).unwrap();
            let idx_b = fibs.iter().position(|&f| f == w[1]).unwrap();
            assert!((idx_a as i64 - idx_b as i64).abs() > 1,
                "Zeckendorf violated: {} and {} are consecutive Fibonacci", w[0], w[1]);
        }
    }

    #[test]
    fn phyllotaxis_positions_spread() {
        let positions: Vec<(f64, f64)> = (0..100).map(|i| phyllotaxis_position(i, 1.0)).collect();
        // All positions should be unique
        for i in 0..positions.len() {
            for j in (i+1)..positions.len() {
                let dx = positions[i].0 - positions[j].0;
                let dy = positions[i].1 - positions[j].1;
                assert!(dx * dx + dy * dy > 1e-10, "Duplicate positions {} and {}", i, j);
            }
        }
    }

    #[test]
    fn phyllotaxis_is_uniform() {
        let cv = phyllotaxis_uniformity(500, 36);
        // Coefficient of variation should be low (< 0.2 for good uniformity)
        assert!(cv < 0.2, "CV = {}, not uniform enough", cv);
    }
}
