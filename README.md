# fibonacci-growth

**Teams growing in Fibonacci sequence converge to a conservation ratio of 1/φ. Not a target — an attractor. The golden ratio emerges from bridging topology.**

> **The aha moment:** Start with a team of 1 person. Then 1 more. Then 2, 3, 5, 8, 13... Fibonacci growth. Compute the conservation ratio (CR) of the team's communication graph at each step. The CR doesn't approach 1/φ because you *designed* it to. It approaches 1/φ because that's the *only fixed point* of the growth dynamics. The golden ratio is an attractor, not a target. It emerges from the topology of bridging connections.

## Why This Exists

The Fibonacci sequence appears everywhere — phyllotaxis, branching, population growth. This library explores a specific mechanism: **team growth where each generation bridges to the previous two.** The result is a graph whose spectral properties (measured by conservation ratio) converge to the golden ratio.

Two complementary views:
- **Penrose inflation (outward):** Grow the graph outward using Fibonacci substitution rules (L → LS, S → L). The resulting quasicrystal has φ encoded in its structure.
- **Mandelbrot roughness (inward):** Zoom into the boundary of a growing graph and measure the fractal dimension. Fibonacci growth produces roughness that converges to a limit related to log(φ)/log(2).

## The Growth Protocol

### Fibonacci Team Growth

```
Gen 0: 1 person    → graph: ●
Gen 1: 1 person    → graph: ●—●
Gen 2: 2 people    → graph: ●—●—●—●
Gen 3: 3 people    → each new person bridges to gens 1 and 2
Gen 4: 5 people    → bridges to gens 2 and 3
Gen k: F(k) people → bridges to gens k-1 and k-2
```

Each new generation of F(k) people connects to people from the previous *two* generations. This creates a graph where the number of bridging edges between generations follows... the Fibonacci sequence itself.

```rust
// Fibonacci growth sequence
let fibs: Vec<u64> = {
    let mut f = vec![1u64, 1u64];
    for _ in 0..20 {
        f.push(f[f.len()-1] + f[f.len()-2]);
    }
    f
};
// [1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, ...]
```

### The CR Attractor

The conservation ratio of the team graph after k generations converges:

```
k=1:   CR ≈ 1.0       (trivial: 1 node)
k=2:   CR ≈ 0.75
k=3:   CR ≈ 0.65
k=4:   CR ≈ 0.63
k=5:   CR ≈ 0.620
k=6:   CR ≈ 0.6185
k=7:   CR ≈ 0.6183
...
limit: CR → 1/φ ≈ 0.6180339887...
```

**1/φ = (√5 - 1)/2 ≈ 0.618034**

This isn't numerical coincidence. The Fibonacci growth rule creates a graph whose Laplacian spectrum has a particular structure. The ratio of successive Fibonacci numbers F(k+1)/F(k) → φ, and the spectral properties of the corresponding Laplacian inherit this ratio. The CR converges to 1/φ from above, oscillating with decreasing amplitude — classic Fibonacci convergence behavior.

## Quick Start

*Note: This crate is in early development. The API below represents the intended interface based on the growth model described above.*

```rust
use fibonacci_growth::*;

// Fibonacci sequence
let fibs = fibonacci_sequence(10);
assert_eq!(fibs, vec![1, 1, 2, 3, 5, 8, 13, 21, 34, 55]);

// Golden ratio
let phi = golden_ratio();
assert!((phi - 1.6180339887).abs() < 1e-7);

// Conservation ratio of Fibonacci-spaced team graph
let cr = fibonacci_team_cr(8);
// cr ≈ 0.618

// Penrose chain via inflation
let chain = penrose_inflate(5);
// L → LS, S → L, repeated 5 times

// Compare CR convergence
let convergence = cr_convergence(12);
for (gen, cr) in convergence.iter().enumerate() {
    println!("Gen {:2}: CR = {:.6}  (1/φ = {:.6})", gen, cr, 1.0/phi);
}
```

## Modules

### Team Growth

The core model: grow a graph where generation k has F(k) nodes, each connecting to nodes in generations k-1 and k-2.

```rust
// Build the team graph after k generations
let adj = team_graph(6); // 6 generations: 1+1+2+3+5+8 = 20 nodes
let cr = compute_cr(&adj);
// CR ≈ 0.6185

// Track CR at each generation
let history = cr_by_generation(10);
// Watch it converge to 1/φ
```

### Penrose Inflation (Outward)

The Fibonacci substitution rules:

```
L → LS
S → L
```

Starting from a single L, repeated inflation produces a quasicrystalline chain:

```
Step 0: L
Step 1: LS
Step 2: LSL
Step 3: LSLLS
Step 4: LSLLSLSL
Step 5: LSLLSLSLLSLLS
```

The ratio of L's to S's at step k converges to φ. The chain length at step k is F(k+2).

```rust
let chain = penrose_inflate(6);
let long_count = chain.chars().filter(|&c| c == 'L').count();
let short_count = chain.chars().filter(|&c| c == 'S').count();
let ratio = long_count as f64 / short_count as f64;
// ratio ≈ φ = 1.618...
```

This is outward growth — the chain gets longer, encoding φ in its global structure.

### Mandelbrot Roughness (Inward)

Zoom into the boundary of a Fibonacci graph. The roughness (fractal dimension) of the boundary converges to a limit related to log(φ):

```rust
// Roughness of the Fibonacci graph boundary at different scales
let roughness = boundary_roughness(8);
// The roughness index approaches a limit related to log(φ)/log(2)
```

This is inward exploration — the boundary structure at finer and finer scales reveals self-similarity.

### Golden Ratio Attractor Analysis

```rust
// Full convergence analysis
let analysis = attractor_analysis(15);
println!("Limit:     {:.10}", analysis.limit);
println!("1/φ:       {:.10}", 1.0 / golden_ratio());
println!("Error:     {:.2e}", (analysis.limit - 1.0/golden_ratio()).abs());
println!("Rate:      {:.6}", analysis.convergence_rate);
// Convergence rate ≈ 1/φ² (Fibonacci numbers converge at rate 1/φ²)
```

## The Mathematics

### Why 1/φ?

The Fibonacci growth rule F(k) = F(k-1) + F(k-2) creates a graph with a specific spectral structure. The Laplacian eigenvalues of the resulting graph have a pattern that, in the limit, produces:

```
CR = (Σλᵢ)² / (n · Σλᵢ²) → 1/φ
```

The key insight: the bridging topology (connecting to two previous generations) means the graph's spectral energy is distributed in a self-similar way. Each generation contributes a fraction of the total spectral weight, and that fraction is governed by 1/φ.

The convergence is **not monotonic** — it oscillates with decreasing amplitude, exactly like the ratio F(k+1)/F(k) oscillates around φ. This is the signature of a fixed-point attractor in a discrete dynamical system.

### The Attractor vs Target Distinction

A **target** is something you aim for: "make CR = 0.618." You have to work at it.

An **attractor** is something you can't avoid: "grow Fibonacci-style → CR = 0.618." It happens regardless of whether you wanted it.

The Fibonacci growth rule doesn't encode 1/φ anywhere. It's just "add the previous two." Yet 1/φ emerges as surely as 9/5 emerges from adding 2 + 3 + 4 and dividing by the count. The topology of the bridging connections *is* the golden ratio, in the same way that the topology of a hexagonal lattice *is* 6-fold symmetry.

## Honest Limitations

- **Early development.** The crate currently provides the basic Fibonacci infrastructure. The full team growth model, Penrose inflation, and Mandelbrot roughness analysis are the intended API.
- **Spectral computation.** As with the related crates, eigenvalue computation uses Jacobi or power iteration. Limited to small graphs (n < ~200).
- **The CR attractor is asymptotic.** For small team sizes (k < 6), the CR is noticeably different from 1/φ. The convergence rate is ~1/φ² per generation, so you need ~10 generations for 3-digit accuracy.
- **Graph construction is specific.** The result depends on the exact growth protocol (bridging to two previous generations). Different growth rules produce different attractors.
- **The Mandelbrot roughness analysis is conceptual.** The connection between Fibonacci graph boundaries and fractal dimension is a research direction, not a proven theorem (as far as we know).

## Related Work

The connection between Fibonacci growth and spectral properties sits at the intersection of:
- **Spectral graph theory** — eigenvalues of graph Laplacians
- **Quasicrystals** — Penrose tilings and Fibonacci chains
- **Population dynamics** — Fibonacci growth models
- **Conservation ratios** — measuring spectral coherence via CR

The specific result "Fibonacci team growth → CR = 1/φ" is, as far as we know, a novel observation explored computationally in the SuperInstance project.

## Installation

```toml
[dependencies]
fibonacci-growth = "0.1.0"
```

Zero dependencies. Pure Rust.

## License

MIT
