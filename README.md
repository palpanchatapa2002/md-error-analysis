# md-error-analysis

**Honest error bars for molecular dynamics time series.**

Consecutive MD frames are correlated, so the textbook standard error σ/√N underestimates the real uncertainty, often by a factor of 5–10. This small tool estimates the true statistical error of any MD observable (energy, density, area per lipid, box volume, RMSD…) in three ways:

| Method | What it gives |
|---|---|
| Flyvbjerg–Petersen block averaging | SEM from the plateau of the blocking curve |
| Integrated autocorrelation time τ_int | statistical inefficiency g = 1 + 2τ_int, SEM = σ√(g/N) |
| Effective sample size | N_eff = N / g independent samples |

## Theory in brief

For a correlated series x₁…x_N, the variance of the mean is

Var(x̄) = (σ²/N)·(1 + 2 Σₜ C(t)) = (σ²/N)·g

where C(t) is the normalised autocorrelation function. Block averaging gets the same quantity without computing C(t): repeatedly average neighbouring pairs of data points; once blocks are longer than the correlation time, the block means become independent and the estimated SEM stops growing (plateau).

## Validation

Running with no input generates an AR(1) process x_t = φx_{t−1} + ε_t (φ = 0.95, N = 10⁵), for which the answer is known exactly:

| Quantity | Exact | This code |
|---|---|---|
| τ_int (frames) | 19.50 | 20.7 |
| SEM | 0.0633 | 0.060 (blocking), 0.066 (ACF) |
| Naive SEM | — | 0.0102 (**~6× too small**) |

![demo](demo.png)

## Usage

```bash
pip install -r requirements.txt

# validation demo
python block_average.py

# GROMACS example: potential energy
echo "Potential" | gmx energy -f md.edr -o potential.xvg
python block_average.py potential.xvg --skip 1000 --label "Potential (kJ/mol)"
```

Options: `--col N` (data column, default 1), `--skip T` (discard equilibration before time T), `--label`.
Output: summary in the terminal and a 3-panel plot `<file>_analysis.png`.

## References

1. H. Flyvbjerg and H. G. Petersen, *J. Chem. Phys.* **91**, 461 (1989).
2. A. D. Sokal, *Monte Carlo Methods in Statistical Mechanics* (lecture notes, 1996).
3. A. Grossfield et al., "Best Practices for Quantification of Uncertainty and Sampling Quality in Molecular Simulations", *Living J. Comput. Mol. Sci.* **1**, 5067 (2018).
