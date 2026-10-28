# md-error-analysis
Statistical error analysis for molecular dynamics. MD frames are correlated, so naive σ/√N underestimates uncertainty (often 5–10×). Computes the true error via Flyvbjerg–Petersen block averaging and integrated autocorrelation time, with statistical inefficiency and effective sample size. Reads GROMACS .xvg; validated on an exact AR(1) case.
