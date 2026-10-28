#!/usr/bin/env python3
"""
Input: a GROMACS .xvg file (e.g. from `gmx energy`), or any 2-column text file (time, value).
Run without a file to see a demo on synthetic AR(1) data with a known exact answer.
"""
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_xvg(path, col=1):
    """Read a GROMACS .xvg (skips # and @ lines) or plain whitespace text file."""
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith(("#", "@")) or not line.strip():
                continue
            rows.append([float(v) for v in line.split()])
    data = np.array(rows)
    return data[:, 0], data[:, col]


def block_average(x):
    """Flyvbjerg-Petersen blocking. Returns block sizes, SEM estimates and their errors."""
    x = np.asarray(x, dtype=float)
    sizes, sems, sem_errs = [], [], []
    block = 1
    while len(x) >= 4:
        n = len(x)
        sem = np.sqrt(x.var(ddof=1) / n)
        sizes.append(block)
        sems.append(sem)
        sem_errs.append(sem / np.sqrt(2 * (n - 1)))
        m = n // 2
        x = 0.5 * (x[0:2 * m:2] + x[1:2 * m:2])  # merge neighbouring pairs
        block *= 2
    return np.array(sizes), np.array(sems), np.array(sem_errs)


def plateau_sem(sizes, sems, errs):
    """First block level whose SEM agrees with the next level within error bars."""
    for i in range(1, len(sems) - 1):
        if abs(sems[i + 1] - sems[i]) < errs[i] + errs[i + 1] and sems[i] >= sems[:i].max():
            return sems[i], sizes[i]
    return sems[-1], sizes[-1]  # no plateau: trajectory probably too short


def autocorrelation(x):
    """Normalised autocorrelation function via FFT."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    n = len(x)
    f = np.fft.rfft(x, n=2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    return acf / acf[0]


def integrated_tau(acf, c=5.0):
    """Integrated autocorrelation time with Sokal's automatic window (window >= c*tau)."""
    tau = 0.5
    for w in range(1, len(acf)):
        tau += acf[w]
        if w >= c * tau:
            break
    return tau


def analyse(t, x, label, out_prefix):
    n = len(x)
    mean, std = x.mean(), x.std(ddof=1)
    naive = std / np.sqrt(n)

    sizes, sems, errs = block_average(x)
    sem_block, blk = plateau_sem(sizes, sems, errs)

    acf = autocorrelation(x)
    tau = integrated_tau(acf)
    g = 1 + 2 * tau
    sem_acf = np.sqrt(g) * naive
    dt = t[1] - t[0] if len(t) > 1 else 1.0

    print(f"\n=== {label} ===")
    print(f"N frames                    : {n}")
    print(f"Mean                        : {mean:.6g}")
    print(f"Naive SEM (sigma/sqrt N)    : {naive:.4g}   <- usually WRONG for MD")
    print(f"Block-averaging SEM         : {sem_block:.4g}   (plateau at block size {blk})")
    print(f"Autocorrelation SEM         : {sem_acf:.4g}")
    print(f"tau_int                     : {tau:.2f} frames = {tau * dt:.4g} time units")
    print(f"Statistical inefficiency g  : {g:.2f}")
    print(f"Independent samples N_eff   : {n / g:.0f}")
    print(f"Naive error underestimated by factor ~{sem_block / naive:.1f}")
    print(f"Report as                   : {mean:.4g} +/- {sem_block:.2g}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    ax[0].plot(t, x, lw=0.5)
    ax[0].axhline(mean, color="k", ls="--")
    ax[0].set(xlabel="time", ylabel=label, title="Time series")

    ax[1].errorbar(np.log2(sizes), sems, yerr=errs, fmt="o-", capsize=3)
    ax[1].axhline(naive, color="r", ls=":", label="naive SEM")
    ax[1].axhline(sem_block, color="g", ls="--", label="plateau SEM")
    ax[1].set(xlabel="log2(block size)", ylabel="SEM", title="Flyvbjerg-Petersen blocking")
    ax[1].legend()

    lag = min(len(acf), int(20 * tau) + 10)
    ax[2].plot(np.arange(lag) * dt, acf[:lag])
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set(xlabel="lag time", ylabel="C(t)", title=f"Autocorrelation (tau_int = {tau:.1f} frames)")

    fig.tight_layout()
    fig.savefig(f"{out_prefix}.png", dpi=150)
    print(f"Plot saved                  : {out_prefix}.png")


def demo():
    """AR(1): x_t = phi*x_{t-1} + noise, exact tau_int = 0.5*(1+phi)/(1-phi)."""
    rng = np.random.default_rng(42)
    phi, n = 0.95, 100_000
    eps = rng.normal(size=n)
    x = np.empty(n)
    x[0] = 0.0
    for i in range(1, n):
        x[i] = phi * x[i - 1] + eps[i]
    exact_tau = 0.5 * (1 + phi) / (1 - phi)
    exact_sem = np.sqrt((1 / (1 - phi**2)) * (1 + phi) / (1 - phi) / n)
    print(f"Demo: AR(1), phi={phi}. Exact tau_int = {exact_tau:.2f}, exact SEM = {exact_sem:.4g}")
    analyse(np.arange(n, dtype=float), x, "AR(1) demo", "demo")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("file", nargs="?", help=".xvg or 2-column text file (omit for demo)")
    p.add_argument("--col", type=int, default=1, help="data column (default 1)")
    p.add_argument("--skip", type=float, default=0.0, help="discard data before this time (equilibration)")
    p.add_argument("--label", default="observable")
    args = p.parse_args()

    if args.file is None:
        demo()
        sys.exit()
    t, x = read_xvg(args.file, args.col)
    keep = t >= args.skip
    analyse(t[keep], x[keep], args.label, args.file.rsplit(".", 1)[0] + "_analysis")
