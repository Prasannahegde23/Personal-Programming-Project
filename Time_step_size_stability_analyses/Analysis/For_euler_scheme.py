#!/usr/bin/env python3
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 1) Critical timestep (in seconds)
dt_crit = 6.22e-9

# 2) Reference time step size
alpha_ref = 0.125

def mean_conc(fn: Path) -> float:
    """
    Read the file `fn`, skip the 3-line header + load-step line,
    and return the mean of the 'H_conc' column.
    """
    df = pd.read_csv(
        fn,
        skiprows=4,
        delim_whitespace=True,
        names=["x","y","disp","H_conc","damage"]
    )
    return df["H_conc"].mean()

def parse_alpha(fn: Path) -> float:
    """
    From filename stems:
      'Stable'         -> alpha = 1.0
      'critic_0_5'     -> alpha = 0.5
      'critic_0_25'    -> alpha = 0.25
      'critic_0_125'   -> alpha = 0.125
    """
    s = fn.stem.lower()
    if "stable" in s:
        return 1.0
    m = re.search(r"critic_(\d+)_?(\d*)", s)
    if not m:
        raise RuntimeError(f"Cannot parse alpha from '{fn.name}'")
    a, b = m.groups()
    return float(f"{a}.{b}") if b else float(a)

def main():
    # 3) List your files here (ensure they exist in this directory)
    files = [
        Path("Stable.txt"),
        Path("critic_0_5.txt"),
        Path("critic_0_25.txt"),
        Path("critic_0_125.txt"),
    ]

    # 4) Read mean values for each file
    data = []
    for fn in files:
        if not fn.exists():
            raise FileNotFoundError(f"File not found: {fn}")
        alpha = parse_alpha(fn)
        dt = alpha * dt_crit
        C  = mean_conc(fn)
        data.append((alpha, dt, C))

    # 5) Extract reference concentration at alpha_ref
    try:
        C_ref = next(C for alpha, dt, C in data if np.isclose(alpha, alpha_ref))
    except StopIteration:
        raise RuntimeError(f"No file found for reference alpha = {alpha_ref}")

    # 6) Build (dt, error) pairs for all non-reference alphas
    pairs = []
    for alpha, dt, C in data:
        if not np.isclose(alpha, alpha_ref):
            err = abs(C - C_ref)
            pairs.append((dt, err))

    # 7) Sort by dt and unpack
    pairs.sort(key=lambda x: x[0])
    dts_plot, errors = zip(*pairs)

    # 8) Fit ln(error) = p ln(dt) + const → get slope p
    p_est, _ = np.polyfit(np.log(dts_plot), np.log(errors), 1)

    # 9) Plot convergence
    plt.figure(figsize=(6,4))
    plt.loglog(dts_plot, errors, "o-")
    for x, e in zip(dts_plot, errors):
        plt.text(x, e, f"{e:.1e}", ha="center", va="bottom", fontsize=8)

    plt.xlabel(r"$\Delta t$ (s)", fontsize=12)
    plt.ylabel(r"Error $|C - C_{\mathrm{ref}}|$ (mol/$m^2)$", fontsize=12)
    plt.title(f"Convergence: $e\\sim(\\Delta t)^p$  (p ≃ {p_est:.2f})", fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

