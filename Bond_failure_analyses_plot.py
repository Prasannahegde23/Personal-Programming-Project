#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

# Threshold below which damage is ignored
DAMAGE_THRESHOLD = 0.0

def read_damage_means(filename, thresh=DAMAGE_THRESHOLD):
    """
    Parses `filename`, skipping its first two header lines.
    Each frame starts with a line "Load step: N".
    Data lines must have ≥5 columns; column[4] is damage.
    Returns two lists: sorted load‐steps, and the corresponding
    mean damage (only including values ≥ thresh).
    """
    dmg_dict = defaultdict(list)

    with open(filename, 'r') as f:
        # skip headers
        next(f)
        next(f)

        current_step = None
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            # new frame?
            if ln.lower().startswith("load step:"):
                try:
                    current_step = int(ln.split(":", 1)[1])
                except ValueError:
                    current_step = None
                continue

            parts = ln.split()
            if current_step is None or len(parts) < 5:
                continue
            try:
                dmg = float(parts[4])
            except ValueError:
                continue

            if dmg >= thresh:
                dmg_dict[current_step].append(dmg)

    # sort and compute means (empty → 0.0)
    steps = sorted(dmg_dict.keys())
    means = [np.mean(dmg_dict[s]) if dmg_dict[s] else 0.0
             for s in steps]
    return steps, means

if __name__ == '__main__':
    files  = ['without_hydrogen.txt', 'with_hydrogen.txt']
    labels = ['Without hydrogen', 'With hydrogen']
    colors = ['magenta', 'blue']       # strong-contrast colors
    markers = ['o', '^']           # dot for first, triangle for second

    plt.figure(figsize=(8,5))
    for fn, lbl, col, mark in zip(files, labels, colors, markers):
        steps, means = read_damage_means(fn)
        plt.plot(
            steps, means,
            marker=mark,
            color=col,
            linestyle='-',
            linewidth=2,
            markersize=6,
            label=lbl
        )

    plt.xlabel('Load Steps')
    plt.ylabel('Mean Damage')
    plt.title('Mean Damage vs Load Step')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
