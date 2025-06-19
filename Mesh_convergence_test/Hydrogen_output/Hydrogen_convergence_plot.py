#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  1 01:47:00 2025

@author: prasannahegde
"""

import numpy as np
import matplotlib.pyplot as plt

# === USER CONFIG ===
file_paths = [
    "mesh_4000.txt",
    "mesh_6000.txt",
    "mesh_8000.txt",
    "mesh_10000.txt",
    "mesh_12000.txt"
]
node_counts       = [4000, 6000, 8000, 10000, 12000]
load_step_target  = 3000

# thresholds for filtering
hyd_threshold  = 1e-7    # ignore nodes with c  < this

plot_fields = ['hyd']  
# ===================

# container
mean_hyd  = []

# --- data extraction ---
for fp in file_paths:
    with open(fp,'r') as f:
        # skip header
        for _ in range(3):
            next(f)
        lines = [ln.strip() for ln in f if ln.strip()]

    c_vals = []
    collecting = False

    for line in lines:
        if line.startswith("Load step"):
            try:
                step = int(line.split(":",1)[1])
            except ValueError:
                collecting = False
            else:
                collecting = (step == load_step_target)
            continue

        if not collecting:
            continue

        parts = line.split()
        # skip non-data rows
        try:
            float(parts[0])
        except:
            continue

        c_vals.append(float(parts[3]))

    # hydrogen: mean(c) of nodes ≥ hyd_threshold
    c = np.array(c_vals)
    mask_c = c >= hyd_threshold
    mean_hyd.append(c[mask_c].mean() if mask_c.any() else 0.0)

# --- prepare plot definitions ---
field_defs = {
    'hyd': {
        'data': mean_hyd,
        'ylabel': 'Mean Hydrogen Concentration (mol/m²)',
        'title': f'Hydrogen Convergence @ Load Step {load_step_target}',
        'ylim': (0.5e-5, 3.0e-5),
        'style': 'r*-',
        'label': 'Hydrogen Concentration'
    },
}

# --- plotting loop ---
for field in plot_fields:
    if field not in field_defs:
        raise ValueError(f"Unknown field '{field}', must be one of {list(field_defs)}")
    props = field_defs[field]

    plt.figure()
    plt.plot(node_counts, props['data'], props['style'], label=props['label'])
    plt.xlabel('Number of nodes')
    plt.ylabel(props['ylabel'])
    plt.title(props['title'])
    if props['ylim'] is not None:
        plt.ylim(*props['ylim'])
    plt.grid(True)
    plt.legend()

plt.show()