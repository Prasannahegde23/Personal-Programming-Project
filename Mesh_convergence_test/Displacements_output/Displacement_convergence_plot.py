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
disp_threshold = 1e-7    # ignore nodes with |u| < this
hyd_threshold  = 1e-7    # ignore nodes with c  < this

# pick which ones to plot: 'disp', 'hyd', 'dmg'
plot_fields = ['disp']  
# ===================

# container
mean_disp = []

# --- data extraction ---
for fp in file_paths:
    with open(fp,'r') as f:
        # skip header
        for _ in range(3):
            next(f)
        lines = [ln.strip() for ln in f if ln.strip()]

    u_vals = []

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

        u_vals.append(float(parts[2]))

    if not u_vals:
        raise RuntimeError(f"No data for load step {load_step_target} in {fp}")

    # displacement: mean(|u|) of nodes ≥ disp_threshold
    u = np.abs(u_vals)
    mask_u = u >= disp_threshold
    mean_disp.append(u[mask_u].mean() if mask_u.any() else 0.0)

# --- prepare plot definitions ---
field_defs = {
    'disp': {
        'data': mean_disp,
        'ylabel': 'Mean |Displacement| (m)',
        'title': f'Displacement Convergence @ Load Step {load_step_target}',
        'ylim': (0.9e-5, 1.1e-5),
        'style': 'o-',
        'label': 'Displacement'
    }
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