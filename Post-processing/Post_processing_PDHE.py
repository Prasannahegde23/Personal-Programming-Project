import numpy as np
import matplotlib.pyplot as plt
import sys
import io
import matplotlib.animation as animation

def read_timesteps(filename):
    """Read the file, skip first 3 header lines, and split into (label, data‐array) for each load step."""
    # Read all lines, drop header
    with open(filename, 'r') as f:
        raw_lines = f.read().splitlines()
    data_lines = raw_lines[3:]  # skip the 3‐line header

    # Rebuild into one string and split on blank lines
    content = "\n".join(data_lines).strip()
    segments = content.split('\n\n')

    labels = []
    timesteps = []
    for seg in segments:
        lines = seg.strip().splitlines()
        if not lines:
            continue
        # first line is "Load step: NNNN"
        labels.append(lines[0].split(":", 1)[1].strip())
        # rest are numeric rows
        data = np.loadtxt(io.StringIO("\n".join(lines[1:])))
        timesteps.append(data)
    return labels, timesteps

def create_animation(labels, timesteps, col, col_name, out_filename):
    """
    Create and save an animation of column `col` vs. X/Y,
    with dynamic color scale based on the data range in that column.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    # compute vmin/vmax from the actual data column instead of fixed ±0.1
    all_vals = np.concatenate([ts[:, col] for ts in timesteps])
    vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)

    # initial scatter
    x0, y0, v0 = timesteps[0][:, 0], timesteps[0][:, 1], timesteps[0][:, col]
    scat = ax.scatter(
        x0, y0,
        c=v0,
        s=9,
        cmap='rainbow',
        marker='o',
        vmin=vmin,
        vmax=vmax
    )
    ax.set_xlabel('X in m')
    ax.set_ylabel('Y in m')
    ax.set_title(f"Load step: {labels[0]}")

    cbar = fig.colorbar(scat, ax=ax, label=col_name)

    def update(frame):
        data = timesteps[frame]
        x, y, v = data[:, 0], data[:, 1], data[:, col]
        scat.set_offsets(np.column_stack((x, y)))
        scat.set_array(v)
        scat.set_clim(vmin, vmax)
        ax.set_title(f"Load step: {labels[frame]}")
        return (scat,)

    ani = animation.FuncAnimation(
        fig, update,
        frames=len(timesteps),
        interval=500,
        blit=False
    )
    ani.save(out_filename, writer='pillow')
    plt.close(fig)


def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py <data_filename>")
        sys.exit(1)

    fname = sys.argv[1]
    labels, timesteps = read_timesteps(fname)

    # displacement is column 2
    create_animation(labels, timesteps, col=2, col_name='Displacement in m', out_filename='Displacement.gif')
    # hydrogen coverage is column 3
    create_animation(labels, timesteps, col=3, col_name='Hydrogen Coverage in $mol/m^2$', out_filename='Hydrogen_Coverage.gif')
    # damage is column 4
    create_animation(labels, timesteps, col=4, col_name='Damage', out_filename='Damage.gif')


if __name__ == '__main__':
    main()
