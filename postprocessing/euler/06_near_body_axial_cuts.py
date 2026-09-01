from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Paths and constants
# ------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"
solution_file = data_dir / "weno5_solution.csv"

figure_dir.mkdir(parents=True, exist_ok=True)

p_inf = 101325.0
chord = 1.0


# ------------------------------------------------------------
# Load solution
# ------------------------------------------------------------

solution = pd.read_csv(solution_file)

x = np.sort(solution["x"].unique())
y = np.sort(solution["y"].unique())


def get_field(variable):
    return (
        solution
        .pivot(index="y", columns="x", values=variable)
        .sort_index()
        .sort_index(axis=1)
        .to_numpy()
    )


mach = get_field("Mach")
pressure = get_field("p") / p_inf
cp = get_field("Cp")


# ------------------------------------------------------------
# Extract horizontal cuts
# ------------------------------------------------------------

y_upper_target = 0.08
y_lower_target = -0.08

upper_index = np.argmin(np.abs(y - y_upper_target))
lower_index = np.argmin(np.abs(y - y_lower_target))

y_upper = y[upper_index]
y_lower = y[lower_index]

x_over_c = x / chord

mach_upper = mach[upper_index, :]
mach_lower = mach[lower_index, :]

pressure_upper = pressure[upper_index, :]
pressure_lower = pressure[lower_index, :]

cp_upper = cp[upper_index, :]
cp_lower = cp[lower_index, :]


print("\nNear-body axial cuts")
print("upper cut y/c =", y_upper / chord)
print("lower cut y/c =", y_lower / chord)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, axes = plt.subplots(
    1,
    3,
    figsize=(13, 4.2),
    sharex=True
)


# Mach
axes[0].plot(
    x_over_c,
    mach_upper,
    linewidth=1.8,
    label=rf"$y/c={y_upper / chord:+.3f}$"
)

axes[0].plot(
    x_over_c,
    mach_lower,
    linewidth=1.8,
    label=rf"$y/c={y_lower / chord:+.3f}$"
)

axes[0].axhline(
    2.5,
    color="0.35",
    linestyle="--",
    linewidth=1.0,
    label="Freestream"
)

axes[0].set_ylabel("Mach number")
axes[0].set_title("Mach")


# Pressure
axes[1].plot(
    x_over_c,
    pressure_upper,
    linewidth=1.8
)

axes[1].plot(
    x_over_c,
    pressure_lower,
    linewidth=1.8
)

axes[1].axhline(
    1.0,
    color="0.35",
    linestyle="--",
    linewidth=1.0
)

axes[1].set_ylabel(r"$p/p_\infty$")
axes[1].set_title("Static pressure")


# Pressure coefficient
axes[2].plot(
    x_over_c,
    cp_upper,
    linewidth=1.8
)

axes[2].plot(
    x_over_c,
    cp_lower,
    linewidth=1.8
)

axes[2].axhline(
    0.0,
    color="0.35",
    linestyle="--",
    linewidth=1.0
)

axes[2].set_ylabel(r"$C_p$")
axes[2].set_title("Pressure coefficient")


# ------------------------------------------------------------
# Common formatting
# ------------------------------------------------------------

for ax in axes:
    ax.axvspan(
        0.0,
        1.0,
        color="0.85",
        alpha=0.35
    )

    ax.axvline(
        0.0,
        color="0.55",
        linestyle=":",
        linewidth=0.8
    )

    ax.axvline(
        0.5,
        color="0.55",
        linestyle=":",
        linewidth=0.8
    )

    ax.axvline(
        1.0,
        color="0.55",
        linestyle=":",
        linewidth=0.8
    )

    ax.set_xlim(-0.25, 1.75)
    ax.set_xlabel(r"$x/c$")
    ax.grid(True, alpha=0.25)
    ax.set_axisbelow(True)


axes[0].legend(
    frameon=True,
    fontsize=8
)

fig.suptitle(
    "Near-Body Axial Flow Profiles\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$",
    fontsize=14
)

fig.tight_layout()


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

output_file = figure_dir / "near_body_axial_cuts.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("Saved:", output_file)

plt.show()
