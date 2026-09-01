from pathlib import Path

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Paths and flow conditions
# ------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"
solution_file = data_dir / "weno5_solution.csv"

figure_dir.mkdir(parents=True, exist_ok=True)

gamma = 1.4
M_inf = 2.5
p_inf = 101325.0
alpha_deg = 5.0

alpha = math.radians(alpha_deg)
theta = math.atan(0.10)

lower_turn = theta + alpha

q_inf = 0.5 * gamma * p_inf * M_inf**2


# ------------------------------------------------------------
# Theta-beta-M relation
# ------------------------------------------------------------

def theta_beta_residual(beta):
    numerator = M_inf**2 * math.sin(beta)**2 - 1.0
    denominator = M_inf**2 * (gamma + math.cos(2.0 * beta)) + 2.0

    rhs = (
        2.0 / math.tan(beta)
        * numerator / denominator
    )

    return rhs - math.tan(lower_turn)


def bisection(left, right, tolerance=1.0e-12):
    f_left = theta_beta_residual(left)

    for _ in range(100):
        middle = 0.5 * (left + right)
        f_middle = theta_beta_residual(middle)

        if abs(f_middle) < tolerance:
            return middle

        if f_left * f_middle < 0.0:
            right = middle
        else:
            left = middle
            f_left = f_middle

    return 0.5 * (left + right)


beta = bisection(
    math.radians(32.0),
    math.radians(33.0)
)


# ------------------------------------------------------------
# Exact post-shock state
# ------------------------------------------------------------

Mn1 = M_inf * math.sin(beta)

pressure_ratio = (
    1.0
    + (2.0 * gamma / (gamma + 1.0))
    * (Mn1**2 - 1.0)
)

Mn2_squared = (
    ((gamma - 1.0) * Mn1**2 + 2.0)
    / (2.0 * gamma * Mn1**2 - (gamma - 1.0))
)

Mn2 = math.sqrt(Mn2_squared)
M2 = Mn2 / math.sin(beta - lower_turn)

Cp2 = (
    (pressure_ratio - 1.0)
    * p_inf / q_inf
)


# ------------------------------------------------------------
# Load CFD solution
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
# Horizontal cut through lower shock
# ------------------------------------------------------------

y_target = -0.30

j = np.argmin(np.abs(y - y_target))
y_cut = y[j]

mach_cut = mach[j, :]
pressure_cut = pressure[j, :]
cp_cut = cp[j, :]


# ------------------------------------------------------------
# Analytical shock location
# ------------------------------------------------------------

# Lower shock angle relative to body x-axis
shock_angle_body = beta - alpha

x_shock_theory = (
    abs(y_cut)
    / math.tan(shock_angle_body)
)


# ------------------------------------------------------------
# Numerical shock location from maximum pressure gradient
# ------------------------------------------------------------

dp_dx = np.gradient(pressure_cut, x)

search = (
    (x > 0.25)
    & (x < 0.80)
)

search_indices = np.where(search)[0]

shock_index = search_indices[
    np.argmax(dp_dx[search])
]

x_shock_cfd = x[shock_index]

shock_angle_body_cfd = math.atan(
    abs(y_cut) / x_shock_cfd
)

beta_cfd = shock_angle_body_cfd + alpha

beta_error = (
    math.degrees(beta_cfd)
    - math.degrees(beta)
)


# ------------------------------------------------------------
# Piecewise analytical states
# ------------------------------------------------------------

x_plot_min = 0.20
x_plot_max = 0.82

theory_mach = np.where(
    x < x_shock_theory,
    M_inf,
    M2
)

theory_pressure = np.where(
    x < x_shock_theory,
    1.0,
    pressure_ratio
)

theory_cp = np.where(
    x < x_shock_theory,
    0.0,
    Cp2
)


# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

print("\nLower oblique-shock verification")
print("cut y/c =", y_cut)
print()

print("Analytical shock")
print("beta =", math.degrees(beta), "deg")
print("shock location x/c =", x_shock_theory)
print("post-shock Mach =", M2)
print("post-shock p/p_inf =", pressure_ratio)
print("post-shock Cp =", Cp2)

print("\nNumerical shock crossing")
print("shock location x/c =", x_shock_cfd)

print("\nShock location difference")
print("delta x/c =", x_shock_cfd - x_shock_theory)

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, axes = plt.subplots(
    1,
    3,
    figsize=(13, 4.3),
    sharex=True
)


# Mach
axes[0].plot(
    x,
    mach_cut,
    linewidth=1.8,
    label="WENO5/HLLC"
)

axes[0].plot(
    x,
    theory_mach,
    "k--",
    linewidth=1.3,
    label="Oblique-shock theory"
)

axes[0].set_ylabel("Mach number")
axes[0].set_title("Mach")


# Static pressure
axes[1].plot(
    x,
    pressure_cut,
    linewidth=1.8
)

axes[1].plot(
    x,
    theory_pressure,
    "k--",
    linewidth=1.3
)

axes[1].set_ylabel(r"$p/p_\infty$")
axes[1].set_title("Static pressure")


# Pressure coefficient
axes[2].plot(
    x,
    cp_cut,
    linewidth=1.8
)

axes[2].plot(
    x,
    theory_cp,
    "k--",
    linewidth=1.3
)

axes[2].set_ylabel(r"$C_p$")
axes[2].set_title("Pressure coefficient")


# Common formatting
for ax in axes:
    ax.axvline(
        x_shock_theory,
        color="black",
        linestyle=":",
        linewidth=1.1
    )

    ax.axvline(
        x_shock_cfd,
        color="0.45",
        linestyle="--",
        linewidth=1.1
    )

    ax.set_xlim(
        x_plot_min,
        x_plot_max
    )

    ax.set_xlabel(r"$x/c$")
    ax.grid(True, alpha=0.25)
    ax.set_axisbelow(True)


axes[0].legend(
    fontsize=8,
    loc="best"
)

fig.suptitle(
    "Lower Leading-Edge Shock Verification\n"
    rf"$y/c={y_cut:.3f},\ "
    rf"M_\infty=2.5,\ "
    rf"\alpha=5^\circ$",
    fontsize=14
)

fig.tight_layout()


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

output_file = figure_dir / "lower_oblique_shock_profiles.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:", output_file)

plt.show()
