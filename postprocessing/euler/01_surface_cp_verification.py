import math
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"

solution_file = data_dir / "weno5_solution.csv"
forces_file = data_dir / "weno5_forces.csv"

figure_dir.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Freestream and geometry
# ---------------------------------------------------------------------

gamma = 1.4
R = 287.0

M_inf = 2.5
p_inf = 101325.0
T_inf = 288.15
alpha_deg = 5.0

chord = 1.0
t_over_c = 0.10

alpha = math.radians(alpha_deg)

rho_inf = p_inf / (R * T_inf)
a_inf = math.sqrt(gamma * R * T_inf)
U_inf = M_inf * a_inf
q_inf = 0.5 * rho_inf * U_inf**2

theta = math.atan(t_over_c)
turn_upper_front = theta - alpha
turn_lower_front = theta + alpha
turn_rear = 2.0 * theta


# ---------------------------------------------------------------------
# Ackeret linearised supersonic theory
# ---------------------------------------------------------------------

theta_ackeret = t_over_c
ackeret_denominator = math.sqrt(M_inf**2 - 1.0)

Cp_UF_ackeret = 2.0 * (theta_ackeret - alpha) / ackeret_denominator
Cp_UR_ackeret = -2.0 * (theta_ackeret + alpha) / ackeret_denominator
Cp_LF_ackeret = 2.0 * (theta_ackeret + alpha) / ackeret_denominator
Cp_LR_ackeret = -2.0 * (theta_ackeret - alpha) / ackeret_denominator

Cl_ackeret = 4.0 * alpha / ackeret_denominator
Cd_ackeret = 4.0 * (alpha**2 + theta_ackeret**2) / ackeret_denominator


# ---------------------------------------------------------------------
# Nonlinear shock-expansion theory
# ---------------------------------------------------------------------

def theta_beta_m_residual(beta, M1, theta_turn, gamma):
    numerator = M1**2 * math.sin(beta)**2 - 1.0
    denominator = M1**2 * (gamma + math.cos(2.0 * beta)) + 2.0
    rhs = (2.0 / math.tan(beta)) * numerator / denominator

    return rhs - math.tan(theta_turn)


def bisection_root(function, left, right, tolerance=1.0e-12, max_iterations=100):
    f_left = function(left)
    f_right = function(right)

    if f_left * f_right > 0.0:
        raise ValueError("Bisection interval does not bracket a root.")

    for _ in range(max_iterations):
        middle = 0.5 * (left + right)
        f_middle = function(middle)

        if abs(f_middle) < tolerance:
            return middle

        if f_left * f_middle < 0.0:
            right = middle
        else:
            left = middle
            f_left = f_middle

    return 0.5 * (left + right)


def prandtl_meyer(M, gamma):
    first_term = math.sqrt((gamma + 1.0) / (gamma - 1.0))
    second_term = math.atan(
        math.sqrt(
            ((gamma - 1.0) / (gamma + 1.0))
            * (M**2 - 1.0)
        )
    )
    third_term = math.atan(math.sqrt(M**2 - 1.0))

    return first_term * second_term - third_term


def d_prandtl_meyer_dM(M, gamma):
    return (
        math.sqrt(M**2 - 1.0)
        / (M * (1.0 + 0.5 * (gamma - 1.0) * M**2))
    )


def expansion_mach(M1, turn_angle, gamma):
    target_nu = prandtl_meyer(M1, gamma) + turn_angle
    M2 = M1 + 0.5

    for _ in range(50):
        residual = prandtl_meyer(M2, gamma) - target_nu
        correction = residual / d_prandtl_meyer_dM(M2, gamma)
        M2 -= correction

        if abs(correction) < 1.0e-12:
            break

    return M2


# Upper-front oblique shock

M1 = M_inf

beta_upper = bisection_root(
    lambda beta: theta_beta_m_residual(
        beta, M1, turn_upper_front, gamma
    ),
    math.radians(24.0),
    math.radians(25.0)
)

Mn1_upper = M1 * math.sin(beta_upper)

pressure_ratio_upper = (
    1.0
    + (2.0 * gamma / (gamma + 1.0))
    * (Mn1_upper**2 - 1.0)
)

p_upper_front = p_inf * pressure_ratio_upper
Cp_UF_nonlinear = (p_upper_front - p_inf) / q_inf

Mn2_upper_squared = (
    ((gamma - 1.0) * Mn1_upper**2 + 2.0)
    / (2.0 * gamma * Mn1_upper**2 - (gamma - 1.0))
)

Mn2_upper = math.sqrt(Mn2_upper_squared)
M2_upper = Mn2_upper / math.sin(beta_upper - turn_upper_front)


# Lower-front oblique shock

beta_lower = bisection_root(
    lambda beta: theta_beta_m_residual(
        beta, M1, turn_lower_front, gamma
    ),
    math.radians(32.0),
    math.radians(33.0)
)

Mn1_lower = M1 * math.sin(beta_lower)

pressure_ratio_lower = (
    1.0
    + (2.0 * gamma / (gamma + 1.0))
    * (Mn1_lower**2 - 1.0)
)

p_lower_front = p_inf * pressure_ratio_lower
Cp_LF_nonlinear = (p_lower_front - p_inf) / q_inf

Mn2_lower_squared = (
    ((gamma - 1.0) * Mn1_lower**2 + 2.0)
    / (2.0 * gamma * Mn1_lower**2 - (gamma - 1.0))
)

Mn2_lower = math.sqrt(Mn2_lower_squared)
M2_lower = Mn2_lower / math.sin(beta_lower - turn_lower_front)


# Rear-panel Prandtl-Meyer expansions

M_upper_rear = expansion_mach(M2_upper, turn_rear, gamma)

pressure_ratio_upper_rear = (
    (1.0 + 0.5 * (gamma - 1.0) * M2_upper**2)
    / (1.0 + 0.5 * (gamma - 1.0) * M_upper_rear**2)
) ** (gamma / (gamma - 1.0))

p_upper_rear = p_upper_front * pressure_ratio_upper_rear
Cp_UR_nonlinear = (p_upper_rear - p_inf) / q_inf

M_lower_rear = expansion_mach(M2_lower, turn_rear, gamma)

pressure_ratio_lower_rear = (
    (1.0 + 0.5 * (gamma - 1.0) * M2_lower**2)
    / (1.0 + 0.5 * (gamma - 1.0) * M_lower_rear**2)
) ** (gamma / (gamma - 1.0))

p_lower_rear = p_lower_front * pressure_ratio_lower_rear
Cp_LR_nonlinear = (p_lower_rear - p_inf) / q_inf


# Integrated nonlinear aerodynamic coefficients

Cx_nonlinear = (
    0.05 * Cp_UF_nonlinear
    - 0.05 * Cp_UR_nonlinear
    + 0.05 * Cp_LF_nonlinear
    - 0.05 * Cp_LR_nonlinear
)

Cy_nonlinear = (
    -0.5 * Cp_UF_nonlinear
    - 0.5 * Cp_UR_nonlinear
    + 0.5 * Cp_LF_nonlinear
    + 0.5 * Cp_LR_nonlinear
)

Cd_nonlinear = (
    Cx_nonlinear * math.cos(alpha)
    + Cy_nonlinear * math.sin(alpha)
)

Cl_nonlinear = (
    -Cx_nonlinear * math.sin(alpha)
    + Cy_nonlinear * math.cos(alpha)
)


# ---------------------------------------------------------------------
# WENO5/HLLC data
# ---------------------------------------------------------------------

forces = pd.read_csv(forces_file)
solution = pd.read_csv(solution_file)

Cd_weno = float(forces["Cd"].iloc[-1])
Cl_weno = float(forces["Cl"].iloc[-1])

Cd_error = 100.0 * (Cd_weno - Cd_nonlinear) / Cd_nonlinear
Cl_error = 100.0 * (Cl_weno - Cl_nonlinear) / Cl_nonlinear


# ---------------------------------------------------------------------
# Console summary
# ---------------------------------------------------------------------

print("\nFreestream")
print("rho_inf =", rho_inf, "kg/m^3")
print("a_inf   =", a_inf, "m/s")
print("U_inf   =", U_inf, "m/s")
print("q_inf   =", q_inf, "Pa")

print("\nGeometry")
print("theta =", math.degrees(theta), "deg")
print("upper LE turn =", math.degrees(turn_upper_front), "deg")
print("lower LE turn =", math.degrees(turn_lower_front), "deg")
print("rear shoulder turn =", math.degrees(turn_rear), "deg")

print("\nShock angles")
print("upper beta =", math.degrees(beta_upper), "deg")
print("lower beta =", math.degrees(beta_lower), "deg")

print("\nNonlinear panel Cp")
print("upper front =", Cp_UF_nonlinear)
print("upper rear  =", Cp_UR_nonlinear)
print("lower front =", Cp_LF_nonlinear)
print("lower rear  =", Cp_LR_nonlinear)

print("\nAerodynamic coefficients")
print("Ackeret              : Cl =", Cl_ackeret, "Cd =", Cd_ackeret)
print("Nonlinear theory     : Cl =", Cl_nonlinear, "Cd =", Cd_nonlinear)
print("WENO5/HLLC           : Cl =", Cl_weno, "Cd =", Cd_weno)
print("WENO vs nonlinear    : dCl =", Cl_error, "%, dCd =", Cd_error, "%")


# ---------------------------------------------------------------------
# Surface Cp extraction
# ---------------------------------------------------------------------

fluid = solution[solution["is_solid"] == 0]
solid = solution[solution["is_solid"] == 1]

upper_surface = []
lower_surface = []

for x in sorted(solid["x"].unique()):
    if x < 0.0 or x > chord:
        continue

    solid_column = solid[solid["x"] == x]

    y_top = solid_column["y"].max()
    y_bottom = solid_column["y"].min()

    upper_candidates = fluid[
        (fluid["x"] == x)
        & (fluid["y"] > y_top)
    ]

    lower_candidates = fluid[
        (fluid["x"] == x)
        & (fluid["y"] < y_bottom)
    ]

    if not upper_candidates.empty:
        upper_cell = upper_candidates.loc[
            (upper_candidates["y"] - y_top).idxmin()
        ]
        upper_surface.append((x / chord, upper_cell["Cp"]))

    if not lower_candidates.empty:
        lower_cell = lower_candidates.loc[
            (y_bottom - lower_candidates["y"]).idxmin()
        ]
        lower_surface.append((x / chord, lower_cell["Cp"]))

upper_surface = pd.DataFrame(upper_surface, columns=["x_over_c", "Cp"])
lower_surface = pd.DataFrame(lower_surface, columns=["x_over_c", "Cp"])

upper_plot = upper_surface[
    (upper_surface["x_over_c"] > 0.04)
    & (upper_surface["x_over_c"] < 0.98)
]

lower_plot = lower_surface[
    (lower_surface["x_over_c"] > 0.04)
    & (lower_surface["x_over_c"] < 0.98)
]


# ---------------------------------------------------------------------
# Figure 1: aerodynamic coefficient comparison
# ---------------------------------------------------------------------

coefficients = ["$C_L$", "$C_D$"]

ackeret_values = [Cl_ackeret, Cd_ackeret]
nonlinear_values = [Cl_nonlinear, Cd_nonlinear]
weno_values = [Cl_weno, Cd_weno]

x = [0, 1]
width = 0.24

fig1, ax1 = plt.subplots(figsize=(8.5, 6.0))

bars_ackeret = ax1.bar(
    [i - width for i in x],
    ackeret_values,
    width,
    label="Ackeret"
)

bars_nonlinear = ax1.bar(
    x,
    nonlinear_values,
    width,
    label="Nonlinear shock-expansion"
)

bars_weno = ax1.bar(
    [i + width for i in x],
    weno_values,
    width,
    label="WENO5/HLLC"
)

for bars in (bars_ackeret, bars_nonlinear, bars_weno):
    ax1.bar_label(bars, fmt="%.5f", padding=3, fontsize=9)

ax1.set_xticks(x)
ax1.set_xticklabels(coefficients)
ax1.set_ylabel("Aerodynamic coefficient")
ax1.set_title(
    "Aerodynamic Coefficient Comparison\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)
ax1.set_ylim(0.0, 0.175)
ax1.grid(axis="y", alpha=0.3)
ax1.set_axisbelow(True)
ax1.legend()

fig1.tight_layout()

force_figure = figure_dir / "force_comparison.png"
fig1.savefig(force_figure, dpi=300, bbox_inches="tight")


# ---------------------------------------------------------------------
# Figure 2: surface pressure coefficient
# ---------------------------------------------------------------------

fig2, ax2 = plt.subplots(figsize=(9.5, 6.0))

ax2.plot(
    upper_plot["x_over_c"],
    upper_plot["Cp"],
    linewidth=1.8,
    label="WENO5/HLLC upper"
)

ax2.plot(
    lower_plot["x_over_c"],
    lower_plot["Cp"],
    linewidth=1.8,
    label="WENO5/HLLC lower"
)

ax2.step(
    [0.0, 0.5, 1.0],
    [Cp_UF_nonlinear, Cp_UR_nonlinear, Cp_UR_nonlinear],
    where="post",
    linestyle="--",
    linewidth=1.6,
    label="Shock-expansion upper"
)

ax2.step(
    [0.0, 0.5, 1.0],
    [Cp_LF_nonlinear, Cp_LR_nonlinear, Cp_LR_nonlinear],
    where="post",
    linestyle="--",
    linewidth=1.6,
    label="Shock-expansion lower"
)

ax2.set_xlim(0.0, 1.0)
ax2.set_xlabel("$x/c$")
ax2.set_ylabel("$C_p$")
ax2.set_title(
    "Surface Pressure Coefficient\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)
ax2.grid(True, alpha=0.3)
ax2.set_axisbelow(True)
ax2.legend()
ax2.invert_yaxis()

fig2.tight_layout()

surface_figure = figure_dir / "surface_cp.png"
fig2.savefig(surface_figure, dpi=300, bbox_inches="tight")


print("\nSaved figures")
print(force_figure)
print(surface_figure)

plt.show()

