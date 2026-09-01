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
alpha_deg = 5.0
t_over_c = 0.10

alpha = math.radians(alpha_deg)
theta = math.atan(t_over_c)

upper_turn = theta - alpha
lower_turn = theta + alpha


# ------------------------------------------------------------
# Oblique-shock theory
# ------------------------------------------------------------

def theta_beta_residual(beta, turn):
    numerator = M_inf**2 * math.sin(beta)**2 - 1.0
    denominator = M_inf**2 * (gamma + math.cos(2.0 * beta)) + 2.0

    rhs = (
        2.0 / math.tan(beta)
        * numerator / denominator
    )

    return rhs - math.tan(turn)


def bisection(turn, left_deg, right_deg):
    left = math.radians(left_deg)
    right = math.radians(right_deg)

    f_left = theta_beta_residual(left, turn)

    for _ in range(100):
        middle = 0.5 * (left + right)
        f_middle = theta_beta_residual(middle, turn)

        if abs(f_middle) < 1.0e-12:
            return middle

        if f_left * f_middle < 0.0:
            right = middle
        else:
            left = middle
            f_left = f_middle

    return 0.5 * (left + right)


beta_upper = bisection(upper_turn, 24.0, 25.0)
beta_lower = bisection(lower_turn, 32.0, 33.0)

# Shock-ray angles relative to the airfoil x-axis
phi_upper_theory = alpha + beta_upper
phi_lower_theory = alpha - beta_lower


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


rho = get_field("rho")
solid = get_field("is_solid")

X, Y = np.meshgrid(x, y)
fluid = solid == 0


# ------------------------------------------------------------
# Density-gradient field
# ------------------------------------------------------------

rho_gradient = rho.copy()
rho_gradient[solid == 1] = np.nan

drho_dy, drho_dx = np.gradient(
    rho_gradient,
    y,
    x
)

grad_rho = np.sqrt(
    drho_dx**2
    + drho_dy**2
)


# ------------------------------------------------------------
# Schlieren display
# ------------------------------------------------------------

reference_gradient = np.nanpercentile(
    grad_rho[fluid],
    99.5
)

schlieren = grad_rho / reference_gradient
schlieren = np.clip(schlieren, 0.0, 1.0)

schlieren = np.nan_to_num(
    schlieren,
    nan=0.0,
    posinf=1.0,
    neginf=0.0
)

schlieren[solid == 1] = 0.0

# Same contrast convention as the final schlieren figure
schlieren_display = schlieren**0.7


# ------------------------------------------------------------
# Detect numerical shock ridges
# ------------------------------------------------------------

def detect_shock(angle, side):
    x_points = []
    y_points = []

    for i, x_pos in enumerate(x):

        if x_pos < 0.15 or x_pos > 1.20:
            continue

        y_expected = math.tan(angle) * x_pos
        band = np.abs(y - y_expected) < 0.10

        if side == "upper":
            band &= y > 0.0
        else:
            band &= y < 0.0

        indices = np.where(
            band & np.isfinite(grad_rho[:, i])
        )[0]

        if len(indices) == 0:
            continue

        j = indices[
            np.argmax(grad_rho[indices, i])
        ]

        x_points.append(x_pos)
        y_points.append(y[j])

    return np.array(x_points), np.array(y_points)


x_upper, y_upper = detect_shock(
    phi_upper_theory,
    "upper"
)

x_lower, y_lower = detect_shock(
    phi_lower_theory,
    "lower"
)


# ------------------------------------------------------------
# CFD shock-line fits
# ------------------------------------------------------------

m_upper, b_upper = np.polyfit(
    x_upper,
    y_upper,
    1
)

m_lower, b_lower = np.polyfit(
    x_lower,
    y_lower,
    1
)

phi_upper_cfd = math.atan(m_upper)
phi_lower_cfd = math.atan(m_lower)

beta_upper_cfd = abs(
    phi_upper_cfd - alpha
)

beta_lower_cfd = abs(
    phi_lower_cfd - alpha
)

upper_error = (
    math.degrees(beta_upper_cfd)
    - math.degrees(beta_upper)
)

lower_error = (
    math.degrees(beta_lower_cfd)
    - math.degrees(beta_lower)
)


# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

print("\nLeading-edge shock-angle verification")

print("\nUpper shock")
print("Theory beta =", math.degrees(beta_upper), "deg")
print("CFD beta    =", math.degrees(beta_upper_cfd), "deg")
print("Difference  =", upper_error, "deg")

print("\nLower shock")
print("Theory beta =", math.degrees(beta_lower), "deg")
print("CFD beta    =", math.degrees(beta_lower_cfd), "deg")
print("Difference  =", lower_error, "deg")


# ------------------------------------------------------------
# Lines for plotting
# ------------------------------------------------------------

# Analytical rays extend from the leading edge
x_theory = np.linspace(0.0, 3.0, 500)

upper_theory_line = (
    math.tan(phi_upper_theory) * x_theory
)

lower_theory_line = (
    math.tan(phi_lower_theory) * x_theory
)

# CFD fits are shown only over the region where data were sampled
x_upper_fit = np.linspace(
    x_upper.min(),
    x_upper.max(),
    250
)

x_lower_fit = np.linspace(
    x_lower.min(),
    x_lower.max(),
    250
)

upper_cfd_line = (
    m_upper * x_upper_fit
    + b_upper
)

lower_cfd_line = (
    m_lower * x_lower_fit
    + b_lower
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(12, 7))

ax.pcolormesh(
    X,
    Y,
    schlieren_display,
    cmap="Greys",
    shading="auto",
    vmin=0.0,
    vmax=1.0,
    rasterized=True
)


# Upper shock: blue family
ax.plot(
    x_theory,
    upper_theory_line,
    color="tab:blue",
    linestyle="--",
    linewidth=1.4,
    label="Upper theory"
)

ax.plot(
    x_upper_fit,
    upper_cfd_line,
    color="tab:blue",
    linewidth=1.7,
    label="Upper CFD fit"
)

ax.scatter(
    x_upper,
    y_upper,
    color="tab:blue",
    s=6,
    zorder=15
)


# Lower shock: orange family
ax.plot(
    x_theory,
    lower_theory_line,
    color="tab:orange",
    linestyle="--",
    linewidth=1.4,
    label="Lower theory"
)

ax.plot(
    x_lower_fit,
    lower_cfd_line,
    color="tab:orange",
    linewidth=1.7,
    label="Lower CFD fit"
)

ax.scatter(
    x_lower,
    y_lower,
    color="tab:orange",
    s=6,
    zorder=15
)


# ------------------------------------------------------------
# Exact diamond geometry
# ------------------------------------------------------------

x_airfoil = [0.0, 0.5, 1.0, 0.5, 0.0]
y_airfoil = [0.0, 0.05, 0.0, -0.05, 0.0]

ax.fill(
    x_airfoil,
    y_airfoil,
    facecolor="white",
    edgecolor="black",
    linewidth=1.2,
    zorder=20
)


# ------------------------------------------------------------
# Angle summary
# ------------------------------------------------------------

angle_text = (
    rf"Upper: $\beta_{{theory}}={math.degrees(beta_upper):.2f}^\circ$, "
    rf"$\beta_{{CFD}}={math.degrees(beta_upper_cfd):.2f}^\circ$, "
    rf"$\Delta\beta={upper_error:+.2f}^\circ$"
    "\n"
    rf"Lower: $\beta_{{theory}}={math.degrees(beta_lower):.2f}^\circ$, "
    rf"$\beta_{{CFD}}={math.degrees(beta_lower_cfd):.2f}^\circ$, "
    rf"$\Delta\beta={lower_error:+.2f}^\circ$"
)

ax.text(
    0.02,
    0.97,
    angle_text,
    transform=ax.transAxes,
    va="top",
    fontsize=9,
    bbox=dict(
        boxstyle="round,pad=0.35",
        facecolor="white",
        edgecolor="0.65",
        alpha=0.92
    )
)


# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

ax.set_xlim(-1.0, 3.0)
ax.set_ylim(-1.5, 1.5)

ax.set_xlabel("$x/c$")
ax.set_ylabel("$y/c$")

ax.set_title(
    "Leading-Edge Shock-Angle Verification\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)

ax.set_aspect(
    "equal",
    adjustable="box"
)

ax.legend(
    loc="lower left",
    fontsize=8,
    frameon=True,
    ncol=2
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fig.tight_layout()

output_file = (
    figure_dir
    / "schlieren_shock_angles.png"
)

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:", output_file)

plt.show()
