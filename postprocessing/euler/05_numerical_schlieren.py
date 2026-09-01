from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"
solution_file = data_dir / "weno5_solution.csv"

figure_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Load solution and reshape to structured fields
# ------------------------------------------------------------

def load_solution():
    solution = pd.read_csv(solution_file)

    x = np.sort(solution["x"].unique())
    y = np.sort(solution["y"].unique())

    X, Y = np.meshgrid(x, y)

    solid = (
        solution
        .pivot(index="y", columns="x", values="is_solid")
        .sort_index()
        .sort_index(axis=1)
        .to_numpy()
    )

    rho = (
        solution
        .pivot(index="y", columns="x", values="rho")
        .sort_index()
        .sort_index(axis=1)
        .to_numpy()
    )

    return X, Y, solid, rho


# ------------------------------------------------------------
# Numerical schlieren from density gradient
# ------------------------------------------------------------

def compute_schlieren(X, Y, rho, solid):
    rho_plot = rho.copy()
    rho_plot[solid == 1] = np.nan

    x_line = X[0, :]
    y_line = Y[:, 0]

    grad_rho_y, grad_rho_x = np.gradient(rho_plot, y_line, x_line)
    grad_rho = np.sqrt(grad_rho_x**2 + grad_rho_y**2)

    fluid = np.isfinite(grad_rho) & (solid == 0)

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

    # Darken mid-level structures slightly for clearer display
    display_gamma = 0.7
    schlieren_display = schlieren**display_gamma

    return schlieren, schlieren_display, reference_gradient


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

def plot_schlieren(X, Y, schlieren_display):
    fig, ax = plt.subplots(figsize=(12, 7))

    image = ax.pcolormesh(
        X,
        Y,
        schlieren_display,
        cmap="Greys",
        shading="auto",
        vmin=0.0,
        vmax=1.0,
        rasterized=True
    )

    x_airfoil = [0.0, 0.5, 1.0, 0.5, 0.0]
    y_airfoil = [0.0, 0.05, 0.0, -0.05, 0.0]

    ax.fill(
        x_airfoil,
        y_airfoil,
        facecolor="white",
        edgecolor="black",
        linewidth=1.1,
        zorder=20
    )

    colorbar = fig.colorbar(
        image,
        ax=ax,
        pad=0.03
    )
    colorbar.set_label("Normalised schlieren strength")

    ax.set_xlim(-1.0, 3.0)
    ax.set_ylim(-1.5, 1.5)

    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$y/c$")

    ax.set_title(
        "Numerical Schlieren: WENO5-JS + HLLC\n"
        r"$M_\infty = 2.5,\ \alpha = 5^\circ$",
        fontweight="bold"
    )

    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()

    output_file = figure_dir / "schlieren.png"

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    print("Numerical schlieren")
    print("Reference density-gradient scale =", np.round(reference_gradient, 8))
    print("Saved:", output_file)

    plt.show()


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

if __name__ == "__main__":
    X, Y, solid, rho = load_solution()

    schlieren, schlieren_display, reference_gradient = compute_schlieren(
        X,
        Y,
        rho,
        solid
    )

    plot_schlieren(X, Y, schlieren_display)
