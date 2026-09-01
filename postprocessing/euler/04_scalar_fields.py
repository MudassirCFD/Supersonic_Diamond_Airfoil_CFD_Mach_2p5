from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm


# ---------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"
solution_file = data_dir / "weno5_solution.csv"

figure_dir.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------
# Common data loading
# ---------------------------------------------------------------------

def load_solution():
    solution = pd.read_csv(solution_file)

    x = np.sort(solution["x"].unique())
    y = np.sort(solution["y"].unique())

    solid = (
        solution
        .pivot(index="y", columns="x", values="is_solid")
        .sort_index()
        .sort_index(axis=1)
        .to_numpy()
    )

    X, Y = np.meshgrid(x, y)

    return solution, X, Y, solid


def get_field(solution, variable):
    field = (
        solution
        .pivot(index="y", columns="x", values=variable)
        .sort_index()
        .sort_index(axis=1)
        .to_numpy()
    )

    return field


# ---------------------------------------------------------------------
# Common plotting function
# ---------------------------------------------------------------------

def plot_field(
    X,
    Y,
    field,
    solid,
    title,
    colorbar_label,
    output_name,
    vmin,
    vmax,
    solid_value,
    cmap="cividis",
    center=None
):
    field = field.copy()

    # Values inside the immersed body are not physical.
    # Give them a harmless value before covering them with the exact airfoil.
    field[solid == 1] = solid_value

    levels = np.linspace(vmin, vmax, 120)

    fig, ax = plt.subplots(figsize=(12, 7))

    if center is None:
        contour = ax.contourf(
            X,
            Y,
            field,
            levels=levels,
            cmap=cmap,
            extend="both",
            antialiased=False
        )
    else:
        norm = TwoSlopeNorm(
            vmin=vmin,
            vcenter=center,
            vmax=vmax
        )

        contour = ax.contourf(
            X,
            Y,
            field,
            levels=levels,
            cmap=cmap,
            norm=norm,
            extend="both",
            antialiased=False
        )

    # Exact diamond geometry
    x_airfoil = [0.0, 0.5, 1.0, 0.5, 0.0]
    y_airfoil = [0.0, 0.05, 0.0, -0.05, 0.0]

    ax.fill(
        x_airfoil,
        y_airfoil,
        facecolor="0.20",
        edgecolor="black",
        linewidth=1.2,
        zorder=20
    )

    colorbar = fig.colorbar(
        contour,
        ax=ax,
        pad=0.02
    )

    colorbar.set_label(colorbar_label)

    ax.set_xlim(-1.0, 3.0)
    ax.set_ylim(-1.5, 1.5)

    ax.set_xlabel("$x/c$")
    ax.set_ylabel("$y/c$")

    ax.set_title(
        title + "\n" +
        r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
    )

    ax.set_aspect("equal", adjustable="box")

    fig.tight_layout()

    output_file = figure_dir / output_name

    fig.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    print("Saved:", output_file)

    plt.close(fig)


# ---------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------

def main():
    gamma = 1.4
    R = 287.0

    M_inf = 2.5
    p_inf = 101325.0
    T_inf = 288.15
    alpha_deg = 5.0

    alpha = np.radians(alpha_deg)

    rho_inf = p_inf / (R * T_inf)
    a_inf = np.sqrt(gamma * R * T_inf)
    U_inf = M_inf * a_inf

    solution, X, Y, solid = load_solution()
    fluid = (solid == 0)

    print("\nLoaded solution")
    print("Rows =", len(solution))
    print("Output folder =", figure_dir)

    # -----------------------------------------------------------------
    # Mach number
    # -----------------------------------------------------------------

    mach = get_field(solution, "Mach")

    print("\nMach field")
    print("Minimum =", solution.loc[solution["is_solid"] == 0, "Mach"].min())
    print("Maximum =", solution.loc[solution["is_solid"] == 0, "Mach"].max())

    plot_field(
        X,
        Y,
        mach,
        solid,
        title="Mach Number Distribution Around the Diamond Airfoil",
        colorbar_label="Mach number",
        output_name="mach.png",
        vmin=1.45,
        vmax=3.00,
        solid_value=M_inf
    )

    # -----------------------------------------------------------------
    # Static pressure
    # -----------------------------------------------------------------

    pressure = get_field(solution, "p") / p_inf

    print("\nPressure field")
    print("Minimum p/p_inf =", pressure[fluid].min())
    print("Maximum p/p_inf =", pressure[fluid].max())

    plot_field(
        X,
        Y,
        pressure,
        solid,
        title="Static Pressure Distribution Around the Diamond Airfoil",
        colorbar_label=r"$p/p_\infty$",
        output_name="pressure.png",
        vmin=pressure[fluid].min(),
        vmax=pressure[fluid].max(),
        solid_value=1.0
    )

    # -----------------------------------------------------------------
    # Density
    # -----------------------------------------------------------------

    density = get_field(solution, "rho") / rho_inf

    print("\nDensity field")
    print("Minimum rho/rho_inf =", density[fluid].min())
    print("Maximum rho/rho_inf =", density[fluid].max())

    plot_field(
        X,
        Y,
        density,
        solid,
        title="Density Distribution Around the Diamond Airfoil",
        colorbar_label=r"$\rho/\rho_\infty$",
        output_name="density.png",
        vmin=density[fluid].min(),
        vmax=density[fluid].max(),
        solid_value=1.0
    )

    # -----------------------------------------------------------------
    # Temperature
    # -----------------------------------------------------------------

    temperature = get_field(solution, "T") / T_inf

    print("\nTemperature field")
    print("Minimum T/T_inf =", temperature[fluid].min())
    print("Maximum T/T_inf =", temperature[fluid].max())

    plot_field(
        X,
        Y,
        temperature,
        solid,
        title="Temperature Distribution Around the Diamond Airfoil",
        colorbar_label=r"$T/T_\infty$",
        output_name="temperature.png",
        vmin=temperature[fluid].min(),
        vmax=temperature[fluid].max(),
        solid_value=1.0
    )

    # -----------------------------------------------------------------
    # Velocity magnitude and components
    # -----------------------------------------------------------------

    u = get_field(solution, "u")
    v = get_field(solution, "v")

    velocity = np.sqrt(u**2 + v**2) / U_inf
    u_ratio = u / U_inf
    v_ratio = v / U_inf

    print("\nVelocity magnitude field")
    print("Minimum |V|/U_inf =", velocity[fluid].min())
    print("Maximum |V|/U_inf =", velocity[fluid].max())

    plot_field(
        X,
        Y,
        velocity,
        solid,
        title="Velocity Magnitude Around the Diamond Airfoil",
        colorbar_label=r"$|\mathbf{V}|/U_\infty$",
        output_name="velocity.png",
        vmin=velocity[fluid].min(),
        vmax=velocity[fluid].max(),
        solid_value=1.0
    )

    print("\nStreamwise velocity field")
    print("Minimum u/U_inf =", u_ratio[fluid].min())
    print("Maximum u/U_inf =", u_ratio[fluid].max())

    plot_field(
        X,
        Y,
        u_ratio,
        solid,
        title="Streamwise Velocity Distribution Around the Diamond Airfoil",
        colorbar_label=r"$u/U_\infty$",
        output_name="u_velocity.png",
        vmin=u_ratio[fluid].min(),
        vmax=u_ratio[fluid].max(),
        solid_value=np.cos(alpha)
    )

    v_limit = np.max(np.abs(v_ratio[fluid]))

    print("\nVertical velocity field")
    print("Maximum |v/U_inf| =", v_limit)

    plot_field(
        X,
        Y,
        v_ratio,
        solid,
        title="Vertical Velocity Distribution Around the Diamond Airfoil",
        colorbar_label=r"$v/U_\infty$",
        output_name="v_velocity.png",
        vmin=-v_limit,
        vmax=v_limit,
        solid_value=np.sin(alpha),
        cmap="RdBu_r",
        center=0.0
    )

    # -----------------------------------------------------------------
    # Speed of sound
    # -----------------------------------------------------------------

    sound_speed = get_field(solution, "a") / a_inf

    print("\nSpeed of sound field")
    print("Minimum a/a_inf =", sound_speed[fluid].min())
    print("Maximum a/a_inf =", sound_speed[fluid].max())

    plot_field(
        X,
        Y,
        sound_speed,
        solid,
        title="Speed of Sound Distribution Around the Diamond Airfoil",
        colorbar_label=r"$a/a_\infty$",
        output_name="sound_speed.png",
        vmin=sound_speed[fluid].min(),
        vmax=sound_speed[fluid].max(),
        solid_value=1.0
    )

    # -----------------------------------------------------------------
    # Pressure coefficient
    # -----------------------------------------------------------------

    cp = get_field(solution, "Cp")
    cp_limit = np.max(np.abs(cp[fluid]))

    print("\nPressure coefficient field")
    print("Minimum Cp =", cp[fluid].min())
    print("Maximum Cp =", cp[fluid].max())

    plot_field(
        X,
        Y,
        cp,
        solid,
        title="Pressure Coefficient Field Around the Diamond Airfoil",
        colorbar_label=r"$C_p$",
        output_name="cp.png",
        vmin=-cp_limit,
        vmax=cp_limit,
        solid_value=0.0,
        cmap="RdBu_r",
        center=0.0
    )

    print("\nAll scalar-field figures completed.")
    print("Saved in:", figure_dir)


if __name__ == "__main__":
    main()
