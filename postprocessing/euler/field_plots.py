"""Flow-field figures for the reference Euler solution."""

import math
import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm

from config import (
    AOA_RAD, CHORD, M_INF, ORANGE, T_OVER_C, U_INF, VY_INF,
    X_LE, Y_CENTER, draw_body, robust_limits, save_figure,
)
from theory import nonlinear_theory

def clean_flow_axes(ax, xs, ys):
    ax.set_aspect("equal")
    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ys[0], ys[-1])
    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$y/c$")
    ax.grid(False)

def add_theory_shock_lines(ax, nonlinear, x0=1.0, length=2.2, color="white"):
    # approximate visible shock rays from trailing edge region, useful as a guide only
    for key, sign in [("os_uf", +1.0), ("os_lf", -1.0)]:
        shock = nonlinear.get(key)
        if shock is None:
            continue

        beta = shock["beta"]
        slope = math.tan(beta - AOA_RAD)

        x = np.array([x0, x0 + length])
        y = sign * slope * (x - x0)

        ax.plot(x, y, "--", color=color, lw=1.25, alpha=0.75, zorder=30)

def contour_plot(xs, ys, field, title, label, filename, output_folder,
                 cmap="viridis", symmetric=False, add_shocks=False):
    x_mesh, y_mesh = np.meshgrid(xs, ys)

    data = np.asarray(field, dtype=float)
    data_plot = np.nan_to_num(data, nan=np.nanmedian(data[np.isfinite(data)]))

    vmin, vmax = robust_limits(data, low=0.8, high=99.2, symmetric=symmetric)

    if symmetric:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    fig, ax = plt.subplots(figsize=(8.3, 4.4))

    contour = ax.contourf(
        x_mesh,
        y_mesh,
        data_plot,
        levels=90,
        cmap=cmap,
        norm=norm,
        extend="both",
    )

    cbar = fig.colorbar(contour, ax=ax, pad=0.018, shrink=0.88)
    cbar.set_label(label)
    cbar.ax.tick_params(labelsize=8)

    draw_body(ax)

    if add_shocks:
        add_theory_shock_lines(ax, nonlinear_theory())

    clean_flow_axes(ax, xs, ys)
    ax.set_title(title)

    save_figure(fig, os.path.join(output_folder, filename))

def plot_all_contours(xs, ys, grids, folders):
    contour_plot(xs, ys, grids["Mach"], "Mach number", "Mach",
                 "mach.png", folders.field, "plasma", add_shocks=True)

    contour_plot(xs, ys, grids["p"], "Pressure", "Pressure [Pa]",
                 "pressure.png", folders.field, "viridis")

    contour_plot(xs, ys, grids["rho"], "Density", r"Density [kg m$^{-3}$]",
                 "density.png", folders.field, "viridis")

    contour_plot(xs, ys, grids["Cp"], "Pressure coefficient", r"$C_p$",
                 "cp.png", folders.field, "RdYlBu_r", symmetric=True)

    contour_plot(xs, ys, grids["T"], "Temperature", "T [K]",
                 "temperature.png", folders.field, "inferno")

    contour_plot(xs, ys, grids["V"], "Velocity magnitude", r"V [m s$^{-1}$]",
                 "velocity.png", folders.field, "viridis")

    contour_plot(xs, ys, grids["u"], "Streamwise velocity", r"u [m s$^{-1}$]",
                 "u_velocity.png", folders.field, "RdBu_r", symmetric=True)

    contour_plot(xs, ys, grids["v"], "Transverse velocity", r"v [m s$^{-1}$]",
                 "v_velocity.png", folders.field, "RdBu_r", symmetric=True)

    contour_plot(xs, ys, grids["a"], "Speed of sound", r"a [m s$^{-1}$]",
                 "sound_speed.png", folders.field, "plasma")

    contour_plot(xs, ys, grids["p0_ratio"], "Isentropic total pressure ratio", r"$p_0/p_{0,\infty}$",
                 "total_pressure_ratio.png", folders.field, "magma")

    contour_plot(xs, ys, grids["q_over_qinf"], "Dynamic pressure ratio", r"$q/q_{\infty}$",
                 "dynamic_pressure_ratio.png", folders.field, "viridis")

def plot_schlieren(xs, ys, grids, folders):
    x_mesh, y_mesh = np.meshgrid(xs, ys)

    raw = np.nan_to_num(grids["gradrho"], nan=1.0)
    raw = np.clip(raw, 1.0, np.nanpercentile(raw, 99.8))
    data = np.log10(raw)

    vmin = np.nanpercentile(data, 6)
    vmax = np.nanpercentile(data, 99.4)

    fig, ax = plt.subplots(figsize=(8.4, 4.3))

    contour = ax.contourf(
        x_mesh,
        y_mesh,
        data,
        levels=np.linspace(vmin, vmax, 90),
        cmap="hot_r",
        extend="both",
    )

    cbar = fig.colorbar(contour, ax=ax, pad=0.018, shrink=0.88)
    cbar.set_label(r"$\log_{10}|\nabla \rho|$")
    cbar.ax.tick_params(labelsize=8)

    draw_body(ax)
    add_theory_shock_lines(ax, nonlinear_theory(), color="black")

    clean_flow_axes(ax, xs, ys)
    ax.set_title("Numerical schlieren: WENO5 HLLC")

    save_figure(fig, os.path.join(folders.field, "schlieren.png"))

def plot_schlieren_shock_angles(xs, ys, grids, folders):
    """Numerical schlieren with nonlinear oblique-shock angle overlay."""
    nonlinear = nonlinear_theory()
    x_mesh, y_mesh = np.meshgrid(xs, ys)

    raw = np.nan_to_num(grids["gradrho"], nan=1.0)
    raw = np.clip(raw, 1.0, np.nanpercentile(raw, 99.85))
    data = np.log10(raw)
    vmin = np.nanpercentile(data, 8)
    vmax = np.nanpercentile(data, 99.3)

    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ax.contourf(
        x_mesh / CHORD,
        y_mesh / CHORD,
        data,
        levels=np.linspace(vmin, vmax, 90),
        cmap="Greys",
        extend="both",
    )

    draw_body(ax)

    x0 = X_LE / CHORD
    x_end = xs[-1] / CHORD
    x_line = np.linspace(x0, x_end, 300)

    upper = nonlinear["os_uf"]
    lower = nonlinear["os_lf"]

    if upper:
        upper_angle = AOA_RAD + upper["beta"]
        y_upper = (Y_CENTER / CHORD) + np.tan(upper_angle) * (x_line - x0)
        ax.plot(x_line, y_upper, ls="--", lw=1.9, color="#00A6FF",
                label=rf"LE upper shock  $\beta={upper['beta_deg']:.1f}^\circ$")

    if lower:
        lower_angle = AOA_RAD - lower["beta"]
        y_lower = (Y_CENTER / CHORD) + np.tan(lower_angle) * (x_line - x0)
        ax.plot(x_line, y_lower, ls="--", lw=1.9, color=ORANGE,
                label=rf"LE lower shock  $\beta={lower['beta_deg']:.1f}^\circ$")

    # Light guide rays from the shoulder expansion corners. They are guides, not exact fan boundaries.
    theta = nonlinear["theta_exact_rad"]
    shoulder_points = [
        (0.5, +0.5 * T_OVER_C, -theta, "#00A6FF"),
        (0.5, -0.5 * T_OVER_C, +theta, ORANGE),
    ]
    for xp, yp, ang, colour in shoulder_points:
        xx = np.linspace(xp, x_end, 220)
        yy = yp + np.tan(ang) * (xx - xp)
        ax.plot(xx, yy, ls=":", lw=1.3, color=colour, alpha=0.85)

    clean_flow_axes(ax, xs / CHORD, ys / CHORD)
    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$y/c$")
    ax.set_title("Numerical schlieren with nonlinear shock angles")
    ax.legend(frameon=True, loc="upper left", fontsize=8)
    save_figure(fig, os.path.join(folders.verification, "schlieren_shock_angles.png"))

def plot_streamlines(xs, ys, grids, folders):
    x_mesh, y_mesh = np.meshgrid(xs, ys)

    u = np.nan_to_num(grids["u"], nan=U_INF)
    v = np.nan_to_num(grids["v"], nan=VY_INF)
    mach = np.nan_to_num(grids["Mach"], nan=M_INF)

    fig, ax = plt.subplots(figsize=(8.4, 4.4))

    vmin, vmax = robust_limits(mach, 0.8, 99.2)
    contour = ax.contourf(x_mesh, y_mesh, mach, levels=90, cmap="plasma",
                          vmin=vmin, vmax=vmax)

    cbar = fig.colorbar(contour, ax=ax, pad=0.018, shrink=0.88)
    cbar.set_label("Mach")

    start_y = np.linspace(ys[0] * 0.88, ys[-1] * 0.88, 32)
    start_points = np.c_[np.full_like(start_y, xs[0]), start_y]

    ax.streamplot(
        xs,
        ys,
        u,
        v,
        start_points=start_points,
        color="white",
        linewidth=0.65,
        density=1.6,
        arrowsize=0.7,
        broken_streamlines=True,
    )

    draw_body(ax)
    clean_flow_axes(ax, xs, ys)
    ax.set_title("Mach number with streamlines")

    save_figure(fig, os.path.join(folders.field, "mach_streamlines.png"))

