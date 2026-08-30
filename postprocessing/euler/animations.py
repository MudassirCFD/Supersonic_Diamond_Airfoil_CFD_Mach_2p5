"""Optional GIF generation from saved solver snapshots."""

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import Normalize, TwoSlopeNorm

from config import (
    BLUE, GREEN, PURPLE, RED, draw_body, numeric_key,
)
from data import (
    find_snapshot_files, load_snapshot_grid, load_snapshot_schlieren,
)
from theory import ackeret, nonlinear_theory

def scan_animation_range(files, field_name, symmetric=False):
    vmin = np.inf
    vmax = -np.inf

    for file in files:
        if field_name == "schlieren":
            _, _, grid, _ = load_snapshot_schlieren(file)
        else:
            _, _, grid = load_snapshot_grid(file, field_name)

        vmin = min(vmin, np.nanpercentile(grid, 1.0))
        vmax = max(vmax, np.nanpercentile(grid, 99.0))

    if symmetric:
        limit = max(abs(vmin), abs(vmax))
        vmin = -limit
        vmax = limit

    return vmin, vmax

def make_field_animation(files, folders, field_name, title, cmap, fps=8, symmetric=False):
    if not files:
        return

    print("creating animation:", field_name)

    xs, ys, first = load_snapshot_grid(files[0], field_name)
    extent = [xs[0], xs[-1], ys[0], ys[-1]]

    vmin, vmax = scan_animation_range(files, field_name, symmetric=symmetric)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))

    first_data = np.nan_to_num(first, nan=np.nanmedian(first[np.isfinite(first)]))

    if symmetric:
        norm = TwoSlopeNorm(vmin=vmin, vcenter=0.0, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)

    image = ax.imshow(
        first_data,
        extent=extent,
        origin="lower",
        cmap=cmap,
        norm=norm,
        interpolation="bilinear",
        aspect="equal",
    )

    cbar = fig.colorbar(image, ax=ax, shrink=0.86, pad=0.018)
    cbar.set_label(field_name)

    draw_body(ax)

    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ys[0], ys[-1])
    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$y/c$")

    step_text = ax.text(0.02, 0.94, "", transform=ax.transAxes, fontsize=9,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.70, pad=2))

    def update(frame_index):
        _, _, grid = load_snapshot_grid(files[frame_index], field_name)
        data = np.nan_to_num(grid, nan=np.nanmedian(grid[np.isfinite(grid)]))
        image.set_data(data)

        step = numeric_key(files[frame_index])
        ax.set_title(f"{title}: step {step}")
        step_text.set_text(f"step {step}")

        return [image, step_text]

    animation = FuncAnimation(fig, update, frames=len(files), interval=1000 / fps, blit=False)

    path = os.path.join(folders.animation, f"{field_name.lower()}_animation.gif")
    animation.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print("saved:", os.path.relpath(path))

def make_schlieren_animation(files, folders, fps=8):
    if not files:
        return

    print("creating animation: schlieren")

    xs, ys, first, label = load_snapshot_schlieren(files[0])
    extent = [xs[0], xs[-1], ys[0], ys[-1]]

    vmin, vmax = scan_animation_range(files, "schlieren", symmetric=False)

    fig, ax = plt.subplots(figsize=(8.2, 4.6))

    image = ax.imshow(
        first,
        extent=extent,
        origin="lower",
        cmap="hot_r",
        vmin=vmin,
        vmax=vmax,
        interpolation="bilinear",
        aspect="equal",
    )

    cbar = fig.colorbar(image, ax=ax, shrink=0.86, pad=0.018)
    cbar.set_label(r"$\log_{10}|\nabla " + label + "|$")

    draw_body(ax)

    ax.set_xlim(xs[0], xs[-1])
    ax.set_ylim(ys[0], ys[-1])
    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$y/c$")

    step_text = ax.text(0.02, 0.94, "", transform=ax.transAxes, fontsize=9,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.70, pad=2))

    def update(frame_index):
        _, _, grid, _ = load_snapshot_schlieren(files[frame_index])
        image.set_data(grid)

        step = numeric_key(files[frame_index])
        ax.set_title(f"Numerical schlieren: WENO5 HLLC")
        step_text.set_text(f"step {step}")

        return [image, step_text]

    animation = FuncAnimation(fig, update, frames=len(files), interval=1000 / fps, blit=False)

    path = os.path.join(folders.animation, "schlieren_animation.gif")
    animation.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print("saved:", os.path.relpath(path))

def make_cl_cd_animation(force_history_df, force_row, folders, fps=10):
    if len(force_history_df) < 2:
        return

    print("creating animation: cl vs cd")

    ack = ackeret()
    nonlinear = nonlinear_theory()

    cd = force_history_df["Cd"].to_numpy(dtype=float)
    cl = force_history_df["Cl"].to_numpy(dtype=float)

    final_cd = float(force_row["Cd"])
    final_cl = float(force_row["Cl"])

    cd_all = np.array([np.nanmin(cd), np.nanmax(cd), ack["Cd"], nonlinear["Cd"], final_cd])
    cl_all = np.array([np.nanmin(cl), np.nanmax(cl), ack["Cl"], nonlinear["Cl"], final_cl])

    cd_pad = 0.10 * max(np.ptp(cd_all), 1.0e-4)
    cl_pad = 0.10 * max(np.ptp(cl_all), 1.0e-4)

    fig, ax = plt.subplots(figsize=(6.4, 5.6))
    ax.set_xlim(np.nanmin(cd_all) - cd_pad, np.nanmax(cd_all) + cd_pad)
    ax.set_ylim(np.nanmin(cl_all) - cl_pad, np.nanmax(cl_all) + cl_pad)
    ax.set_xlabel(r"$C_d$")
    ax.set_ylabel(r"$C_l$")
    ax.set_title(r"Animated $C_l$ versus $C_d$ convergence")

    ax.scatter(ack["Cd"], ack["Cl"], s=90, marker="x", color=RED, lw=2.0, label="Ackeret")
    ax.scatter(nonlinear["Cd"], nonlinear["Cl"], s=80, marker="D", color=PURPLE, edgecolor="black", lw=0.5, label="nonlinear theory")
    ax.scatter(final_cd, final_cl, s=130, marker="*", color=GREEN, edgecolor="black", lw=0.5, label="final CFD")

    line, = ax.plot([], [], color=BLUE, lw=2.0, label="CFD trajectory")
    point, = ax.plot([], [], "o", color=BLUE, markersize=6)

    ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True, loc="best")

    def update(frame_index):
        line.set_data(cd[:frame_index + 1], cl[:frame_index + 1])
        point.set_data([cd[frame_index]], [cl[frame_index]])

        iteration = int(force_history_df["iteration"].iloc[frame_index])
        ax.set_title(rf"$C_l$ versus $C_d$ convergence | iteration {iteration}")

        return [line, point]

    animation = FuncAnimation(fig, update, frames=len(force_history_df), interval=1000 / fps, blit=False)

    path = os.path.join(folders.animation, "cl_vs_cd_animation.gif")
    animation.save(path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print("saved:", os.path.relpath(path))

def make_animations(input_folder, folders, force_history_df, force_row, fps=8, stride=1):
    files = find_snapshot_files(input_folder, stride=stride)
    if not files:
        print("No snapshot files found. Animations skipped.")
        return

    for field_name, title, cmap, symmetric in [
        ("Mach", "Mach number evolution", "plasma", False),
        ("Cp", r"$C_p$ evolution", "RdYlBu_r", True),
        ("p", "Pressure evolution", "viridis", False),
        ("rho", "Density evolution", "viridis", False),
    ]:
        try:
            make_field_animation(files, folders, field_name, title, cmap, fps=fps, symmetric=symmetric)
        except Exception as exc:
            print(f"animation skipped for {field_name}: {exc}")

    try:
        make_schlieren_animation(files, folders, fps=fps)
    except Exception as exc:
        print(f"schlieren animation skipped: {exc}")

    try:
        make_cl_cd_animation(force_history_df, force_row, folders, fps=10)
    except Exception as exc:
        print(f"Cl-Cd animation skipped: {exc}")

