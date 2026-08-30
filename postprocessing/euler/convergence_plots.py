"""Residual and aerodynamic-coefficient convergence plots."""

import os

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from config import BLUE, GREEN, GREY, ORANGE, PURPLE, RED, save_figure
from theory import ackeret, nonlinear_theory

def moving_average_reflect(y, window):
    y = np.asarray(y, dtype=float)
    if len(y) < 3:
        return y

    window = int(max(3, window))
    if window % 2 == 0:
        window += 1

    if len(y) < window:
        window = len(y) if len(y) % 2 == 1 else len(y) - 1

    if window < 3:
        return y

    pad = window // 2
    y_pad = np.pad(y, pad_width=pad, mode="edge")
    kernel = np.ones(window) / window
    return np.convolve(y_pad, kernel, mode="valid")

def plot_residual(residual_df, folders):
    fig, ax = plt.subplots(figsize=(9.5, 4.6))

    iteration = residual_df["iteration"].to_numpy()
    window = max(51, len(iteration) // 250)

    curves = [
        ("residual", "absolute residual", BLUE),
        ("normalized_residual", "normalised residual", ORANGE),
        ("solution_change", "solution change", GREEN),
    ]

    for column, label, colour in curves:
        if column not in residual_df.columns:
            continue

        y = residual_df[column].to_numpy(dtype=float)
        smooth = moving_average_reflect(y, window)

        stride = max(1, len(y) // 2500)
        ax.semilogy(iteration[::stride], y[::stride], color=colour, alpha=0.18, lw=0.7)
        ax.semilogy(iteration, smooth, color=colour, lw=2.3, label=label)

    ax.set_xlabel("iteration")
    ax.set_ylabel("residual")
    ax.set_title("Convergence history")
    ax.grid(True, which="both", color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True, loc="best")

    save_figure(fig, os.path.join(folders.convergence, "residual_history.png"))

def plot_force_history(force_history_df, force_row, folders):
    ack = ackeret()
    nonlinear = nonlinear_theory()

    fig, axes = plt.subplots(2, 1, figsize=(8.8, 6.4), sharex=True)

    axes[0].plot(force_history_df["iteration"], force_history_df["Cd"], color=BLUE, lw=2.1, label="CFD")
    axes[0].axhline(ack["Cd"], linestyle="--", color="0.25", label="Ackeret")
    axes[0].axhline(nonlinear["Cd"], linestyle=":", color="0.25", lw=2.0, label="nonlinear OS plus PM")
    axes[0].axhline(float(force_row["Cd"]), linestyle="-.", color=GREEN, lw=1.7, label="final force file")
    axes[0].set_ylabel(r"$C_d$")
    axes[0].set_title("Drag coefficient convergence")
    axes[0].grid(True, color="0.86", linestyle="--", linewidth=0.7)
    axes[0].legend(frameon=True, loc="best", ncol=2)

    axes[1].plot(force_history_df["iteration"], force_history_df["Cl"], color=BLUE, lw=2.1, label="CFD")
    axes[1].axhline(ack["Cl"], linestyle="--", color="0.25", label="Ackeret")
    axes[1].axhline(nonlinear["Cl"], linestyle=":", color="0.25", lw=2.0, label="nonlinear OS plus PM")
    axes[1].axhline(float(force_row["Cl"]), linestyle="-.", color=GREEN, lw=1.7, label="final force file")
    axes[1].set_ylabel(r"$C_l$")
    axes[1].set_xlabel("iteration")
    axes[1].set_title("Lift coefficient convergence")
    axes[1].grid(True, color="0.86", linestyle="--", linewidth=0.7)

    save_figure(fig, os.path.join(folders.convergence, "force_history.png"))

def plot_cl_vs_cd(force_history_df, force_row, folders):
    ack = ackeret()
    nonlinear = nonlinear_theory()

    cd = force_history_df["Cd"].to_numpy(dtype=float)
    cl = force_history_df["Cl"].to_numpy(dtype=float)

    final_cd = float(force_row["Cd"])
    final_cl = float(force_row["Cl"])

    # Overview plot
    fig, ax = plt.subplots(figsize=(6.7, 5.8))
    ax.plot(cd, cl, color=BLUE, lw=2.0, label="CFD trajectory")
    ax.scatter(cd[0], cl[0], s=50, color=GREY, label="start", zorder=4)
    ax.scatter(final_cd, final_cl, s=160, marker="*", color=GREEN, edgecolor="black", lw=0.6, label="final CFD", zorder=5)
    ax.scatter(ack["Cd"], ack["Cl"], s=95, marker="x", color=RED, lw=2.2, label="Ackeret", zorder=5)
    ax.scatter(nonlinear["Cd"], nonlinear["Cl"], s=90, marker="D", color=PURPLE, edgecolor="black", lw=0.5, label="nonlinear theory", zorder=5)
    ax.set_xlabel(r"$C_d$")
    ax.set_ylabel(r"$C_l$")
    ax.set_title(r"Full $C_l$ versus $C_d$ convergence path")
    ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True, loc="upper left")
    save_figure(fig, os.path.join(folders.convergence, "cl_vs_cd_overview.png"))

    # Converged-region plot
    tail = max(20, len(cd) // 5)
    cd_tail = cd[-tail:]
    cl_tail = cl[-tail:]

    cd_targets = np.array([final_cd, ack["Cd"], nonlinear["Cd"], *cd_tail])
    cl_targets = np.array([final_cl, ack["Cl"], nonlinear["Cl"], *cl_tail])

    cd_pad = max(0.0025, 0.16 * np.ptp(cd_targets))
    cl_pad = max(0.0040, 0.16 * np.ptp(cl_targets))

    fig, ax = plt.subplots(figsize=(6.8, 5.6))
    ax.plot(cd_tail, cl_tail, color=BLUE, lw=2.0, label="final CFD trajectory")
    ax.scatter(final_cd, final_cl, s=160, marker="*", color=GREEN, edgecolor="black", lw=0.6, label="final CFD", zorder=5)
    ax.scatter(ack["Cd"], ack["Cl"], s=95, marker="x", color=RED, lw=2.2, label="Ackeret", zorder=5)
    ax.scatter(nonlinear["Cd"], nonlinear["Cl"], s=90, marker="D", color=PURPLE, edgecolor="black", lw=0.5, label="nonlinear theory", zorder=5)
    ax.set_xlim(np.nanmin(cd_targets) - cd_pad, np.nanmax(cd_targets) + cd_pad)
    ax.set_ylim(np.nanmin(cl_targets) - cl_pad, np.nanmax(cl_targets) + cl_pad)
    ax.set_xlabel(r"$C_d$")
    ax.set_ylabel(r"$C_l$")
    ax.set_title(r"Converged $C_l$ versus $C_d$ comparison")
    ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True, loc="best")

    inset = inset_axes(ax, width="37%", height="37%", loc="lower right")
    inset.plot(cd, cl, color=BLUE, lw=1.3)
    inset.scatter(cd[0], cl[0], s=20, color=GREY)
    inset.scatter(final_cd, final_cl, s=60, marker="*", color=GREEN, edgecolor="black", lw=0.4)
    inset.set_title("full path", fontsize=7)
    inset.tick_params(labelsize=6)
    inset.grid(True, color="0.88", linestyle="--", linewidth=0.5)

    save_figure(fig, os.path.join(folders.convergence, "cl_vs_cd_converged.png"))

