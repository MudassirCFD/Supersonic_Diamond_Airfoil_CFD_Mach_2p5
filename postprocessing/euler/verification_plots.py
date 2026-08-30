"""Plots used for analytical verification of the Euler result."""

import math
import os

import matplotlib.pyplot as plt
import numpy as np

from config import (
    BLUE, CHORD, GAS_R, GREEN, M_INF, ORANGE, P_INF, PURPLE, RED,
    RHO_INF, T_INF, save_figure,
)
from data import extract_surface_cp, lower_shock_profile_table, sample_line_at_y
from theory import ackeret, nonlinear_theory, nu_pm

def plot_riemann_fan_comparison(folders):
    """Conceptual figure: why HLLC is sharper than Rusanov for Euler waves."""
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.0), sharey=True)

    cases = [
        (axes[0], "Rusanov: two wave envelope", False),
        (axes[1], "HLLC: three wave structure", True),
    ]

    for ax, title, hllc in cases:
        ax.set_xlim(-1.5, 1.5)
        ax.set_ylim(0.0, 1.05)
        ax.axvline(0.0, color="0.55", ls="--", lw=0.9)
        ax.grid(True, color="0.88", linestyle="--", linewidth=0.7)
        ax.set_xlabel("x")
        ax.set_title(title, fontweight="bold")
        ax.set_xticks([-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5])
        ax.set_yticks([])

        t_top = 0.88
        s_l = -1.2
        s_r = 1.2
        s_star = 0.24

        ax.fill_between([s_l, 0.0], [t_top, 0.0], [0.0, 0.0], color=BLUE, alpha=0.10)
        ax.fill_between([0.0, s_r], [0.0, t_top], [0.0, 0.0], color=RED, alpha=0.08)

        ax.annotate("", xy=(s_l, t_top), xytext=(0.0, 0.0),
                    arrowprops=dict(arrowstyle="-|>", color=BLUE, lw=2.2))
        ax.annotate("", xy=(s_r, t_top), xytext=(0.0, 0.0),
                    arrowprops=dict(arrowstyle="-|>", color=RED, lw=2.2))
        ax.text(-0.62, 0.47, r"$S_L$", color=BLUE, fontsize=11, fontweight="bold")
        ax.text(0.60, 0.47, r"$S_R$", color=RED, fontsize=11, fontweight="bold")
        ax.text(0.0, 0.02, "interface", ha="center", color="0.45", fontsize=8)

        if hllc:
            ax.annotate("", xy=(s_star, t_top), xytext=(0.0, 0.0),
                        arrowprops=dict(arrowstyle="-|>", color=GREEN, lw=2.2))
            ax.text(s_star + 0.05, 0.47, r"$S^*$", color=GREEN, fontsize=11, fontweight="bold")
            ax.text(-0.55, 0.68, r"$U_L^*$", color=BLUE, fontsize=11, fontweight="bold")
            ax.text(0.50, 0.68, r"$U_R^*$", color=RED, fontsize=11, fontweight="bold")
            ax.text(0.02, 0.93,
                    "contact wave retained\nless numerical diffusion",
                    transform=ax.transAxes, ha="left", va="top", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.92, pad=3))
        else:
            ax.text(-0.55, 0.68, r"$(U_L + U_R)/2$", color="0.45", fontsize=10)
            ax.text(0.02, 0.93,
                    "single maximum signal speed\nrobust but more diffusive",
                    transform=ax.transAxes, ha="left", va="top", fontsize=8,
                    bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.92, pad=3))

    axes[0].set_ylabel("t")
    fig.suptitle("Rusanov versus HLLC: Riemann fan structure", fontsize=15, fontweight="bold")
    save_figure(fig, os.path.join(folders.theory, "rusanov_vs_hllc_riemann_fan.png"))

def plot_lower_oblique_shock_profiles(xs, ys, grids, folders, y_target=-0.35 * CHORD):
    """Horizontal cut through the lower leading-edge shock with Rankine-Hugoniot reference."""
    table = lower_shock_profile_table(xs, ys, grids, y_target=y_target)
    y_actual = float(table["y"].iloc[0])

    x = table["x_over_c"].to_numpy()
    nonlinear = nonlinear_theory()
    shock = nonlinear["os_lf"]

    rh = {
        "Mach": shock["M2"] if shock else np.nan,
        "p": P_INF * shock["p2p1"] if shock else np.nan,
        "rho": RHO_INF * shock["rho2rho1"] if shock else np.nan,
    }
    rh["T"] = rh["p"] / (rh["rho"] * GAS_R) if np.isfinite(rh["rho"]) else np.nan

    panels = [
        ("Mach_CFD", "Mach number", "Mach", M_INF, rh["Mach"]),
        ("pressure_CFD_Pa", "Pressure [Pa]", "Pressure [Pa]", P_INF, rh["p"]),
        ("density_CFD_kg_m3", r"Density [kg/m$^3$]", r"Density [kg/m$^3$]", RHO_INF, rh["rho"]),
        ("temperature_CFD_K", "Temperature [K]", "Temperature [K]", T_INF, rh["T"]),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.0), sharex=True)
    axes = axes.ravel()

    for ax, (col, title, ylabel, free_value, rh_value) in zip(axes, panels):
        y = table[col].to_numpy()
        ax.plot(x, y, color=BLUE, lw=2.0, label="CFD")
        ax.axhline(free_value, color=GREEN, ls="--", lw=1.1, label="freestream")
        ax.axhline(rh_value, color=RED, ls="--", lw=1.1, label="Rankine Hugoniot")
        ax.axvspan(0.0, 1.0, color=BLUE, alpha=0.07)
        for xpos in [0.0, 0.5, 1.0]:
            ax.axvline(xpos, color="0.78", ls="--", lw=0.6)
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.grid(True, color="0.88", linestyle="--", linewidth=0.7)

    axes[2].set_xlabel(r"$x/c$")
    axes[3].set_xlabel(r"$x/c$")
    axes[1].legend(frameon=True, loc="best", fontsize=8)
    fig.suptitle(rf"Lower oblique shock profile at $y/c={y_actual / CHORD:.3f}$", fontsize=14, fontweight="bold")
    save_figure(fig, os.path.join(folders.verification, "lower_oblique_shock_profiles.png"))

def plot_near_body_cuts(xs, ys, grids, folders):
    y_upper_target = +0.08 * CHORD
    y_lower_target = -0.08 * CHORD

    panels = [
        ("Mach", "Mach", M_INF, "Mach near-body axial cuts"),
        ("p", "Pressure [Pa]", P_INF, "Pressure near-body axial cuts"),
        ("Cp", r"$C_p$", 0.0, r"$C_p$ near-body axial cuts"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.9))

    for ax, (key, ylabel, reference, title) in zip(axes, panels):
        upper, y_upper = sample_line_at_y(xs, ys, grids[key], y_upper_target)
        lower, y_lower = sample_line_at_y(xs, ys, grids[key], y_lower_target)

        ax.plot(xs / CHORD, upper, color=BLUE, lw=2.2, label=rf"$y/c={y_upper / CHORD:+.2f}$")
        ax.plot(xs / CHORD, lower, color=ORANGE, lw=2.2, label=rf"$y/c={y_lower / CHORD:+.2f}$")
        ax.axhline(reference, color="0.35", linestyle=":", lw=1.2)

        for xpos, label in [(0.0, "LE"), (0.5, "mid"), (1.0, "TE")]:
            ax.axvline(xpos, color="0.75", lw=0.7, ls="--")
            if key == "Cp":
                ax.text(xpos + 0.01, 0.95, label, transform=ax.get_xaxis_transform(),
                        fontsize=7, color="0.35", va="top")

        ax.axvspan(0.0, 1.0, alpha=0.08, color=BLUE)
        ax.set_xlabel(r"$x/c$")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)

    axes[0].legend(frameon=True, loc="best")
    save_figure(fig, os.path.join(folders.verification, "near_body_axial_cuts.png"))

def plot_surface_cp(df, folders):
    surface = extract_surface_cp(df)

    if surface.empty:
        print("warning: no surface Cp extracted")
        return

    mask = (surface["x_over_c"] > 0.04) & (surface["x_over_c"] < 0.98)
    surface = surface.loc[mask].copy()

    ack = ackeret()
    nonlinear = nonlinear_theory()

    x = surface["x_over_c"].to_numpy()
    cp_upper = surface["Cp_upper"].to_numpy()
    cp_lower = surface["Cp_lower"].to_numpy()

    x_step = np.array([0.0, 0.5, 0.5, 1.0])

    upper_nl = np.array([nonlinear["Cp_uf"], nonlinear["Cp_uf"], nonlinear["Cp_ur"], nonlinear["Cp_ur"]])
    lower_nl = np.array([nonlinear["Cp_lf"], nonlinear["Cp_lf"], nonlinear["Cp_lr"], nonlinear["Cp_lr"]])

    upper_ack = np.array([ack["Cp_uf"], ack["Cp_uf"], ack["Cp_ur"], ack["Cp_ur"]])
    lower_ack = np.array([ack["Cp_lf"], ack["Cp_lf"], ack["Cp_lr"], ack["Cp_lr"]])

    fig, ax = plt.subplots(figsize=(8.8, 5.0))

    ax.plot(x, cp_upper, color=BLUE, lw=2.6, label="Upper CFD")
    ax.plot(x, cp_lower, color=ORANGE, lw=2.6, label="Lower CFD")

    ax.step(x_step, upper_nl, where="post", color=BLUE, ls="--", lw=1.9, label="Upper nonlinear theory")
    ax.step(x_step, lower_nl, where="post", color=ORANGE, ls="--", lw=1.9, label="Lower nonlinear theory")

    ax.step(x_step, upper_ack, where="post", color="0.35", ls=":", lw=2.0, label="Ackeret")
    ax.step(x_step, lower_ack, where="post", color="0.35", ls=":", lw=2.0)

    ax.axvline(0.5, color="0.45", lw=1.0, ls="--")

    ax.invert_yaxis()
    ax.set_xlabel(r"$x/c$")
    ax.set_ylabel(r"$C_p$")
    ax.set_title("Surface pressure coefficient")
    ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(ncol=2, frameon=True, loc="upper left", fontsize=8)

    save_figure(fig, os.path.join(folders.verification, "surface_cp.png"))

def plot_force_comparison(force_row, folders):
    ack = ackeret()
    nonlinear = nonlinear_theory()

    labels = [r"$C_d$", r"$C_l$"]
    cfd = np.array([float(force_row["Cd"]), float(force_row["Cl"])])
    ack_vals = np.array([ack["Cd"], ack["Cl"]])
    nl_vals = np.array([nonlinear["Cd"], nonlinear["Cl"]])

    x = np.arange(len(labels))
    width = 0.25

    fig, ax = plt.subplots(figsize=(7.2, 4.4))

    ax.bar(x - width, ack_vals, width=width, label="Ackeret", edgecolor="black", color="#d9d9d9")
    ax.bar(x, nl_vals, width=width, label="nonlinear OS plus PM", edgecolor="black", color=PURPLE, alpha=0.75)
    ax.bar(x + width, cfd, width=width, label="CFD", edgecolor="black", color=GREEN, alpha=0.85)

    for i, value in enumerate(cfd):
        err = 100.0 * (value - nl_vals[i]) / nl_vals[i]
        ax.text(x[i] + width, value, f"{value:.4f}\n{err:+.1f}%", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("coefficient")
    ax.set_title("Force coefficients versus theory")
    ax.grid(True, axis="y", color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True)

    save_figure(fig, os.path.join(folders.verification, "force_comparison.png"))

def plot_prandtl_meyer(folders):
    nonlinear = nonlinear_theory()

    mach_values = np.linspace(1.001, 5.0, 700)
    nu_values_deg = np.degrees([nu_pm(mach) for mach in mach_values])

    states = [
        ("Freestream", M_INF, BLUE),
        ("Upper shock", nonlinear["os_uf"]["M2"] if nonlinear["os_uf"] else M_INF, ORANGE),
        ("Lower shock", nonlinear["os_lf"]["M2"] if nonlinear["os_lf"] else M_INF, GREEN),
        ("Upper rear PM", nonlinear["pm_ur"]["M2"], RED),
        ("Lower rear PM", nonlinear["pm_lr"]["M2"], PURPLE),
    ]

    fig, ax = plt.subplots(figsize=(7.8, 4.8))

    ax.plot(mach_values, nu_values_deg, color=BLUE, linewidth=2.4, label=r"$\nu(M)$")

    for name, mach_value, color in states:
        nu_value = math.degrees(nu_pm(mach_value))
        ax.scatter(mach_value, nu_value, s=70, color=color, edgecolor="black", linewidth=0.4, zorder=4)

    table_text = "\n".join(f"{name:<14s}  M = {float(mach):.3f}" for name, mach, _ in states)

    ax.text(
        0.03,
        0.96,
        table_text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        family="monospace",
        bbox=dict(facecolor="white", edgecolor="0.70", linewidth=0.7, alpha=0.94, pad=5),
    )

    ax.set_xlabel("Mach number")
    ax.set_ylabel(r"Prandtl Meyer angle, $\nu$ [deg]")
    ax.set_title("Prandtl Meyer expansion states")
    ax.grid(True, color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=False, loc="lower right")

    save_figure(fig, os.path.join(folders.theory, "prandtl_meyer_states.png"))

def plot_theory_cp_bar(folders):
    ack = ackeret()
    nonlinear = nonlinear_theory()

    labels = ["upper front", "upper rear", "lower front", "lower rear"]
    ack_vals = [ack["Cp_uf"], ack["Cp_ur"], ack["Cp_lf"], ack["Cp_lr"]]
    nl_vals = [nonlinear["Cp_uf"], nonlinear["Cp_ur"], nonlinear["Cp_lf"], nonlinear["Cp_lr"]]

    x = np.arange(len(labels))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.bar(x - width / 2, ack_vals, width, label="Ackeret", color="#d9d9d9", edgecolor="black")
    ax.bar(x + width / 2, nl_vals, width, label="nonlinear OS plus PM", color=PURPLE, alpha=0.75, edgecolor="black")
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=15)
    ax.set_ylabel(r"$C_p$")
    ax.set_title("Theoretical surface pressure states")
    ax.grid(True, axis="y", color="0.86", linestyle="--", linewidth=0.7)
    ax.legend(frameon=True)

    save_figure(fig, os.path.join(folders.theory, "theory_cp_states.png"))

