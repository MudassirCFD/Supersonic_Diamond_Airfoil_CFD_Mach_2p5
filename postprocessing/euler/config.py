"""Reference condition and small plotting helpers for the Euler case."""

from dataclasses import dataclass
import math
import os
import re
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


GAMMA = 1.4
GAS_R = 287.0

M_INF = 2.5
P_INF = 101325.0
T_INF = 288.15
AOA_DEG = 5.0

CHORD = 1.0
T_OVER_C = 0.10
X_LE = 0.0
Y_CENTER = 0.0

AOA_RAD = math.radians(AOA_DEG)

A_INF = math.sqrt(GAMMA * GAS_R * T_INF)
V_INF = M_INF * A_INF
U_INF = V_INF * math.cos(AOA_RAD)
VY_INF = V_INF * math.sin(AOA_RAD)
RHO_INF = P_INF / (GAS_R * T_INF)
Q_INF = 0.5 * RHO_INF * V_INF * V_INF
P0_INF = P_INF * (1.0 + 0.5 * (GAMMA - 1.0) * M_INF ** 2) ** (GAMMA / (GAMMA - 1.0))
T0_INF = T_INF * (1.0 + 0.5 * (GAMMA - 1.0) * M_INF ** 2)

THETA_ACK = T_OVER_C
BETA_ACK = math.sqrt(M_INF * M_INF - 1.0)

BLUE = "tab:blue"
ORANGE = "tab:orange"
GREEN = "tab:green"
RED = "tab:red"
PURPLE = "tab:purple"
GREY = "0.45"

plt.rcParams.update({
    "figure.dpi": 130,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
})


@dataclass
class OutputFolders:
    root: str
    field: str
    verification: str
    convergence: str
    theory: str
    data: str
    animation: str


def make_output_folders(root):
    folders = OutputFolders(
        root=root,
        field=os.path.join(root, "figures_field"),
        verification=os.path.join(root, "figures_verification"),
        convergence=os.path.join(root, "figures_convergence"),
        theory=os.path.join(root, "figures_theory"),
        data=os.path.join(root, "data_extracts"),
        animation=os.path.join(root, "animations"),
    )

    for folder in folders.__dict__.values():
        os.makedirs(folder, exist_ok=True)

    return folders


def save_figure(fig, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        try:
            fig.tight_layout()
        except Exception:
            pass
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print("saved:", os.path.relpath(path))


def numeric_key(path):
    numbers = re.findall(r"\d+", os.path.basename(path))
    return int(numbers[-1]) if numbers else 0


def diamond_xy():
    thickness = T_OVER_C * CHORD
    return (
        [X_LE, X_LE + 0.5 * CHORD, X_LE + CHORD, X_LE + 0.5 * CHORD, X_LE],
        [Y_CENTER, Y_CENTER + 0.5 * thickness, Y_CENTER, Y_CENTER - 0.5 * thickness, Y_CENTER],
    )


def draw_body(ax, zorder=20):
    x, y = diamond_xy()
    ax.fill(x, y, color="white", zorder=zorder)
    ax.plot(x, y, color="black", linewidth=1.1, zorder=zorder + 1)


def robust_limits(data, low=0.8, high=99.2, symmetric=False):
    clean = np.asarray(data)
    clean = clean[np.isfinite(clean)]

    if clean.size == 0:
        return -1.0, 1.0

    vmin, vmax = np.percentile(clean, [low, high])

    if symmetric:
        limit = max(abs(vmin), abs(vmax))
        return -limit, limit

    if vmax <= vmin:
        pad = max(abs(vmin) * 0.05, 1.0e-12)
        return vmin - pad, vmax + pad

    return vmin, vmax
