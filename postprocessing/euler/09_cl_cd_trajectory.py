from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"
history_file = data_dir / "weno5_force_history.csv"

figure_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Analytical reference values
# ------------------------------------------------------------

Cl_ackeret = 0.15234482685494613
Cd_ackeret = 0.03075202535387683

Cl_nonlinear = 0.1564852721947345
Cd_nonlinear = 0.03171356478110243


# ------------------------------------------------------------
# Load force history
# ------------------------------------------------------------

history = pd.read_csv(history_file)

iteration = history["iteration"]
Cl = history["Cl"]
Cd = history["Cd"]

Cl_final = float(Cl.iloc[-1])
Cd_final = float(Cd.iloc[-1])

window = 12

Cl_recent = Cl.tail(window)
Cd_recent = Cd.tail(window)


# ------------------------------------------------------------
# Convergence summary
# ------------------------------------------------------------

print("\nAerodynamic coefficient trajectory")
print("history samples =", len(history))
print("final iteration =", int(iteration.iloc[-1]))

print("\nFinal WENO5/HLLC state")
print("Cl =", Cl_final)
print("Cd =", Cd_final)

print("\nFinal 12 recorded samples")
print("Cl mean  =", Cl_recent.mean())
print("Cd mean  =", Cd_recent.mean())
print("Cl range =", Cl_recent.max() - Cl_recent.min())
print("Cd range =", Cd_recent.max() - Cd_recent.min())


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9.5, 6))

ax.plot(
    Cd,
    Cl,
    linewidth=1.5,
    label="WENO5/HLLC trajectory"
)

ax.scatter(
    Cd.iloc[0],
    Cl.iloc[0],
    s=45,
    marker="o",
    label="Initial recorded state",
    zorder=10
)

ax.scatter(
    Cd_final,
    Cl_final,
    s=80,
    marker="*",
    label="Final WENO5/HLLC",
    zorder=12
)

ax.scatter(
    Cd_ackeret,
    Cl_ackeret,
    s=55,
    marker="s",
    label="Ackeret theory",
    zorder=11
)

ax.scatter(
    Cd_nonlinear,
    Cl_nonlinear,
    s=55,
    marker="D",
    label="Shock-expansion theory",
    zorder=11
)


# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

ax.set_xlabel("$C_D$")
ax.set_ylabel("$C_L$")

ax.set_title(
    "Aerodynamic Coefficient Convergence Trajectory\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)

ax.grid(True, alpha=0.25)
ax.set_axisbelow(True)

ax.legend(
    fontsize=8,
    loc="best"
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fig.tight_layout()

output_file = figure_dir / "cl_vs_cd_overview.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:", output_file)

plt.show()
