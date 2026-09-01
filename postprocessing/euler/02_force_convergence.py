from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt


# Paths
repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"
figure_dir = Path(__file__).resolve().parent / "figures"

force_history_file = data_dir / "weno5_force_history.csv"
forces_file = data_dir / "weno5_forces.csv"

figure_dir.mkdir(parents=True, exist_ok=True)


# Load solver data
history = pd.read_csv(force_history_file)
forces = pd.read_csv(forces_file)

iteration = history["iteration"]
Cl = history["Cl"]
Cd = history["Cd"]

Cl_final = float(forces["Cl"].iloc[-1])
Cd_final = float(forces["Cd"].iloc[-1])


# Late-iteration convergence check
window = 12

Cl_recent = Cl.tail(window)
Cd_recent = Cd.tail(window)

Cl_mean = Cl_recent.mean()
Cd_mean = Cd_recent.mean()

Cl_std = Cl_recent.std()
Cd_std = Cd_recent.std()

Cl_drift = 100.0 * (Cl_recent.iloc[-1] - Cl_recent.iloc[0]) / Cl_mean
Cd_drift = 100.0 * (Cd_recent.iloc[-1] - Cd_recent.iloc[0]) / Cd_mean


print("\nForce convergence")
print("history samples =", len(history))
print("final iteration =", int(iteration.iloc[-1]))

print("\nLift coefficient")
print("final =", Cl_final)
print("last 12 mean =", Cl_mean)
print("last 12 std  =", Cl_std)
print("last 12 drift =", Cl_drift, "%")

print("\nDrag coefficient")
print("final =", Cd_final)
print("last 12 mean =", Cd_mean)
print("last 12 std  =", Cd_std)
print("last 12 drift =", Cd_drift, "%")


# Force convergence plot
fig, ax = plt.subplots(figsize=(9.5, 6.0))

ax.plot(
    iteration,
    Cl,
    linewidth=1.5,
    label="$C_L$"
)

ax.plot(
    iteration,
    Cd,
    linewidth=1.5,
    label="$C_D$"
)

ax.axhline(
    Cl_final,
    linestyle="--",
    linewidth=1.0,
    label=f"Final $C_L$ = {Cl_final:.6f}"
)

ax.axhline(
    Cd_final,
    linestyle="--",
    linewidth=1.0,
    label=f"Final $C_D$ = {Cd_final:.6f}"
)

ax.set_xlabel("Iteration")
ax.set_ylabel("Aerodynamic coefficient")

ax.set_title(
    "Aerodynamic Force Convergence\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)

ax.grid(True, alpha=0.3)
ax.set_axisbelow(True)
ax.legend()

fig.tight_layout()

output_file = figure_dir / "force_history.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved figure")
print(output_file)

plt.show()
