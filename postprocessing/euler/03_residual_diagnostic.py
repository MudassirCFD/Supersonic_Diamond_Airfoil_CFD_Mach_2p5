from pathlib import Path

import math
import pandas as pd
import matplotlib.pyplot as plt


# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------

repo_root = Path(__file__).resolve().parents[2]
data_dir = repo_root / "results" / "euler" / "reference" / "data"

figure_dir = Path(__file__).resolve().parent / "figures"
diagnostic_dir = figure_dir / "diagnostics"

residual_file = data_dir / "weno5_residual.csv"

diagnostic_dir.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Load residual history
# ------------------------------------------------------------

history = pd.read_csv(residual_file)

iteration = history["iteration"]
normalized_residual = history["normalized_residual"]
solution_change = history["solution_change"]


# ------------------------------------------------------------
# Diagnostic statistics
# ------------------------------------------------------------

initial_residual = float(normalized_residual.iloc[0])
final_residual = float(normalized_residual.iloc[-1])

initial_change = float(solution_change.iloc[0])
final_change = float(solution_change.iloc[-1])

if initial_residual > 0.0 and final_residual > 0.0:
    residual_reduction = math.log10(
        initial_residual / final_residual
    )
else:
    residual_reduction = float("nan")


# ------------------------------------------------------------
# Console summary
# ------------------------------------------------------------

print("\nResidual diagnostic")
print("history samples =", len(history))
print("final iteration =", int(iteration.iloc[-1]))

print("\nNormalised residual")
print("initial =", initial_residual)
print("final   =", final_residual)
print("reduction =", residual_reduction, "orders")

print("\nSolution change")
print("initial =", initial_change)
print("final   =", final_change)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

fig, ax = plt.subplots(figsize=(9.5, 6))

ax.semilogy(
    iteration,
    normalized_residual,
    linewidth=1.4,
    label="Normalised residual"
)

ax.semilogy(
    iteration,
    solution_change,
    linewidth=1.4,
    label="Solution change"
)

ax.axhline(
    5.0e-6,
    color="0.35",
    linestyle="--",
    linewidth=1.0,
    label=r"Solution-change target $5\times10^{-6}$"
)


# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------

ax.set_xlabel("Iteration")
ax.set_ylabel("Convergence measure")

ax.set_title(
    "Euler Solver Residual Diagnostic\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ$"
)

ax.grid(
    True,
    which="both",
    alpha=0.25
)

ax.set_axisbelow(True)

ax.legend(
    fontsize=8,
    loc="best"
)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

fig.tight_layout()

output_file = diagnostic_dir / "residual_history.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nSaved:", output_file)

plt.show()
