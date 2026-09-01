import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# Paths
figure_dir = Path(__file__).resolve().parent / "figures"
figure_dir.mkdir(parents=True, exist_ok=True)


# Verified analytical panel states
panels = [
    "Upper front",
    "Upper rear",
    "Lower front",
    "Lower rear"
]

Cp_ackeret = np.array([
    0.011114742666923902,
    -0.16345956952187002,
    0.16345956952187002,
    -0.011114742666923902
])

Cp_nonlinear = np.array([
    0.011030317307464748,
    -0.1228864285912023,
    0.2153103388455268,
    -0.009858814562242172
])


# Plot
x = np.arange(len(panels))
width = 0.34

fig, ax = plt.subplots(figsize=(9.5, 6))

bars_ackeret = ax.bar(
    x - width / 2,
    Cp_ackeret,
    width,
    label="Ackeret"
)

bars_nonlinear = ax.bar(
    x + width / 2,
    Cp_nonlinear,
    width,
    label="Shock-expansion"
)


# Numerical values
ax.bar_label(
    bars_ackeret,
    fmt="%.4f",
    padding=3,
    fontsize=8
)

ax.bar_label(
    bars_nonlinear,
    fmt="%.4f",
    padding=3,
    fontsize=8
)


# Formatting
ax.axhline(
    0.0,
    color="black",
    linewidth=0.9
)

ax.set_xticks(x)
ax.set_xticklabels(panels)

ax.set_ylabel("$C_p$")

ax.set_title(
    "Analytical Panel Pressure States\n"
    r"$M_\infty = 2.5,\ \alpha = 5^\circ,\ t/c = 0.10$"
)

ax.grid(
    axis="y",
    alpha=0.25
)

ax.set_axisbelow(True)
ax.legend()


# Save
fig.tight_layout()

output_file = figure_dir / "theory_cp_states.png"

fig.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

print("\nAnalytical panel Cp states")

for panel, cp_a, cp_nl in zip(
    panels,
    Cp_ackeret,
    Cp_nonlinear
):
    difference = cp_nl - cp_a

    print(
        panel,
        "| Ackeret =", cp_a,
        "| nonlinear =", cp_nl,
        "| difference =", difference
    )

print("\nSaved:", output_file)

plt.show()
