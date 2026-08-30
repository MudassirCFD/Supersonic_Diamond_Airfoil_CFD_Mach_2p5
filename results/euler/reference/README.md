# WENO5-HLLC Euler reference results

Reference case for the two-dimensional Mach 2.5 diamond-airfoil Euler calculation.

## Case

| Quantity | Value |
|---|---:|
| Mach number | 2.5 |
| Angle of attack | 5 deg |
| Gamma | 1.4 |
| Gas constant | 287 J/(kg K) |
| Freestream pressure | 101325 Pa |
| Freestream temperature | 288.15 K |
| Chord | 1 m |
| Thickness ratio | 0.10 |
| Cartesian grid | 720 x 360 |
| Domain | x = [-1, 3], y = [-1.5, 1.5] |

## Reference coefficients

| Method | Cd | Cl |
|---|---:|---:|
| WENO5-JS + HLLC | 0.031606 | 0.156778 |
| Ackeret theory | 0.030752 | 0.152345 |
| Nonlinear shock-expansion theory | 0.031714 | 0.156485 |

Relative difference of the Euler solution:

| Reference | Cd | Cl |
|---|---:|---:|
| Ackeret | +2.78% | +2.91% |
| Nonlinear shock-expansion | -0.34% | +0.19% |

The nonlinear reference uses oblique-shock relations on the forward panels and Prandtl-Meyer expansion relations across the two shoulders.

## Archived result files

The compact reference dataset contains:

- `weno5_forces.csv`
- `weno5_force_history.csv`
- `weno5_residual.csv`
- `surface_cp.csv`
- `panel_forces.csv`
- `weno5_diagnostics.csv`

Large full-field CSV files and snapshot sequences are not stored in the normal Git repository. They can be regenerated using the public solver.

## Post-processing

The corresponding post-processing scripts are available under:

`postprocessing/euler/`

They generate field contours, numerical schlieren, shock-angle checks, surface pressure comparisons, force and residual histories, analytical verification plots and optional animations.
