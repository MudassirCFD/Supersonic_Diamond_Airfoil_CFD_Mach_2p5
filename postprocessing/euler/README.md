# Euler post processing

These scripts reproduce the figures used to examine and verify the Mach 2.5 Euler reference solution.

The numerical data are stored in

```text
results/euler/reference/data/
```

The scripts should be run from the repository root.

## Python requirements

```text
numpy
pandas
matplotlib
```

## Scripts

| Script | Purpose |
|---|---|
| `01_surface_cp_verification.py` | Compares numerical surface pressure and integrated forces with Ackeret and nonlinear shock expansion theory |
| `02_force_convergence.py` | Checks the lift and drag history and the final force window |
| `03_residual_diagnostic.py` | Examines the residual and solution change history |
| `04_scalar_fields.py` | Produces Mach, pressure, density, temperature, velocity, speed of sound and pressure coefficient fields |
| `05_numerical_schlieren.py` | Produces numerical schlieren from the density gradient |
| `06_near_body_axial_cuts.py` | Samples the numerical field close to the upper and lower surfaces |
| `07_lower_shock_verification.py` | Compares the lower shock crossing and post shock state with oblique shock theory |
| `08_shock_angle_verification.py` | Extracts the leading edge shock angles from the numerical density gradient field |
| `09_cl_cd_trajectory.py` | Shows the path of the aerodynamic coefficients towards the final state |
| `10_theory_cp_states.py` | Compares the four analytical panel pressure states from Ackeret and nonlinear theory |

## Running the scripts

From the repository root:

```powershell
py .\postprocessing\euler\01_surface_cp_verification.py
```

or, for example,

```powershell
py .\postprocessing\euler\08_shock_angle_verification.py
```

The generated figures are written to

```text
postprocessing/euler/figures/
```

The residual history is kept separately in

```text
postprocessing/euler/figures/diagnostics/
```

because it is used as a numerical diagnostic rather than the main evidence for convergence.

## Main verification figures

The strongest checks of the Euler reference solution are:

```text
schlieren_shock_angles.png
lower_oblique_shock_profiles.png
surface_cp.png
force_comparison.png
force_history.png
```

They examine different parts of the same physical chain:

```text
flow turning
    ↓
shock and expansion states
    ↓
surface pressure
    ↓
lift and wave drag
```

The shock angle is measured from a fitted density gradient ridge in `08_shock_angle_verification.py`.

The horizontal lower shock cut in `07_lower_shock_verification.py` has a different purpose. It checks the local shock crossing and the change in Mach number, pressure and pressure coefficient across the shock. A single horizontal cut is not used to define the global shock angle.

## Surface pressure

The airfoil is represented using a Cartesian immersed boundary.

For this reason, the numerical surface pressure used in `01_surface_cp_verification.py` is taken from the nearest external fluid cell above or below the immersed surface.

It should therefore be interpreted as a nearest fluid cell pressure estimate rather than an exact body fitted wall value.

## Convergence

The final lift and drag are stationary over the final 12 recorded force samples.

The residual and solution change histories remain oscillatory at late iterations. They are retained and plotted honestly in `03_residual_diagnostic.py`, but they are not used alone to claim convergence.

The Euler reference solution is judged from the combined evidence of stationary forces, correct shock geometry, correct post shock states, expansion states and surface pressure loading.
