# Euler reference data

This folder contains the numerical data used for the final Euler verification of the Mach 2.5 diamond airfoil case.

The solution was produced with the custom WENO5 JS + HLLC finite volume solver in `src/euler_solver/`.

## Reference case

| Quantity | Value |
|---|---:|
| Freestream Mach number | 2.5 |
| Angle of attack | 5° |
| Specific heat ratio, `γ` | 1.4 |
| Gas constant, `R` | 287 J/(kg K) |
| Freestream pressure | 101325 Pa |
| Freestream temperature | 288.15 K |
| Chord, `c` | 1 m |
| Thickness ratio, `t/c` | 0.10 |
| Grid | 720 × 360 |
| Domain in `x/c` | -1 to 3 |
| Domain in `y/c` | -1.5 to 1.5 |
| Final iteration | 40000 |

The corresponding freestream values are

```math
\rho_\infty = 1.22523\ \mathrm{kg/m^3},
```

```math
a_\infty = 340.263\ \mathrm{m/s},
```

```math
U_\infty = 850.657\ \mathrm{m/s},
```

and

```math
q_\infty
=
\frac{1}{2}\rho_\infty U_\infty^2
=
443296.875\ \mathrm{Pa}.
```

## Data files

The numerical data are stored in `data/`.

| File | Contents |
|---|---|
| `weno5_solution.csv` | Final 720 × 360 flow field |
| `weno5_forces.csv` | Final force components and aerodynamic coefficients |
| `weno5_force_history.csv` | Recorded lift and drag history |
| `weno5_residual.csv` | Residual and solution change history |

### `weno5_solution.csv`

The final flow field contains 259200 rows, corresponding to the full 720 × 360 Cartesian grid.

Columns:

```text
x
y
is_solid
rho
u
v
p
Cp
T
a
Mach
```

where:

```text
x, y       Cartesian coordinates
is_solid   0 for fluid, 1 for immersed solid
rho        density
u, v       Cartesian velocity components
p          static pressure
Cp         pressure coefficient
T          static temperature
a          local speed of sound
Mach       local Mach number
```

The pressure coefficient is

```math
C_p
=
\frac{p-p_\infty}{q_\infty}.
```

Because the airfoil is represented by a Cartesian immersed boundary, surface pressure used in the verification plots is taken from the nearest external fluid cells. It is therefore a nearest fluid cell pressure estimate rather than an exact body fitted wall value.

### `weno5_forces.csv`

Columns:

```text
Fx
Fy
Drag
Lift
Cd
Cl
```

The final aerodynamic coefficients are

```math
C_D = 0.03160637093,
```

```math
C_L = 0.15677811494.
```

For this Euler calculation the aerodynamic force comes from pressure only. There is no viscous skin friction contribution.

The Cartesian force components are resolved into freestream drag and lift using

```math
D
=
F_x\cos\alpha
+
F_y\sin\alpha,
```

```math
L
=
-F_x\sin\alpha
+
F_y\cos\alpha.
```

The force coefficients are normalised by the freestream dynamic pressure and chord.

## Analytical verification

The numerical force result is compared with both Ackeret linear theory and nonlinear shock expansion theory.

| Method | `Cd` | `Cl` |
|---|---:|---:|
| Ackeret theory | 0.030752 | 0.152345 |
| Nonlinear shock expansion theory | 0.031714 | 0.156485 |
| WENO5 JS + HLLC | 0.031606 | 0.156778 |

Relative to nonlinear shock expansion theory,

```math
\Delta C_D = -0.338\%,
```

```math
\Delta C_L = +0.187\%.
```

The leading edge shock angles provide an independent check of the resolved wave field:

| Shock | Theory | Numerical | Difference |
|---|---:|---:|---:|
| Upper | 24.091° | 24.552° | +0.461° |
| Lower | 32.531° | 32.349° | -0.182° |

The analytical panel pressure coefficients are

| Panel | `Cp` |
|---|---:|
| Upper forward | +0.011030 |
| Upper rear | -0.122886 |
| Lower forward | +0.215310 |
| Lower rear | -0.009859 |

These quantities are checked independently in the post processing scripts rather than relying on the final lift and drag alone.

## Force history

`weno5_force_history.csv` contains 801 recorded samples up to iteration 40000.

The final 12 recorded samples give

```text
Cl mean  = 0.156778114941
Cd mean  = 0.0316063709305

Cl range = 0
Cd range = 0
```

at the precision stored in the file.

The integrated aerodynamic loading is therefore stationary over the final recorded force window.

## Residual history

`weno5_residual.csv` contains one entry for each of the 40000 iterations.

At the final iteration,

```text
normalised residual = 0.0321148475769
solution change     = 0.000665051256075
```

The normalised residual has fallen by approximately 1.49 orders of magnitude, but the residual and solution change remain oscillatory at late iterations.

For that reason, the residual history is retained as a numerical diagnostic. It is not presented as evidence that a prescribed residual tolerance was reached.

The reference solution is instead assessed using the combined evidence from:

* stationary lift and drag
* analytical shock angles
* post shock pressure and Mach number
* Prandtl Meyer expansion states
* surface pressure loading
* the resolved off body wave structure

## Reproducing the figures

The scripts used to read these files and reproduce the Euler figures are stored in

```text
postprocessing/euler/
```

For example, from the repository root:

```powershell
py .\postprocessing\euler\01_surface_cp_verification.py
py .\postprocessing\euler\02_force_convergence.py
py .\postprocessing\euler\04_scalar_fields.py
py .\postprocessing\euler\08_shock_angle_verification.py
```

The full set of scripts covers the surface pressure, force history, scalar fields, numerical schlieren, axial cuts, lower shock state and shock angle verification.
