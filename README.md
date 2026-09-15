# Supersonic Diamond Airfoil CFD at Mach 2.5

**From compressible-flow theory to a custom C++ shock-capturing solver and evidence-driven wall-resolved SA-RANS development**

This project studies a two-dimensional diamond airfoil at Mach 2.5 and 5° angle of attack.

The geometry creates a clear system of compression waves, oblique shocks and Prandtl-Meyer expansions. This makes it a useful case for checking whether a numerical solution recovers the physics predicted by classical compressible-flow theory [1].

A custom C++ finite-volume Euler solver was developed using HLLC intercell fluxes [2,3], WENO5-JS reconstruction [4] and SSP-RK3 time integration [5]. The airfoil is represented on a Cartesian grid using an immersed-boundary treatment [8].

The verified Euler solution provides the inviscid reference, while the wall-resolved SA-RANS branch is being rebuilt around a redesigned hybrid mesh after the original structured RANS mesh failed the adopted convergence and mesh-quality criteria.

## Numerical development of the wave field

<p align="center">
  <img src="./animations/weno5_hllc_schlieren_smoke.gif"
       alt="Numerical schlieren evolution of the WENO5-JS HLLC solution"
       width="900">
</p>

<p align="center">
  <em>Evolution of the numerical schlieren field from the initial transient towards the established Mach 2.5 shock-expansion structure.</em>
</p>

The solution should not be accepted only because lift and drag become stable.

The shock angles, Rankine-Hugoniot states, Prandtl-Meyer expansion states, surface pressure, off-body wave structure and integrated aerodynamic forces were checked against independent analytical references [1].

> **Main engineering question:**  
> Does the same aerodynamic conclusion survive when the modelling fidelity is increased?

---

## 1. Why this problem?

A diamond airfoil is simple geometrically, but the flow physics are not.

At supersonic speed, each change in surface direction creates either a compression or an expansion.

For this case:

- the leading edges generate attached oblique shocks;
- the mid-chord corners generate expansion waves;
- the upper and lower surfaces experience different pressure states because the airfoil is at 5° incidence;
- those pressure differences create lift, wave drag and pitching moment.

The important advantage is that these states can be calculated independently using oblique-shock and Prandtl-Meyer theory [1].

This provides an independent reference outside the CFD solver.

The solver therefore has to do more than produce a smooth-looking flow field. It has to recover the correct wave system and produce aerodynamic forces that are consistent with that physics.

---

## 2. Reference case

| Quantity | Value |
|---|---:|
| Geometry | 2D diamond airfoil |
| Chord, `c` | 1.0 m |
| Thickness ratio, `t/c` | 0.10 |
| Freestream Mach number | 2.5 |
| Angle of attack | 5° |
| Ratio of specific heats, `γ` | 1.4 |
| Gas constant, `R` | 287 J kg⁻¹ K⁻¹ |
| Freestream pressure | 101325 Pa |
| Freestream temperature | 288.15 K |
| Freestream velocity | 850.657 m/s |
| Euler reference grid | 720 × 360 |
| RANS chord Reynolds number | approximately 5.82 × 10⁷ |

---

## 3. Solver development

### 3.1 Euler formulation

The inviscid branch solves the two-dimensional compressible Euler equations in conservative form,

```math
\frac{\partial \mathbf{U}}{\partial t}
+
\frac{\partial \mathbf{F}}{\partial x}
+
\frac{\partial \mathbf{G}}{\partial y}
=0,
```
### 3.2 Why Rusanov was replaced by HLLC
Rusanov was useful as a robust starting point.

It uses a single maximum signal speed and adds relatively strong numerical dissipation. This makes it simple and stable, but the same dissipation can smear shocks and contact structures [2].

The earlier compressible-flow studies showed this behaviour directly. Rusanov remained robust, while HLLC recovered discontinuous wave structure more accurately.

For this airfoil, that difference matters.

The pressure change across each wave contributes directly to the aerodynamic force. If the numerical method smears the shock too strongly, it can also change the surface pressure and therefore the predicted wave drag and lift.

HLLC was therefore adopted.

---

### 3.3 HLLC flux

HLLC restores the intermediate contact wave that is missing from the simpler two-wave HLL representation [2,3].

The approximate wave system contains

```math
S_L,\qquad S_*,\qquad S_R.
```
### 3.4 Barth-Jespersen reconstruction

Before moving to WENO5-JS, an HLLC solver using Barth-Jespersen limited reconstruction [6] was developed.

The limiter allows higher-order reconstruction in smooth regions, but reduces the reconstruction close to strong gradients where non-physical oscillations can appear.

This provided an important intermediate solver:

**HLLC + Barth-Jespersen**

This branch was used to check the HLLC flux, pressure field, force integration and convergence behaviour before adding the more expensive WENO5-JS reconstruction.

### 3.5 Why WENO5-JS?

The flow contains two very different numerical regions.

High-order accuracy is required in smooth regions.

Across shocks, the solution is discontinuous and a normal high-order reconstruction can create non-physical oscillations.

WENO5-JS handles this by reconstructing the solution from several candidate stencils. In smooth regions, the stencils combine to recover fifth-order accuracy. Close to a discontinuity, the nonlinear weights reduce the influence of stencils that cross the shock [4].

For this problem, that matters because the numerical method must preserve:

- the leading-edge compression waves;
- the sharp pressure rise across the oblique shocks;
- the mid-chord expansion structure;
- the pressure distribution that produces wave drag and lift.

The WENO5-JS formulation was therefore adopted when the improvement in local wave resolution justified the additional computational cost.

The final Euler solver uses:

```text
Finite-volume conservation
        ↓
WENO5-JS reconstruction
        ↓
HLLC intercell flux
        ↓
SSP-RK3 time integration
```

### 3.6 Build and run

The Euler solver is written in C++17 and uses only the standard library.

Compile from the repository root with:

```bash
g++ -O3 -std=c++17 src/euler_solver/main.cpp -o diamond_solver
```

Run the reference case with:

```bash
./diamond_solver \
    --max-iters 40000 \
    --min-iters 40000 \
    --output results/euler/reference
```

The default solver configuration already contains the reference condition described in Section 2: Mach 2.5, 5° angle of attack, `t/c = 0.10` and the `720 × 360` Cartesian grid.

Available runtime options can be viewed with:

```bash
./diamond_solver --help
```

---

## 4. WENO5-JS + HLLC flow field

The final Euler solution recovers the expected asymmetric supersonic wave system.

At the leading edge, both surfaces turn the flow towards the body and create attached oblique shocks.

Because the airfoil is at 5° angle of attack, the effective turning angle is different on the two sides:

- the lower forward surface produces the stronger compression;
- the upper forward surface produces the weaker compression.

The pressure rise is therefore much larger on the lower forward panel.

At the mid-chord corners, the surfaces turn away from the local flow. The compressed flow then expands through Prandtl-Meyer fans before leaving the trailing edge.

This creates the main aerodynamic loading pattern:

**strong lower-surface compression + upper-surface pressure reduction → positive lift**

while the streamwise component of the pressure loading produces wave drag.

The assessment does not begin with the final `C_D` and `C_L`.

The physical wave system is checked first.

> **A correct force coefficient with the wrong shock or expansion structure would not be an acceptable solution.**

The next checks therefore compare the numerical wave field directly with analytical compressible-flow theory.

---

## 5. Verification against compressible-flow theory

Several independent checks were used because no single CFD quantity is enough to prove that the full solution is physically correct.

The first check is the leading-edge shock geometry.

### 5.1 Leading-edge shock verification

For the symmetric diamond geometry,

```text
thickness ratio, t/c = 0.10
half angle, δ         = atan(0.10)
                      ≈ 5.71 deg
```

At `α = 5°`, the upper and lower forward panels experience different effective compression turns:

```text
upper surface turning angle ≈ 0.71 deg
lower surface turning angle ≈ 10.71 deg
```

The lower surface therefore produces the substantially stronger leading-edge compression.

For an attached oblique shock, the shock angle `β` is obtained from the nonlinear `θ-β-M` relation,

```math
\tan\theta
=
2\cot\beta
\left[
\frac{M_\infty^2\sin^2\beta-1}
{M_\infty^2(\gamma+\cos 2\beta)+2}
\right].
```

Using the weak-shock solution gives:

| Surface | Turning angle, `θ` | Theory shock angle, `β` |
|---|---:|---:|
| Upper forward panel | 0.71° | 24.09° |
| Lower forward panel | 10.71° | 32.53° |

The corresponding Rankine-Hugoniot states are:

| Quantity | Upper shock | Lower shock |
|---|---:|---:|
| `θ` | 0.71° | 10.71° |
| `β` | 24.09° | 32.53° |
| `Mn1` | 1.020 | 1.344 |
| `p2/p1` | 1.048 | 1.942 |
| `ρ2/ρ1` | 1.034 | 1.593 |
| `T2/T1` | 1.014 | 1.219 |
| `M2` | 2.470 | 2.056 |
| `p2` | 106.2 kPa | 196.8 kPa |

The upper shock is weak because the local flow turning is only about `0.71°`. The lower surface turns the flow by approximately `10.71°`, producing a much stronger pressure rise and a larger reduction in Mach number.

The analytical shock-angle prediction was then compared directly with the numerical solution.

The numerical shock location was extracted from maxima in the density-gradient field. A straight line was fitted through the detected ridge of each leading-edge shock, giving:

| Shock | Theory `β` | WENO5-JS + HLLC `β` | Difference |
|---|---:|---:|---:|
| Upper leading edge | 24.09° | 24.55° | +0.46° |
| Lower leading edge | 32.53° | 32.35° | -0.18° |

<p align="center">
  <img src="postprocessing/euler/figures/schlieren_shock_angles.png"
       alt="Leading-edge shock-angle verification from the WENO5-JS HLLC density-gradient field"
       width="900">
</p>

The extracted upper shock differs from theory by approximately `+0.46°`, while the lower shock differs by approximately `-0.18°`.

Both numerical shock angles therefore lie within **0.5°** of the nonlinear oblique-shock prediction.

The comparison is based on the resolved density-gradient ridge rather than a manually drawn line. The remaining difference is consistent with the finite numerical shock thickness, Cartesian-grid resolution and the finite spatial interval over which the numerical ridge is fitted.

The stronger lower shock was also checked using an off-body horizontal cut through the numerical solution at approximately `y/c = -0.304`.

<p align="center">
  <img src="postprocessing/euler/figures/lower_oblique_shock_profiles.png"
       alt="Lower oblique-shock Mach pressure and pressure-coefficient verification"
       width="900">
</p>

Across this cut, the numerical solution is compared with the ideal Rankine-Hugoniot post-shock state in terms of Mach number, pressure ratio and pressure coefficient.

The analytical shock is a mathematical discontinuity, while the finite-volume calculation captures the transition over a finite number of cells. Away from this captured shock thickness, the numerical solution approaches the theoretical post-shock state closely.

The leading-edge verification therefore checks both **where the shock forms** and **whether the state behind it is physically correct**:

```text
local flow turning
        ↓
theoretical shock angle
        ↓
numerically extracted shock angle
        ↓
post-shock Mach and pressure state
        ↓
surface-pressure loading
```

> **The WENO5-JS + HLLC solution recovers the upper and lower leading-edge shock angles to within 0.5° of nonlinear oblique-shock theory while also reproducing the expected post-shock states.**

### 5.2 Expansion-state verification

At the mid-chord corners, both surfaces turn through

```math
\Delta\theta = 2\delta \approx 11.42^\circ,
```

producing Prandtl-Meyer expansion fans.

Although the geometric turning is the same on both sides, the flow entering each expansion is different because of the unequal leading-edge shocks established in Section 5.1.

| Quantity | Upper surface | Lower surface |
|---|---:|---:|
| Mach before expansion | 2.470 | 2.056 |
| Expansion angle | 11.42° | 11.42° |
| Mach after expansion | 3.004 | 2.509 |
| Pressure ratio across expansion | 0.441 | 0.493 |
| Rear-panel pressure | 46.85 kPa | 96.95 kPa |

The upper expansion accelerates the flow to approximately **Mach 3.00**, producing the strongest pressure reduction in the solution. The lower rear-panel pressure remains much higher because its expansion begins from the strongly compressed lower-surface state.

These analytical states provide the reference for the numerical pressure loading examined next.

> **The solver must preserve the different states created by the leading-edge shocks as the flow passes through the two expansion fans.**

### 5.3 Surface-pressure verification

The analytical shock and expansion states give the following panel pressure coefficients:

| Panel | Analytical `Cp` |
|---|---:|
| Upper forward | +0.0110 |
| Lower forward | +0.2153 |
| Upper rear | -0.1229 |
| Lower rear | -0.0099 |

The numerical surface-pressure distribution is compared with both Ackeret theory and the nonlinear shock-expansion solution.

<p align="center">
  <img src="postprocessing/euler/figures/surface_cp.png"
       alt="Surface pressure coefficient verification for the Mach 2.5 diamond airfoil"
       width="900">
</p>

The WENO5-JS + HLLC solution reproduces the expected pressure-loading pattern: strong compression on the lower forward panel and the largest pressure reduction on the upper rear panel.

Because the airfoil is represented using a Cartesian immersed boundary, the numerical values are taken from the nearest external fluid cells rather than body-fitted wall faces. They should therefore be interpreted as a **nearest-fluid-cell pressure proxy**, particularly close to the sharp corners where the analytical solution changes discontinuously.

Away from these local corner regions, the numerical pressure levels follow the nonlinear shock-expansion prediction closely.

> **The resolved surface pressure provides the link between the local wave structure and the integrated lift and wave drag.**

### 5.4 Integrated-force verification

The final analytical check is whether the verified shock-expansion pressure field produces the correct integrated aerodynamic forces.

For the reference case, the final WENO5-JS + HLLC Euler solution gives

```math
C_D = 0.031606,
\qquad
C_L = 0.156778.
```

The numerical result is compared with both Ackeret linear theory and the nonlinear shock-expansion solution:

| Method | `Cd` | `Cl` |
|---|---:|---:|
| Ackeret linear theory | 0.030752 | 0.152345 |
| Nonlinear shock-expansion theory | 0.031714 | 0.156485 |
| WENO5-JS + HLLC Euler | **0.031606** | **0.156778** |

Relative to Ackeret linear theory,

```math
\Delta C_D
=
\frac{C_{D,\mathrm{CFD}}-C_{D,\mathrm{Ackeret}}}
{C_{D,\mathrm{Ackeret}}}\times100
\approx +2.78\%,
```

```math
\Delta C_L
=
\frac{C_{L,\mathrm{CFD}}-C_{L,\mathrm{Ackeret}}}
{C_{L,\mathrm{Ackeret}}}\times100
\approx +2.91\%.
```

Ackeret theory provides a useful linearised reference, but the present airfoil contains finite compression and expansion turning angles. The nonlinear shock-expansion solution is therefore the more appropriate reference for the final quantitative comparison.

Relative to nonlinear shock-expansion theory,

```math
\Delta C_D
=
\frac{0.03160637093-0.03171356478}
{0.03171356478}\times100
\approx -0.338\%,
```

```math
\Delta C_L
=
\frac{0.15677811494-0.15648527219}
{0.15648527219}\times100
\approx +0.187\%.
```

The WENO5-JS + HLLC solution therefore agrees with nonlinear shock-expansion theory to approximately **0.34% in drag** and **0.19% in lift**.

<p align="center">
  <img src="postprocessing/euler/figures/force_comparison.png"
       alt="Comparison of lift and drag coefficients from Ackeret theory, nonlinear shock-expansion theory and WENO5-JS HLLC Euler CFD"
       width="800">
</p>

This force agreement is not treated as an isolated verification result. The same numerical solution has already been examined through the complete compressible-flow chain:

- leading-edge shock angles;
- Rankine-Hugoniot post-shock states;
- Prandtl-Meyer expansion states;
- off-body shock structure;
- surface-pressure loading.

The relationship between these checks is

```text
shock and expansion geometry
            ↓
local thermodynamic states
            ↓
surface-pressure distribution
            ↓
integrated lift and wave drag
```

The integrated aerodynamic coefficients are therefore the final consequence of the resolved wave system and pressure loading rather than standalone numerical targets.

> **The WENO5-JS + HLLC solution reproduces the nonlinear shock-expansion prediction to approximately 0.34% in drag and 0.19% in lift, while independently recovering the wave structure and pressure states responsible for those forces.**
> 
## 6. Numerical convergence and robustness

Agreement with theory is only useful if the numerical solution itself is sufficiently settled.

For the final WENO5-JS + HLLC calculation, convergence was monitored using both the residual history and the integrated aerodynamic forces.

### 6.1 Convergence assessment

The final WENO5 JS + HLLC calculation was continued to **40,000 iterations**.

The integrated aerodynamic coefficients became stationary well before the end of the run:

```math
C_D = 0.03160637093,
\qquad
C_L = 0.15677811494.
```

<p align="center">
  <img src="postprocessing/euler/figures/force_history.png"
       alt="Lift and drag coefficient history for the WENO5 JS HLLC calculation"
       width="900">
</p>

The final 12 recorded force samples are identical at the precision stored in the force history:

| Quantity | Final value | Range over final 12 samples |
|---|---:|---:|
| `Cd` | 0.03160637093 | 0 |
| `Cl` | 0.15677811494 | 0 |

The cellwise residual and solution change histories remain oscillatory at late iterations and do not provide standalone evidence of convergence. Their history is therefore retained in `figures/diagnostics/` as a numerical diagnostic rather than used as the primary convergence result.

Acceptance of the Euler reference is based on the stationary aerodynamic coefficients together with the independently verified shock geometry, pressure states and surface loading.

> **The final aerodynamic loading is stationary, while the residual history is reported separately and interpreted as a diagnostic rather than a convergence claim.**

### 6.2 Robustness to reconstruction method

Before the final WENO5-JS formulation, the HLLC solver was tested with Barth-Jespersen limited reconstruction.

Both reconstruction approaches produced essentially the same integrated aerodynamic loading, while WENO5-JS resolved the local shock-expansion structure more sharply.

| Solver branch | Main role | Integrated loading |
|---|---|---|
| HLLC + Barth-Jespersen | Intermediate verification branch | Consistent with final solution |
| WENO5-JS + HLLC | Final reference solver | `Cd = 0.031606`, `Cl = 0.156778` |

This provides a useful robustness check: changing the reconstruction altered the local resolution of discontinuities without changing the aerodynamic conclusion.

> **WENO5-JS was retained because it improved wave resolution, not because it produced a different lift or drag result.**
### 6.3 Sharp edge and immersed-boundary limitation

The diamond airfoil has ideal zero-radius leading and trailing edges.

On the Cartesian immersed-boundary grid, these corners cannot be represented as exact body-fitted surface points. The local pressure and gradient fields very close to the corners are therefore more sensitive to the discrete body representation than the flow away from the surface.

For this reason, isolated values taken directly at a sharp corner are not used as standalone evidence. The solution is assessed using the surrounding shock and expansion structure, surface-pressure behaviour away from the vertices, analytical state checks and integrated aerodynamic forces.

## 7. Increasing the modelling fidelity: wall-resolved SA-RANS

The Euler branch establishes the inviscid shock-expansion reference, but it cannot represent skin friction, boundary-layer development or viscous interaction with the pressure field.

The same Mach 2.5, `α = 5°` diamond configuration is therefore being investigated using compressible wall-resolved RANS in OpenFOAM with the Spalart-Allmaras turbulence model [7].

The investigation was extended beyond solver convergence alone. Domain sensitivity, mass conservation, numerical settings and mesh quality were examined separately, which ultimately led to a redesign of the RANS meshing strategy before any production result was accepted.

> **The RANS stage is treated as a verification problem in its own right rather than as a direct extension of the Euler solution.**

### 7.1 Initial RANS assessment

The first wall-resolved SA-RANS calculation used a structured body-fitted H-grid.

The near-wall treatment was strong, with an average `y+ ≈ 0.33` and `P99(y+) ≈ 0.57`. However, continuation beyond `34,000` iterations showed that the viscous part of the solution had not reached a stationary state.

The pressure contribution had largely stabilised, but viscous drag and wall shear continued to evolve beyond the adopted convergence limits.

The calculation was therefore not accepted as the final RANS reference.

This triggered a wider investigation of the numerical setup rather than simply extending the run further.

### 7.2 Why the original H-grid was abandoned

The convergence issue was not traced to domain size or mass conservation, so the mesh itself was examined in more detail.

The structured H-grid showed a broad determinant-quality problem rather than a small number of isolated bad cells. More than `23%` of the D1 medium grid fell below the adopted determinant threshold.

The low-quality cells followed the global structured mapping from the wall towards the far field, showing that the near-wall resolution was being propagated too aggressively through the complete domain.

Two further structured `blockMesh` redesigns were tested, but both retained or worsened the same underlying problem.

The H-grid strategy was therefore abandoned rather than refined further.

> **The problem was treated as a topology issue, not simply as a lack of cells.**

### 7.3 Transition to a hybrid Gmsh topology

After the structured redesigns failed, the meshing strategy was changed rather than patched again.

The replacement approach uses a hybrid topology so that each region of the flow can be resolved according to its own numerical requirement:

```text
structured wall-resolved strips
        +
explicit sharp-edge treatment
        +
unstructured outer field
        +
local wake refinement
        +
shock-aligned refinement corridors
```
This removes the need for the near-wall spacing to control the complete far-field mesh.

The outer transition, wake region and shock-refinement strategy were developed and checked separately before being combined with the wall treatment.

> **The new mesh is built around the physics of the problem rather than a single global structured mapping.**

### 7.4 Wall-resolved strip and sharp-edge development

The new wall treatment was developed separately from the outer mesh so that the near-wall requirements could be verified directly.

The revised standard wall strip uses:

```text
first wall spacing   = 1.2e-6 m
tangential cells     = 5500 per panel
wall-normal cells    = 90
normal progression   = 1.0743
```
The first-cell height was increased only after the previous RANS solution showed sufficient y+ margin to remain wall resolved.

Local OpenFOAM tests confirmed that the interior wall strip could satisfy the required quality criteria without reproducing the earlier global H-grid problem.

The exact sharp leading and trailing edges were then treated separately because the corner interface introduced different numerical constraints from the normal wall region.

Several local corner variants were rejected before the final V4 treatment removed the interpolation-weight and volume-ratio failures at the actual sharp-edge cap.

> **The wall strip and sharp-edge regions were therefore designed and verified as separate mesh problems rather than forcing one spacing strategy to satisfy both.**

### 7.5 Current RANS status

The accepted wall, corner, wake, shock and outer-field strategies are now being integrated into the full replacement mesh.

Local tests have already been used to remove the dominant wall-strip and sharp-edge quality problems before committing to another production CFD run.

The complete mesh will only be promoted to the new RANS reference after the full OpenFOAM audit confirms:

```text
one connected region
correct physical patches
frontAndBack = empty
zero determinant failures
zero low interpolation-weight failures
zero low volume-ratio failures
acceptable non-orthogonality and skewness
verified wall spacing and wall-face count
```
Only after this mesh passes the global quality gate will the new wall-resolved SA-RANS calculation begin.

> **No aerodynamic result from the redesigned mesh will be accepted before both the mesh and the resulting RANS solution pass their respective verification gates.**

## 8. Verification path to aerodynamic design

The redesigned RANS methodology will not move directly into a design study once the mesh is generated.

The replacement hybrid mesh must first pass the complete OpenFOAM quality audit and then produce a stationary wall-resolved SA-RANS reference solution under the locked convergence criteria.

Only after that reference has been established will a formal three-grid convergence study be carried out.

The intended progression is:

```text
integrated hybrid mesh certification
        ↓
converged RANS reference
        ↓
G3 / G2 / G1 grid-convergence study
        ↓
Richardson extrapolation and GCI
        ↓
freeze production CFD methodology
        ↓
controlled angle-of-attack study
        ↓
aerodynamic design and optimisation
```
The objective is to ensure that later aerodynamic trends are produced by changes in operating condition or geometry rather than by unresolved numerical uncertainty.
> **The design study begins only after the RANS mesh, convergence behaviour and grid sensitivity have been demonstrated independently.**

### 8.1 Formal grid convergence and GCI

Once the redesigned RANS reference case has passed the mesh-quality and convergence gates, a new three-grid family will be generated from the corrected topology.

The coarse, medium and fine meshes will retain the same:

- physical domain;
- wall treatment;
- sharp-edge strategy;
- wake and shock refinement logic;
- solver settings;
- force and moment definitions.

The objective is to quantify numerical uncertainty rather than simply show that the coefficients change only slightly with refinement.

The study will use Richardson extrapolation and the Grid Convergence Index to assess quantities such as

```math
C_D,\qquad
C_{D,p},\qquad
C_{D,v},\qquad
C_L,\qquad
C_m.
```
Surface pressure and wall quantities will also be checked to confirm that global force agreement is supported by consistent local flow behaviour.

> **Grid independence will only be claimed from a mesh family that shares the corrected topology and demonstrates systematic convergence.**

### 8.2 Freeze the production CFD methodology

After the corrected mesh family has passed the grid-convergence study, the accepted numerical setup will be frozen before any design changes are introduced.

The production methodology will retain the same:

- computational domain;
- meshing strategy;
- wall resolution;
- sharp-edge treatment;
- turbulence model;
- numerical schemes;
- boundary conditions;
- force and moment definitions;
- convergence criteria;
- post-processing procedure.

This prevents changes in aerodynamic performance from being confused with changes in the CFD setup.

> **Once the numerical methodology has been verified, the solver and mesh setup stop being additional design variables.**

### 8.3 Controlled angle-of-attack study

The first production study will be a controlled angle-of-attack sweep using the frozen RANS methodology.

The purpose is to track how incidence changes the complete shock, pressure and boundary-layer system rather than looking only at the final force coefficients.

For each angle of attack, the study will examine:

- leading-edge shock angle and strength;
- expansion behaviour;
- surface-pressure redistribution;
- viscous drag and wall shear;
- lift, drag and pitching moment;
- changes in aerodynamic efficiency.

The same mesh strategy, solver settings and convergence criteria will be retained across the sweep so that the observed trends can be attributed to the operating condition.

> **The angle-of-attack study will be used to connect changes in the wave structure and viscous response directly to the resulting aerodynamic performance.**

### 8.4 Aerodynamic design and optimisation

Once the angle-of-attack study has established the baseline aerodynamic trends, the same verified RANS methodology will be used for controlled geometry changes.

The objective will be to improve aerodynamic performance without losing sight of the underlying flow physics.

Candidate designs will therefore be assessed through both the integrated coefficients and the changes responsible for them:

- shock strength and position;
- expansion behaviour;
- surface-pressure distribution;
- viscous drag and wall shear;
- lift-to-drag ratio;
- pitching-moment behaviour.

Any geometry change will be compared against the same frozen reference methodology so that the effect of the design itself can be separated from numerical variation.

> **The optimisation stage will only begin after the numerical uncertainty of the baseline method has been established and the main aerodynamic mechanisms are understood.**

---

## References

[1] J. D. Anderson Jr., *Modern Compressible Flow: With Historical Perspective*, 3rd ed., McGraw-Hill, 2003.

[2] E. F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics: A Practical Introduction*, 3rd ed., Springer, 2009.

[3] E. F. Toro, M. Spruce and W. Speares, “Restoration of the contact surface in the HLL-Riemann solver,” *Shock Waves*, Vol. 4, pp. 25–34, 1994. https://doi.org/10.1007/BF01414629

[4] G.-S. Jiang and C.-W. Shu, “Efficient implementation of weighted ENO schemes,” *Journal of Computational Physics*, Vol. 126, No. 1, pp. 202–228, 1996. https://doi.org/10.1006/jcph.1996.0130

[5] C.-W. Shu and S. Osher, “Efficient implementation of essentially non-oscillatory shock-capturing schemes,” *Journal of Computational Physics*, Vol. 77, No. 2, pp. 439–471, 1988. https://doi.org/10.1016/0021-9991(88)90177-5

[6] T. J. Barth and D. C. Jespersen, “The design and application of upwind schemes on unstructured meshes,” AIAA Paper 89-0366, 27th Aerospace Sciences Meeting, 1989. https://doi.org/10.2514/6.1989-366

[7] P. R. Spalart and S. R. Allmaras, “A one-equation turbulence model for aerodynamic flows,” AIAA Paper 92-0439, 30th Aerospace Sciences Meeting and Exhibit, 1992. https://doi.org/10.2514/6.1992-439

[8] R. Mittal and G. Iaccarino, “Immersed Boundary Methods,” *Annual Review of Fluid Mechanics*, Vol. 37, pp. 239–261, 2005. https://doi.org/10.1146/annurev.fluid.37.061903.175743
