# Euler solver source layout

The original single-file solver has been split by numerical role without changing the algorithms.

- `types.hpp` — data structures and case configuration
- `flow.hpp` — flow-state conversions, thermodynamics and analytical helpers
- `geometry.hpp` — diamond geometry and wall geometry helpers
- `fluxes.hpp` — Euler fluxes, HLLC, LLF fallback and WENO5-JS reconstruction
- `immersed_boundary.hpp` — immersed-boundary face classification and wall normals
- `solver.hpp` — residual assembly, pseudo-time step, SSP-RK3 support and far-field damping
- `forces.hpp` — pressure sampling and aerodynamic force integration
- `io.hpp` — CSV output, convergence checks and command-line parsing
- `main.cpp` — case setup and solver loop

Build from this folder with:

```bash
g++ -std=c++17 -O3 main.cpp -o euler_solver
```
