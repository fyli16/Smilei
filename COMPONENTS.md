# Smilei — Codebase Components & Implementation Guide

This document summarizes the main components of the [Smilei](https://smileipic.github.io/Smilei/)
particle-in-cell (PIC) code as laid out in this repository, and then investigates **where to
implement two custom features** needed for our Alfvén-wave instability studies:

1. **Absorbing "mask" damping regions** of finite width at both ends of the box (x-min / x-max),
   used to artificially damp the waves and mimic absorbing boundary conditions.
2. **Alfvén-wave injection** via an oscillating current at a fixed `x` location inside the plasma,
   placed close to the left mask so that the left-going wave is quickly absorbed and the
   right-going wave is used for the physics study.

Target: **1D and 2D Cartesian geometry** with the standard **Yee** FDTD Maxwell solver.

---

## 1. Top-level layout

| Path | Purpose |
|------|---------|
| [makefile](makefile) | Build system entry point |
| [smilei.sh](smilei.sh) | Convenience launch script |
| [src/](src) | All C++ source code (the simulation engine) |
| [happi/](happi) | Python post-processing / analysis package |
| [benchmarks/](benchmarks) | Reference input namelists (`tst*.py`) used for validation |
| [validation/](validation) | Automated non-regression test harness |
| [doc/](doc) | Sphinx documentation sources |
| [tools/](tools) | Auxiliary tools |
| [scripts/](scripts) | Helper scripts |

Input files ("namelists") are **Python scripts**. The C++ engine embeds a Python interpreter and
reads simulation parameters from named blocks (e.g. `Main(...)`, `Species(...)`, `Antenna(...)`).

---

## 2. Core C++ components (`src/`)

### Simulation driver
- [src/Smilei.cpp](src/Smilei.cpp) — `main()` and the **top-level PIC time loop**. Orchestrates
  particle push, current deposition, antennas, Maxwell solve, boundary conditions, diagnostics.
- [src/Params/](src/Params) — parses and stores all global parameters read from the namelist
  ([Params.cpp](src/Params/Params.cpp)).
- [src/Python/](src/Python) — the embedded-Python layer. `pyinit.py` defines the namelist blocks
  (including the `Antenna` class), `pyprofiles.py` defines profile helpers, `PyTools` bridges
  Python values into C++.

### Fields (electromagnetics)
- [src/Field/](src/Field) — the low-level grid containers (`Field1D`, `Field2D`, `Field3D`).
- [src/ElectroMagn/](src/ElectroMagn) — the EM field manager. `ElectroMagn` (base) plus
  `ElectroMagn1D` / `ElectroMagn2D` / `ElectroMagn3D` / `ElectroMagnAM` hold `Ex..Bz`, currents
  `Jx/Jy/Jz`, and implement `saveMagneticFields`, `centerMagneticFields`, **`applyAntenna`**,
  `initAntennas`, `applyExternalFields`, and **`boundaryConditions`**.
  `Laser.cpp` handles laser injection at boundaries.
- [src/ElectroMagnSolver/](src/ElectroMagnSolver) — the FDTD Maxwell solvers. Each timestep splits
  into a **Maxwell–Ampère** update (E) and a **Maxwell–Faraday** update (B):
  - `MA_Solver{1,2,3}D_norm` — Ampère (E) update, standard.
  - **`MF_Solver1D_Yee`, `MF_Solver2D_Yee`, `MF_Solver3D_Yee`** — the **Yee** Faraday (B) update
    (this is the solver we use).
  - Other variants: `M4`, `Lehe`, `Bouchard`, `Cowan`, `Grassi`, `Friedman`, `PML_*` (perfectly
    matched layers), `PXR_*` (spectral/PICSAR).
  - `SolverFactory.h` selects the solver from the namelist.
- [src/ElectroMagnBC/](src/ElectroMagnBC) — **field boundary conditions**, one object per box side.
  Registered in `ElectroMagnBC_Factory.h`. Existing types:
  - `*_SM` — Silver–Müller (absorbing / injecting).
  - `*_refl` — reflective.
  - `*_PML` — perfectly matched layers.
  - **`ElectroMagnBC2D_Trans_Damping`** — a **damping-layer** BC that multiplies the fields by a
    position-dependent coefficient inside a boundary layer. **This is the closest existing model
    to our mask damping regions** and the best template to copy.

### Particles
- [src/Particles/](src/Particles) — particle data containers.
- [src/Species/](src/Species) — a plasma species (`Species.cpp`), its initialization, and the
  per-species loop ordering.
- [src/Pusher/](src/Pusher) — particle equations of motion (Boris, Vay, Higuera–Cary…).
- [src/Interpolator/](src/Interpolator) — gathers fields at particle positions.
- [src/Projector/](src/Projector) — deposits charge/current onto the grid (charge-conserving).
- [src/ParticleBC/](src/ParticleBC) — particle boundary conditions.
- [src/ParticleInjector/](src/ParticleInjector) — continuous particle injection from boundaries.

### Additional physics modules
- [src/Collisions/](src/Collisions), [src/Ionization/](src/Ionization),
  [src/Radiation/](src/Radiation), [src/MultiphotonBreitWheeler/](src/MultiphotonBreitWheeler),
  [src/Merging/](src/Merging) — binary collisions, field/collisional ionization, radiation
  reaction, QED pair production, particle merging.

### Parallelization & domain
- [src/Patch/](src/Patch) — domain decomposition into **patches**. `VectorPatch.cpp` is the key
  orchestration layer: it holds every patch and drives `solveMaxwell`, **`applyAntennas`**,
  `finalizeSyncAndBCFields`, diagnostics, and MPI synchronization.
- [src/DomainDecomposition/](src/DomainDecomposition), [src/SmileiMPI/](src/SmileiMPI) — MPI layout
  and communication.
- [src/MovWindow/](src/MovWindow) — moving simulation window.

### Support
- [src/Diagnostic/](src/Diagnostic) — all diagnostics (fields, scalars, probes, particle binning…).
- [src/Checkpoint/](src/Checkpoint) — dump / restart.
- [src/Profiles/](src/Profiles) — space/time profile evaluation used by species, lasers, antennas.
- [src/Tools/](src/Tools) — utilities (error handling, timers, `TITLE`, etc.).

---

## 3. The field time loop (where our changes plug in)

In [src/Smilei.cpp](src/Smilei.cpp) the relevant per-timestep ordering is:

```
dynamics()          // push particles, deposit J
sumDensities()      // assemble total Jx/Jy/Jz on the grid
applyAntennas()     // <-- external currents (antenna) ADDED to J   [Smilei.cpp ~L533]
solveMaxwell()      // Ampère (E) then Faraday (B)                  [VectorPatch.cpp ~L965]
finalizeSyncAndBCFields()  // EMfields->boundaryConditions()        [VectorPatch.cpp ~L1167]
```

Key call chain details:
- **Antennas**: [VectorPatch::applyAntennas](src/Patch/VectorPatch.cpp) →
  [ElectroMagn::applyAntenna](src/ElectroMagn/ElectroMagn.cpp) adds
  `intensity * antennaField(i)` into the current field each step. Time dependence comes from the
  antenna's `time_profile`, spatial shape from its `space_profile`.
- **Maxwell**: [VectorPatch::solveMaxwell](src/Patch/VectorPatch.cpp) calls
  `MaxwellAmpereSolver_` then `MaxwellFaradaySolver_` (our Yee solvers).
- **Boundary conditions**: [ElectroMagn::boundaryConditions](src/ElectroMagn/ElectroMagn.cpp)
  loops over `emBoundCond[]` and calls each `apply()`. For 1D there are 2 sides (x-min, x-max);
  for 2D there are 4 (x-min, x-max, y-min, y-max). Objects are created in
  [ElectroMagnBC_Factory.h](src/ElectroMagnBC/ElectroMagnBC_Factory.h).

---

## 4. Feature 1 — Alfvén-wave injection (fixed-x current antenna)

**Recommendation: no C++ change is required — use an `Antenna` block in the namelist.**

An antenna is an extra current `J` added to the grid every timestep, at whatever spatial location
we choose, with any time dependence. This is exactly a localized wave driver. A current sheet at a
fixed `x = x0` radiates waves in **both** ±x directions; by placing `x0` just to the right of the
left mask, the left-going wave is absorbed almost immediately and the right-going wave survives —
precisely the setup we want.

For an Alfvén wave (transverse perturbation of `B` and velocity), drive a **transverse current**
(`Jy` and/or `Jz`), oscillating at the desired frequency `ω`, localized around `x0`:

```python
import numpy as np

x0     = 5.0          # injection location (just right of the left mask)
w      = 0.3          # spatial half-width of the driver
omega  = ...          # Alfvén wave angular frequency
J0     = ...          # driver amplitude

def envelope_x(x):
    return np.exp(-((x - x0)**2) / (2*w**2))

# Option A: separate space and time profiles (scalar time profile)
Antenna(
    field         = "Jz",
    space_profile = lambda x: J0 * envelope_x(x),      # 1D: f(x); 2D: f(x,y)
    time_profile  = lambda t: np.sin(omega * t),
)

# Option B: a single space-time profile for a true traveling-wave phase in 2D
# Antenna(
#     field = "Jz",
#     space_time_profile = lambda x, y, t: J0*envelope_x(x)*np.sin(omega*t - ky*y),
# )
```

- In **2D** the `space_profile` takes `(x, y)` and the `space_time_profile` takes `(x, y, t)`.
- Use `Jy`/`Jz` (transverse) for shear Alfvén-type perturbations; the polarization depends on the
  background `B0` orientation in your setup.
- Reference for the interface: [doc/Sphinx/Use/namelist.rst](doc/Sphinx/Use/namelist.rst)
  (the *Antenna* section), Python class in [src/Python/pyinit.py](src/Python/pyinit.py), parsing in
  [src/ElectroMagn/ElectroMagnFactory.h](src/ElectroMagn/ElectroMagnFactory.h).

If a namelist-only antenna proves insufficient (e.g. you need a directional/one-way launcher that
cancels the left-going branch, or injection tied to particle velocities), the C++ hook to extend is
[ElectroMagn::applyAntenna](src/ElectroMagn/ElectroMagn.cpp) together with
[VectorPatch::applyAntennas](src/Patch/VectorPatch.cpp).

---

## 5. Feature 2 — Mask damping regions at both x ends

**Recommendation: implement a new field boundary-condition class that damps a finite-width layer at
x-min and x-max, modeled on [ElectroMagnBC2D_Trans_Damping](src/ElectroMagnBC/ElectroMagnBC2D_Trans_Damping.cpp).**

That existing class already demonstrates the exact technique we want — multiply `E` and `B` by a
smooth position-dependent coefficient `coeff[j]` (0 at the wall, →1 at the inner edge) inside a
layer of `N` cells — but it operates on the **y** boundary and its width/strength are hardcoded.
We need the same along **x**, for both 1D and 2D, with width and damping strength read from the
namelist.

### Where to add code

1. **New BC classes** (mirror the naming used by the factory):
   - `src/ElectroMagnBC/ElectroMagnBC1D_Damping.{h,cpp}`
   - `src/ElectroMagnBC/ElectroMagnBC2D_Damping.{h,cpp}`

   Each derives from `ElectroMagnBC1D` / `ElectroMagnBC2D`, precomputes the damping profile
   `coeff[i]` over the mask width in the constructor, and in `apply()` multiplies the field arrays
   by `coeff` inside the mask. Note the Yee grid staggering: apply the coefficient to every field
   component present (`Ex,Ey,Ez,Bx,By,Bz` in 2D; `Ey,Ez,By,Bz` in 1D — `Ex`/`Bx` are trivial in
   1D), indexing primal vs dual arrays as done in the Yee solvers
   ([MF_Solver1D_Yee.cpp](src/ElectroMagnSolver/MF_Solver1D_Yee.cpp),
   [MF_Solver2D_Yee.cpp](src/ElectroMagnSolver/MF_Solver2D_Yee.cpp)).
   Use `patch->isXmin()` / `patch->isXmax()` to act only on the boundary patches, exactly like
   `isYmin()` / `isYmax()` in the reference class.

2. **Register them** in
   [src/ElectroMagnBC/ElectroMagnBC_Factory.h](src/ElectroMagnBC/ElectroMagnBC_Factory.h): add a
   new option (e.g. `"damping"`) in the `1Dcartesian` and `2Dcartesian` branches, so
   `EM_boundary_conditions = [["damping","damping"], ...]` selects it for the x sides.

3. **Namelist parameters** for mask width and damping coefficient. Add fields to
   [src/Params/Params.cpp](src/Params/Params.cpp) / `Params.h` (read via `PyTools`) so the width
   (`ny_l`-equivalent) and strength (`cdamp`) are configurable rather than hardcoded as they are in
   the Trans_Damping reference.

### How it runs
The new `apply()` is invoked automatically through
[ElectroMagn::boundaryConditions](src/ElectroMagn/ElectroMagn.cpp), which is called each step from
[VectorPatch::finalizeSyncAndBCFields](src/Patch/VectorPatch.cpp) (single-decomposition Cartesian
path) after the Yee Maxwell solve. No change to the solver or main loop is needed.

### Alternatives considered
- **PML** (`ElectroMagnBC2D_PML`) is a physically rigorous absorber and is already available for
  x-boundaries in 2D via `EM_boundary_conditions = [["PML", ...]]`. For a plain absorbing boundary
  it may be simpler than writing new code. However, our requirement is an explicit *finite-width
  interior mask* whose width/strength we tune for Alfvén-wave damping, so the custom damping-layer
  class gives the direct control we want and matches the "two mask regions" design.
- A **volumetric mask applied outside the BC system** (a multiplicative array applied to `E`/`B`
  right after the Faraday update in [VectorPatch::solveMaxwell](src/Patch/VectorPatch.cpp)) is also
  possible and would decouple the masks from the box edges, but it is more invasive; the BC-class
  route reuses existing infrastructure and is the recommended first implementation.

---

## 6. Summary of touch-points

| Task | Change type | Files |
|------|-------------|-------|
| Alfvén-wave injection at fixed x | **Namelist only** | `Antenna(...)` block in the input `.py` |
| (optional) advanced injection | C++ | [ElectroMagn::applyAntenna](src/ElectroMagn/ElectroMagn.cpp), [VectorPatch::applyAntennas](src/Patch/VectorPatch.cpp) |
| Mask damping (1D/2D) | **New C++ BC + factory + params** | new `ElectroMagnBC{1,2}D_Damping.{h,cpp}`, [ElectroMagnBC_Factory.h](src/ElectroMagnBC/ElectroMagnBC_Factory.h), [Params.cpp](src/Params/Params.cpp) |
| Reference to copy for damping | — | [ElectroMagnBC2D_Trans_Damping.cpp](src/ElectroMagnBC/ElectroMagnBC2D_Trans_Damping.cpp) |
| Yee solver (unchanged, for staggering reference) | — | [MF_Solver1D_Yee.cpp](src/ElectroMagnSolver/MF_Solver1D_Yee.cpp), [MF_Solver2D_Yee.cpp](src/ElectroMagnSolver/MF_Solver2D_Yee.cpp) |
