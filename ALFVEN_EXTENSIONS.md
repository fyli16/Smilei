# Alfvén-wave study extensions to Smilei

This document summarizes two custom additions made to Smilei for Alfvén-wave
instability studies with the standard **Yee** Maxwell solver (1D and 2D Cartesian):

1. A **`"damping"` electromagnetic boundary condition** — an interior absorbing
   "mask"/sponge layer of finite thickness at both ends of the box along `x`.
2. A **`FieldInjector` module** — an additive "soft source" that injects an `E` or
   `B` field profile at a fixed location *inside* the plasma, launching a propagating
   wave (with a customizable transverse profile in 2D).

Both were implemented on the branch `feature/alfven-damping-mask`.

---

## 1. `"damping"` boundary condition (absorbing mask)

### Purpose
Place two masks of thickness `x_mask` inside the box, `[0, x_mask]` and
`[xmax - x_mask, xmax]`, that progressively damp the **transverse wave fields**
(`Ey, Ez, By, Bz`) so waves are absorbed before reaching the box wall. The uniform
guide field `Bx` (`B0 ∥ x`) and the parallel field `Ex` are **left untouched**, so the
mean magnetic field and the electrostatic plasma response are preserved.

### Behavior
Inside a mask, every field component that is damped is multiplied each time step by a
smooth parabolic coefficient
```
coeff(x) = 1 - EM_damping_coefficient * s(x)^2
```
where `s` is the normalized depth into the mask: `s = 0` at the inner edge (no damping,
`coeff = 1`) and `s = 1` at the box wall (`coeff = 1 - EM_damping_coefficient`). The mask
depends only on `x`, so in 2D the coefficient is uniform along `y`. The coefficient is
computed from each cell's **global** `x`, so a mask may span several MPI patches.

### Namelist parameters (`Main` block)

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `EM_boundary_conditions` | list | — | use `"damping"` on the x-min / x-max sides |
| `EM_damping_thickness` | float | `0.` | mask thickness `x_mask` (code length units); must be > 0 |
| `EM_damping_coefficient` | float | `1.` | max damping strength in `[0,1]` (`1` = fields fully removed at the wall) |

Constraints (validated at runtime): thickness > 0, coefficient in `[0,1]`, the two masks
must not overlap (`2 * thickness < xmax`), and `"damping"` is currently supported only
along `x`.

### Files
- New BC classes: `src/ElectroMagnBC/ElectroMagnBC1D_Damping.{h,cpp}`,
  `src/ElectroMagnBC/ElectroMagnBC2D_Damping.{h,cpp}`
- Registered in `src/ElectroMagnBC/ElectroMagnBC_Factory.h`
- Parameters in `src/Params/Params.{h,cpp}`, defaults in `src/Python/pyinit.py`
- OpenPMD metadata recognition in `src/Params/OpenPMDparams.cpp`
- Documentation in `doc/Sphinx/Use/namelist.rst`

### Notes
- `"damping"` is not an "open" boundary in the Silver–Müller sense; the box wall behind
  the mask is effectively closed, but fields are ~0 there thanks to the mask.
- You still set the **particle** boundary conditions separately (e.g. periodic or
  thermalizing) in each `Species`.

---

## 2. `FieldInjector` (additive soft field source)

### Purpose
Inject an Alfvén wave by directly adding an `E` or `B` field profile at a fixed `x`
plane inside the plasma, rather than driving a current with an `Antenna`. Because the
source is **additive and not reset**, the injected field enters Maxwell's equations and
**propagates self-consistently** as a wave. Placing the plane just to the right of the
left mask means the left-going branch is absorbed quickly while the right-going wave is
used for the physics.

### Behavior
Each time step the chosen field (one of `Ex, Ey, Ez, Bx, By, Bz`) is updated as
```
field(x, y, t) += profile(x, y, t)      # evaluated at the current simulation time
```
just before the Maxwell solve. Injection targets the solver-evolved field (e.g. `Bz`,
not the time-centered `Bz_m`), so it both propagates and is felt by the particle pusher.
In 2D/3D the profile depends on the transverse coordinates, so the field pattern along
`y` (and `z`) at the injection plane is fully customizable.

### Namelist block

```python
FieldInjector(
    field   = "Bz",                      # one of Ex, Ey, Ez, Bx, By, Bz
    profile = lambda x, y, t: ...,       # (x,t) in 1D, (x,y,t) in 2D, (x,y,z,t) in 3D
)
```

| Parameter | Type | Meaning |
|-----------|------|---------|
| `field` | string | injected field: `Ex`, `Ey`, `Ez`, `Bx`, `By`, `Bz` |
| `profile` | profile | space–time profile added to the field each step |

### Files
- Struct + apply method: `src/ElectroMagn/ElectroMagn.{h,cpp}`
  (`applyFieldInjectors`, reusing the additive `applyPrescribedField`, `mode=3`)
- Namelist parsing: `src/ElectroMagn/ElectroMagnFactory.h`
- Loop wiring: `src/Patch/VectorPatch.{h,cpp}`, `src/Smilei.cpp`
  (applied with the antennas, inside the OpenMP region, before `solveMaxwell`)
- Python class + registration: `src/Python/pyinit.py`, `src/Python/pycontrol.py`
- Documentation in `doc/Sphinx/Use/namelist.rst`

### Notes
- It is a **soft source**: the effective radiated amplitude depends on the source
  amplitude and the time step — calibrate empirically (same behavior as an `Antenna`).
- For a clean Alfvén wave, inject the transverse magnetic perturbation (`Bz` and/or `By`)
  and, if desired, the transverse electric field with the Alfvénic polarization
  `δE = -v_A · δB`.
- Requires a Yee-type solver; not compatible with spectral / multiple-decomposition
  solvers (guarded at runtime).

---

## 3. Example namelist — 1D

A uniform magnetized plasma with a guide field `B0 x̂`, absorbing masks at both ends,
and a `Bz` Alfvén driver just right of the left mask.

```python
import math

# --- physical / numerical parameters (normalized units) ---
B0      = 0.1                     # uniform guide field along x
x_mask  = 20.                     # mask thickness at each end
Lx      = 200.                    # box length
dx      = 0.1
dt      = 0.95 * dx               # CFL
t_sim   = 400.

omega_inj = 0.05                  # driver angular frequency
x_inj     = x_mask + 5.           # injection plane, just right of the left mask
w_inj     = 1.0                   # driver spatial half-width
dB        = 1e-3 * B0             # driver amplitude

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,
    cell_length = [dx],
    grid_length = [Lx],
    number_of_patches = [16],
    timestep = dt,
    simulation_time = t_sim,
    # absorbing masks at x-min and x-max:
    EM_boundary_conditions = [["damping", "damping"]],
    EM_damping_thickness   = x_mask,
    EM_damping_coefficient = 1.0,
)

# uniform guide field B0 along x
ExternalField(field="Bx", profile=B0)

# background magnetized plasma (ions + electrons)
def n0(x):
    return 1.0

for name, mass, charge in [("ions", 100., 1.), ("electrons", 1., -1.)]:
    Species(
        name = name,
        position_initialization = "regular",
        momentum_initialization = "maxwell-juettner",
        particles_per_cell = 64,
        mass = mass,
        charge = charge,
        number_density = n0,
        temperature = [1e-4]*3,
        boundary_conditions = [["periodic"]],   # particle BCs (independent of field BCs)
    )

# Alfvén wave soft source: transverse B perturbation at a fixed x
FieldInjector(
    field = "Bz",
    profile = lambda x, t: dB * math.exp(-((x - x_inj)**2)/(2*w_inj**2)) * math.sin(omega_inj*t),
)

DiagFields(every=200, fields=["Ey","Ez","By","Bz","Bx"])
```

---

## 4. Example namelist — 2D

Same idea in 2D, with a transverse (`y`) profile at the injection plane. Masks act along
`x`; `y` is periodic.

```python
import numpy as np

B0      = 0.1
x_mask  = 20.
Lx, Ly  = 200., 50.
dx = dy = 0.2
dt      = 0.95 / np.sqrt(1/dx**2 + 1/dy**2)   # CFL
t_sim   = 400.

omega_inj = 0.05
x_inj     = x_mask + 5.
w_inj     = 1.0
ky        = 2*np.pi / Ly          # one transverse wavelength across the box
dB        = 1e-3 * B0

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    cell_length = [dx, dy],
    grid_length = [Lx, Ly],
    number_of_patches = [16, 8],
    timestep = dt,
    simulation_time = t_sim,
    # masks along x, periodic along y:
    EM_boundary_conditions = [ ["damping", "damping"], ["periodic", "periodic"] ],
    EM_damping_thickness   = x_mask,
    EM_damping_coefficient = 1.0,
)

# uniform guide field B0 along x
ExternalField(field="Bx", profile=B0)

for name, mass, charge in [("ions", 100., 1.), ("electrons", 1., -1.)]:
    Species(
        name = name,
        position_initialization = "regular",
        momentum_initialization = "maxwell-juettner",
        particles_per_cell = 32,
        mass = mass,
        charge = charge,
        number_density = 1.0,
        temperature = [1e-4]*3,
        boundary_conditions = [ ["periodic"], ["periodic"] ],
    )

# Alfvén wave soft source with a custom transverse (y) profile at a fixed x plane
def inject_Bz(x, y, t):
    envelope_x = np.exp(-((x - x_inj)**2)/(2*w_inj**2))
    return dB * envelope_x * np.sin(omega_inj*t - ky*y)

FieldInjector(field="Bz", profile=inject_Bz)

DiagFields(every=200, fields=["Ey","Ez","By","Bz","Bx"])
```

> In 2D you can shape the field along `y` arbitrarily via the injector `profile`
> (here a single transverse wavelength `sin(ω t − k_y y)`; replace with any function of
> `y` for a custom transverse structure). Add a second `FieldInjector` for the transverse
> electric field if you want to impose the full Alfvénic polarization.

---

## 5. Quick reference

| Feature | Namelist entry | Key knobs |
|---------|----------------|-----------|
| Absorbing masks | `EM_boundary_conditions = [["damping","damping"], ...]` | `EM_damping_thickness`, `EM_damping_coefficient` |
| Guide field | `ExternalField(field="Bx", profile=B0)` | `B0` |
| Wave injection | `FieldInjector(field="Bz", profile=...)` | injection plane `x_inj`, `omega_inj`, amplitude, transverse profile |

Both features assume the **Yee** solver (`maxwell_solver = "Yee"`, the default) and were
tested for 1D and 2D Cartesian geometry.
