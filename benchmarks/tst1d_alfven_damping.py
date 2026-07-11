#  1D Alfven-wave test with interior absorbing masks and in-domain field injection.
#
#  Setup:
#   - uniform guide field  B0  along x  (supports parallel-propagating Alfven waves)
#   - uniform magnetized plasma (ions + electrons, reduced mass ratio)
#   - two "damping" masks of thickness x_mask at both ends of the box
#     (absorbing boundary conditions for the waves)
#   - a FieldInjector driving Bz at a fixed x plane, just right of the left mask:
#     the left-going branch is absorbed by the left mask, the right-going wave
#     is used for the physics study.

import math

# --- guide field and geometry ---
B0      = 0.1          # uniform guide field along x
x_mask  = 20.          # mask thickness at each end
Lx      = 200.         # box length
dx      = 0.1
dt      = 0.95 * dx    # CFL
t_sim   = 400.

# --- Alfven driver ---
omega_inj = 0.05       # driver angular frequency
x_inj     = x_mask + 5.   # injection plane, just right of the left mask
w_inj     = 1.0        # driver spatial half-width
dB        = 1.e-3 * B0 # driver amplitude

Main(
    geometry = "1Dcartesian",
    interpolation_order = 2,

    cell_length = [dx],
    grid_length = [Lx],

    number_of_patches = [ 20 ],

    timestep = dt,
    simulation_time = t_sim,

    # absorbing masks at x-min and x-max:
    EM_boundary_conditions = [ ["damping", "damping"] ],
    EM_damping_thickness   = x_mask,
    EM_damping_coefficient = 1.0,

    print_every = 200,
)

# uniform guide field B0 along x
ExternalField(field="Bx", profile=B0)

# background magnetized plasma
for name, mass, charge in [("ions", 100., 1.), ("electrons", 1., -1.)]:
    Species(
        name = name,
        position_initialization = "regular",
        momentum_initialization = "maxwell-juettner",
        particles_per_cell = 64,
        mass = mass,
        charge = charge,
        number_density = 1.0,
        temperature = [1.e-4]*3,
        boundary_conditions = [ ["periodic"] ],
    )

# Alfven wave soft source: transverse B perturbation at a fixed x
FieldInjector(
    field = "Bz",
    profile = lambda x, t: dB * math.exp(-((x - x_inj)**2)/(2*w_inj**2)) * math.sin(omega_inj*t),
)

DiagScalar(every=200)

DiagFields(
    every = 200,
    fields = ["Ex","Ey","Ez","Bx","By","Bz"],
)

DiagProbe(
    every = 100,
    origin = [0.],
    corners = [[Lx]],
    number = [1000],
    fields = ["Ey","Ez","By","Bz"],
)
