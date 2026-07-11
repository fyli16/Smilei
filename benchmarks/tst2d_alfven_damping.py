#  2D Alfven-wave test with interior absorbing masks and in-domain field injection.
#
#  Setup:
#   - uniform guide field  B0  along x
#   - uniform magnetized plasma (ions + electrons, reduced mass ratio)
#   - two "damping" masks of thickness x_mask at both ends along x
#   - a FieldInjector driving Bz at a fixed x plane (just right of the left mask),
#     with a customizable transverse (y) profile.
#   - y is periodic.

import numpy as np

# --- guide field and geometry ---
B0       = 0.1
x_mask   = 20.
Lx, Ly   = 200., 50.
dx = dy  = 0.2
dt       = 0.95 / np.sqrt(1./dx**2 + 1./dy**2)   # CFL
t_sim    = 400.

# --- Alfven driver ---
omega_inj = 0.05
x_inj     = x_mask + 5.       # injection plane, just right of the left mask
w_inj     = 1.0              # driver spatial half-width along x
ky        = 2.*np.pi / Ly    # one transverse wavelength across the box
dB        = 1.e-3 * B0

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,

    cell_length = [dx, dy],
    grid_length = [Lx, Ly],

    number_of_patches = [ 20, 10 ],

    timestep = dt,
    simulation_time = t_sim,

    # masks along x, periodic along y:
    EM_boundary_conditions = [ ["damping", "damping"], ["periodic", "periodic"] ],
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
        particles_per_cell = 32,
        mass = mass,
        charge = charge,
        number_density = 1.0,
        temperature = [1.e-4]*3,
        boundary_conditions = [ ["periodic"], ["periodic"] ],
    )

# Alfven wave soft source with a custom transverse (y) profile at a fixed x plane
def inject_Bz(x, y, t):
    envelope_x = np.exp(-((x - x_inj)**2)/(2*w_inj**2))
    return dB * envelope_x * np.sin(omega_inj*t - ky*y)

FieldInjector(field="Bz", profile=inject_Bz)

DiagScalar(every=200)

DiagFields(
    every = 200,
    fields = ["Ex","Ey","Ez","Bx","By","Bz"],
)

DiagProbe(
    every = 200,
    origin  = [0.,  Ly/2.],
    corners = [[Lx, Ly/2.]],
    number  = [1000],
    fields = ["Ey","Ez","By","Bz"],
)
