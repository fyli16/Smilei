#  2D kinetic/inertial Alfven-wave test with DIVERGENCE-SAFE electric-field injection.
#
#  Why inject E (not B) at finite perpendicular scale:
#   The Yee solver only changes B through Faraday (dB/dt = -curl E), which is
#   divergence-free by construction, and there is NO divergence cleaning. A direct
#   B soft-source with finite transverse structure (d(dBy)/dy != 0) would seed a
#   growing div(B) != 0 error. Injecting E instead lets the code generate a
#   divergence-free B automatically -- including the compressional Bx that a
#   finite-k_perp wave must carry.
#
#  Geometry: guide field B0 along x, propagation mainly along x, finite k_perp along y.
#   - "shear" KAW  : drive Ey (in-plane transverse E) -> out-of-plane Bz, div(B)=0, Bx-free
#   - "circular"   : also drive Ez (out-of-plane E) 90 deg out of phase -> in-plane By and
#                    the accompanying Bx together, divergence-free.
#
#  y is periodic (sets a clean k_perp); x uses the interior "damping" masks.

import numpy as np
twopi = 2*np.pi

# ---------------- plasma parameters (same normalization as inp.py) ----------------
mi      = 100.0        # ion mass (unit me)
wpiwci  = 20.
Ti_Te   = 1./4
beta    = 5e-4
beta_e  = beta / (1+Ti_Te)
beta_i  = beta / (1+1/Ti_Te)
Te1     = 0.25*beta_e*mi/wpiwci**2
Ti1     = 0.25*beta_i*mi/wpiwci**2
Te      = [Te1]*3
Ti      = [Ti1]*3
n0      = mi
ppc     = 100
B0      = mi/wpiwci

# ---------------- geometry ----------------
box  = [128., 8.]           # di ; Ly finite so we can resolve k_perp
res  = [32, 32]             # cells per di
number_of_patches = [128, 8]

wci_tmax = 250.
tmax = wci_tmax * wpiwci

dr = [1./i for i in res]
dt_CFL = 1./np.sqrt(sum(i**2 for i in res))
rest = (1./dt_CFL)/0.98
dt = 1./rest
Nt = int(rest)
dim = len(box)
for i in range(dim):
    fac = box[i]*res[i]/number_of_patches[i]
    if fac != fac//1: fac = int(fac)+1
    box[i] = float(fac)*number_of_patches[i]/res[i]

dump_freq  = Nt*wpiwci
flush_freq = 10*dump_freq

Main(
    geometry = "2Dcartesian",
    interpolation_order = 2,
    grid_length = box,
    cell_length = dr,
    number_of_patches = number_of_patches,
    timestep = dt,
    simulation_time = tmax,
    maxwell_solver = 'Yee',
    EM_boundary_conditions = [ ['damping'], ['periodic'] ],
    EM_damping_thickness   = 0.2*box[0],
    EM_damping_coefficient = 0.1,
    solve_poisson = False,
    print_every = int(Nt*wpiwci),
)

CurrentFilter(model="binomial", passes=[10], kernelFIR=[0.25, 0.5, 0.25])

# uniform guide field along x
ExternalField(field="Bx", profile=B0)

Species(name='ion', position_initialization='random', momentum_initialization='maxwell-juettner',
        temperature=Ti, particles_per_cell=ppc, mass=mi, charge=1.0, number_density=n0,
        boundary_conditions=[['thermalize'],['periodic']], thermal_boundary_temperature=Ti, time_frozen=0.0)

Species(name='elec', position_initialization='random', momentum_initialization='maxwell-juettner',
        temperature=Te, particles_per_cell=ppc, mass=1.0, charge=-1.0, number_density=n0,
        boundary_conditions=[['thermalize'],['periodic']], thermal_boundary_temperature=Te, time_frozen=0.0)

# ================= Alfven-wave injection via the ELECTRIC field =================
dE       = 0.02 * B0        # driver amplitude (E in units where c=1 -> comparable to dB)
w0wci    = 0.1*twopi        # driver frequency in units of wci
w0wpi    = w0wci/wpiwci
t_upramp = 50.*wpiwci
x_inj    = 30.0             # injection plane (just right of the left mask at 0.2*128=25.6 di)

# perpendicular wavenumber: m_perp transverse wavelengths across Ly (finite k_perp -> KAW)
m_perp = 1
ky     = twopi * m_perp / box[1]

def ramp(t):
    return np.sin(0.5*np.pi*t/t_upramp)**2 if t < t_upramp else 1.

# single x-cell source (Ey is x-primal -> node at integer*dx)
def _column(x):
    return np.abs(x - x_inj) < 0.5*dr[0]

# ---- clean SHEAR KAW: drive perpendicular in-plane Ey -> out-of-plane Bz (no Bx) ----
def inject_Ey(x, y, t):
    return dE * ramp(t) * _column(x) * np.cos(ky*y) * np.sin(w0wpi*t)

FieldInjector(field="Ey", profile=inject_Ey)

# ---- CIRCULAR variant: uncomment to also drive out-of-plane Ez 90 deg out of phase.
#      Ez generates By AND the divergence-safe Bx together; with Ey above this gives a
#      rotating (By, Bz) circular wave. Flip the sign of the cos term to swap handedness.
# def inject_Ez(x, y, t):
#     return dE * ramp(t) * _column(x) * np.cos(ky*y) * np.cos(w0wpi*t)
# FieldInjector(field="Ez", profile=inject_Ez)

# ---------------- diagnostics ----------------
DiagScalar(every=int(Nt), precision=8)

DiagFields(
    every = [0, dump_freq],
    flush_every = flush_freq,
    fields = ["Ex","Ey","Ez","Bx","By","Bz",
              "Rho_elec","Rho_ion","Jx","Jy","Jz"],
)

# probe along x at mid-y to track propagation / polarization
from numpy import s_
DiagProbe(
    every = int(Nt),
    origin  = [0.,      box[1]/2.],
    corners = [[box[0], box[1]/2.]],
    number  = [int(box[0]*res[0])],
    fields = ["Ex","Ey","Ez","Bx","By","Bz"],
)
