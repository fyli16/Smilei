import numpy as np
twopi = 2*np.pi
twopi2 = twopi**2.

mi = 100.0  # ion mass, unit me
wpiwci = 20

Ti_Te = 1./4  # ratio Ti/Te
beta = 5e-4 # total beta
beta_e = beta / (1+Ti_Te)
beta_i = beta / (1+1/Ti_Te)
Te1 = 0.25*beta_e*mi/wpiwci**2.  # Te (one component)
Ti1 = 0.25*beta_i*mi/wpiwci**2.  # Ti (one component)
Te = [Te1, Te1, Te1]
Ti = [Ti1, Ti1, Ti1]

n0 = mi # plasma density in units of nr
ppc = 100 # particles per cell, same for e, i


# ---------- basic setup -------------------
box = [128, 1]  # di
res = [32, 32]   # resolution: cell number per unit
number_of_patches = [256, 4]

wci_tmax = 250. # max run time in units of wci^-1
tmax = wci_tmax * wpiwci  # max run time, wpi^-1

#********** DO NOT CHANGE **********************
dr =[ 1./i for i in res ]
dt_CFL = 1./np.sqrt(sum(i**2 for i in res))  # Courant condition
rest = (1./dt_CFL)/0.98
dt =1./rest
# timesteps to finish one time unit
Nt = int(rest)
# sim. dimension
dim = len(box)
# adjust box size (cells) to be divided by number_of_pathces
for i in range(dim):
	fac = float(box[i])*res[i]/number_of_patches[i]
	if fac!=fac//1: fac = int(fac)+1
	box[i] = float(fac)*number_of_patches[i]/res[i]

# total number of cells
cells = list( np.array(box)*np.array(res) )
#************************************************

dump_freq  = Nt*wpiwci  # dump every wci^-1  
flush_freq = 10*dump_freq


# ------------ global settings --------------------
Main(
  geometry = "2Dcartesian",
  interpolation_order = 2,
  grid_length  = box,
  cell_length = dr,
  # each must be power of 2; total must be greater than MPI process number
  number_of_patches = number_of_patches,
  #patch_arrangement = "linearized_ZYX",
  timestep = dt,
  simulation_time = tmax,
  maxwell_solver = 'Yee',
  EM_boundary_conditions = [
    ['damping'], ['periodic'],
    # ['silver-muller'], ['periodic'],
    # ['silver-muller'],
  ],
  EM_damping_thickness = 0.2*box[0],
  EM_damping_coefficient = 0.1,
#   number_of_pml_cells=[[64, 64]],
  solve_poisson = False,  # "False" to have a neutralizing background
  # reference_angular_frequency_SI = 3.e8/0.4e-6, #speed of light / laser wavelegth
  print_every = int(Nt*wpiwci), # print every wci^-1
  # random_seed = smilei_mpi_rank
)

CurrentFilter(
     model = "binomial",
     passes = [10],
     #kernelFIR = [0.1, 0.2, 0.3, 0.2, 0.1],
     #kernelFIR = [0.05, 0.1, 0.2, 0.3, 0.2, 0.1, 0.05],
     kernelFIR = [0.25, 0.5, 0.25],
)

B0 = mi/wpiwci

ExternalField(
    field = "Bx",
    profile = B0,
)


# ion
Species(name = 'ion',position_initialization = 'random',
        momentum_initialization = 'maxwell-juettner', 
        temperature = Ti, particles_per_cell = ppc, 
        mass = mi, charge = 1.0, number_density = n0, 
        boundary_conditions = [['thermalize'],['periodic']], 
        thermal_boundary_temperature = Ti, time_frozen=0.0,
)

# electron
Species(name = 'elec',position_initialization = 'random', 
        momentum_initialization = 'maxwell-juettner', 
        temperature = Te, particles_per_cell = ppc, 
        mass = 1.0, charge = -1.0, number_density = n0, 
        boundary_conditions = [['thermalize'],['periodic']], 
        thermal_boundary_temperature = Te, time_frozen=0.0,
)

# Alfvén wave soft source: transverse B perturbation at a fixed x

dB_B0 = 0.02*B0
w0wci = 0.1*twopi   # w0/wci = 0.6283
w0wpi = w0wci/wpiwci
t_upramp = 50.*wpiwci  # wci^-1
theta = 0
i_inj = int(30.0 / dr[0])  # injection cell index
x_inj = (i_inj+0.5) * dr[0]  # injection cell location

def inject_Bz(x, y, t):
    rt = np.sin(0.5 * np.pi * t / t_upramp)**2 if t < t_upramp else 1. # ramp in time
    column = np.abs(x-x_inj) < 0.5 * dr[0]
    return dB_B0 * rt * column * np.sin(w0wpi*t+theta)

FieldInjector(
    field = "Bz",
    profile = inject_Bz
)



#--------------------------------------------------
# ---------------- diagnostics --------------------
#--------------------------------------------------
DiagScalar(
    every = int(Nt),
    precision = 8,
)

#---------Field --------------
from numpy import s_
# central x-y plane
DiagFields(
  every = [0, dump_freq],
  flush_every = flush_freq,
#   subgrid = s_[:, :, int(cells[2]/2.)],
  fields = ["Ex", 'Ey', 'Ez', 'Bx', 'By', 'Bz', 
            'Rho_elec', 'Rho_ion',
            'Jx', 'Jy', 'Jz',
            'Jx_elec', 'Jy_elec', 'Jz_elec',
            'Jx_ion', 'Jy_ion', 'Jz_ion',
            ],
)
