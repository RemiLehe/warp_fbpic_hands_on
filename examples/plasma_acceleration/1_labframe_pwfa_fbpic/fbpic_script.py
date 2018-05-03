"""
This is a typical input script that runs a simulation of
beam-driven plasma wakefield acceleration using FBPIC, in the lab frame.

Usage
-----
- Modify the parameters below to suit your needs
- Type "python fbpic_script.py" in a terminal
"""

# -------
# Imports
# -------
import numpy as np
from scipy.constants import c
# Import the relevant structures in FBPIC
from fbpic.main import Simulation
from fbpic.openpmd_diag import FieldDiagnostic, ParticleDiagnostic
from fbpic.lpa_utils.bunch import add_elec_bunch

# ----------
# Parameters
# ----------

# Whether to use the GPU
use_cuda = False

# The simulation box
Nz = 200         # Number of gridpoints along z
zmax = 30.e-6    # Right end of the simulation box (meters)
zmin = -10.e-6   # Left end of the simulation box (meters)
Nr = 100         # Number of gridpoints along r
rmax = 20.e-6    # Length of the box along r (meters)
Nm = 1           # Number of modes used

# The simulation timestep
dt = (zmax-zmin)/Nz/c   # Timestep (seconds)
N_step = 400     # Number of iterations to perform

# The plasma particles
plasma_zmin = 0.      # Position of the beginning of the plasma (meters)
plasma_zmax = 1.e-2   # Position of the end of the plasma (meters)
plasma_rmin = 0.      # Minimal radial position of the plasma (meters)
plasma_rmax = 20.e-6  # Maximal radial position of the plasma (meters)
plasma_n = 1.e18*1.e6 # Density (electrons.meters^-3)
plasma_nz = 2         # Number of particles per cell along z
plasma_nr = 2         # Number of particles per cell along r
plasma_nt = 4         # Number of particles per cell along theta
# The plasma density profile
ramp_length = 50.e-6
def plasma_dens_func( z, r ) :
    """Linear ramp followed by flat profile"""
    # Allocate relative density
    n = np.ones_like(z)
    # Make linear ramp
    n = np.where( z<ramp_length, z/ramp_length, n )
    return(n)

# The driver beam (flat-top)
beam_gamma = 2.e3     # Corresponds to 1 GeV
beam_n = 5.e18*1.e6  # Beam density
beam_zmin = 15.e-6
beam_zmax = 20.e-6
beam_rmin = 0.e-6
beam_rmax = 2.e-6
beam_nz = 2           # Number of particles per cell along z
beam_nr = 2           # Number of particles per cell along r
beam_nt = 4           # Number of particles per cell along theta

# The diagnostics and the checkpoints/restarts
diag_period = 20    # Period of the diagnostics in number of timesteps


# ---------------------------
# Carrying out the simulation
# ---------------------------

# Initialize the simulation object
sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt,
    plasma_zmin, plasma_zmax, plasma_rmin, plasma_rmax,
    plasma_nz, plasma_nr, plasma_nt, plasma_n,
    dens_func=plasma_dens_func,
    zmin=zmin, boundaries='open', use_cuda=use_cuda )
# At this point, the code initialized one electron species for the plasma

# Add the driver electron bunch
add_elec_bunch( sim, beam_gamma, beam_n,
    beam_zmin, beam_zmax, beam_rmin, beam_rmax,
    beam_nz, beam_nr, beam_nt )

# Configure the moving window
sim.set_moving_window( v=c )

# Add diagnostics
sim.diags = [
    FieldDiagnostic( diag_period, sim.fld, comm=sim.comm ),
    ParticleDiagnostic( diag_period,
        {"plasma electrons" : sim.ptcl[0],
         "driver beam" : sim.ptcl[1]},
        comm=sim.comm ) ]

### Run the simulation
sim.step( N_step )
