"""
This is a typical input script that runs a simulation of
laser-wakefield acceleration using FBPIC.

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
from fbpic.lpa_utils.laser import add_laser_pulse, GaussianLaser
from fbpic.lpa_utils.boosted_frame import BoostConverter
from fbpic.openpmd_diag import FieldDiagnostic, \
                  BoostedFieldDiagnostic, BoostedParticleDiagnostic
# ----------
# Parameters
# ----------
use_cuda = False

# The simulation box
Nz = 400         # Number of gridpoints along z
zmax = 0.e-6     # Length of the box along z (meters)
zmin = -20.e-6
Nr = 50          # Number of gridpoints along r
rmax = 75.e-6   # Length of the box along r (meters)
Nm = 2           # Number of modes used
# Boosted frame
gamma_boost = 10.
# The simulation timestep
dt = min( rmax/(2*gamma_boost*Nr), (zmax-zmin)/Nz/c )  # Timestep (seconds)
# (See the section Advanced use > Running boosted-frame simulation
# of the FBPIC documentation for an explanation of the above calculation of dt)
N_step = 301     # Number of iterations to perform

# Boosted frame converter
boost = BoostConverter(gamma_boost)

# The laser (conversion to boosted frame is done inside 'add_laser_pulse')
a0 = 4.          # Laser amplitude
w0 = 25.e-6      # Laser waist
tau = 15.e-15    # Laser duration
z0 = -10.e-6     # Laser centroid
zfoc = 0.e-6     # Focal position
lambda0 = 0.8e-6 # Laser wavelength

# The density profile
ramp_up = 500.e-6

# The particles of the plasma
p_zmin = 0.e-6   # Position of the beginning of the plasma (meters)
p_zmax = 5.e-3
p_rmin = 0.      # Minimal radial position of the plasma (meters)
p_rmax = 130.e-6 # Maximal radial position of the plasma (meters)
n_e = 2.e18*1.e6 # The density in the lab frame (electrons.meters^-3)
p_nz = 2         # Number of particles per cell along z
p_nr = 2         # Number of particles per cell along r
p_nt = 6         # Number of particles per cell along theta

# Density profile
# Convert parameters to boosted frame
# (NB: the density is converted inside the Simulation object)
ramp_up,  = boost.static_length( [ ramp_up ] )
# Define the density function
def dens_func( z, r ):
    """
    Parabolic density profile
    """
    # Allocate relative density
    n = np.ones_like(z)
    # Make ramp up
    n = np.where( z<ramp_up, z/ramp_up, n )
    n = np.where( z<0, 0, n )
    return(n)

# The diagnostics in the boosted frame
diag_period = 50
# The diagnostics in the lab frame
Ntot_snapshot_lab = 20
dt_snapshot_lab = 6*(zmax-zmin)/c

# ---------------------------
# Carrying out the simulation
# ---------------------------

# Velocity of the Galilean frame (for suppression of the NCI)
v_comoving = - c * np.sqrt( 1. - 1./gamma_boost**2 )
# Initialize the simulation object
sim = Simulation( Nz, zmax, Nr, rmax, Nm, dt,
    p_zmin, p_zmax, p_rmin, p_rmax, p_nz, p_nr, p_nt, n_e,
    dens_func=dens_func, zmin=zmin, initialize_ions=True,
    v_comoving=v_comoving, gamma_boost=gamma_boost,
    boundaries='open', use_cuda=use_cuda )

# Add a laser to the fields of the simulation
laser_profile = GaussianLaser( a0, w0, tau, z0, zf=zfoc )
add_laser_pulse( sim, laser_profile, gamma_boost=gamma_boost )

# Configure the moving window
sim.set_moving_window( v=c )

# Add a field diagnostic
sim.diags = [
    FieldDiagnostic(diag_period, sim.fld, sim.comm ),
    BoostedFieldDiagnostic( zmin, zmax, c,
                dt_snapshot_lab, Ntot_snapshot_lab, gamma_boost,
                period=diag_period, fldobject=sim.fld, comm=sim.comm)
    ]

### Run the simulation
sim.step( N_step )
