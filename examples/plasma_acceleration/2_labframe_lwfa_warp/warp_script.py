"""
This is a typical input script that runs a simulation of
laser-wakefield acceleration using Warp in 2D / Circ / 3D.

Usage
-----
- Modify the parameters below to suit your needs
- Type "python warp_script.py" in a terminal
"""
# Import warp-specific packages
from warp.init_tools import *

# -----------------------------------------------------------------------------
# Parameters (Modify the values below to suit your needs)
# -----------------------------------------------------------------------------

# General parameters
# ------------------
# Dimension of simulation ("3d", "circ", "2d")
dim = "2d"
# Number of azimuthal modes beyond m=0, for "circ" (not used for "2d" and "3d")
circ_m = 1
# Total number of timesteps in the simulation
N_steps = 500

# Simulation box
# --------------
# Number of grid cells in the longitudinal direction
Nz = 200
# Number of grid cells in transverse direction (represents Nr in "circ")
Nx = 100
# Number of grid cells in the 3rd dimension (not used for "2d" and "circ")
Ny = 100
# Dimension of the box in longitudinal direction (meters)
zmin = -15.e-6
zmax = 5.e-6
# Dimension of the box in transverse direction (box ranges from -xmax to xmax)
xmax = 15.e-6
# Dimension of the box in 3rd direction (not used for "2d" and "circ")
ymax = 15.e-6

# Field boundary conditions (longitudinal and transverse respectively)
f_boundz  = openbc
f_boundxy = openbc
if dim == "circ":
    f_boundxy = dirichlet
# Particles boundary conditions (longitudinal and transverse respectively)
p_boundz  = absorb
p_boundxy = absorb

# Diagnostics
# -----------
# Period of diagnostics (in number of timesteps)
diag_period = 20

# Numerical parameters
# --------------------
# Field solver (0:Yee, 1:Karkkainen on EF,B, 3:Lehe)
stencil = 1
# Particle shape (1:linear, 2:quadratic, 3:cubic)
depos_order = 1
# Gathering mode (1:from cell centers, 4:from Yee mesh)
efetch = 1
# Particle pusher (0:Boris, 1:Vay)
particle_pusher = 1

# Current smoothing parameters
# ----------------------------
# Turn current smoothing on or off (0:off; 1:on)
use_smooth = 1
# Number of passes of smoother and compensator in each direction (x, y, z)
npass_smooth = array([[ 1 , 0 ], [ 1 , 0 ], [ 1 , 0 ]])
# Smoothing coefficients in each direction (x, y, z)
alpha_smooth = array([[ 0.5, 3./2], [ 0.5, 3./2], [0.5, 3./2]])
# Stride in each direction (x, y, z)
stride_smooth = array([[ 1 , 1 ], [ 1 , 1 ], [ 1 , 1 ]])

# Laser parameters
# ----------------
# Position of the antenna (meters)
laser_source_z = 0.e-6
# Polarization angle with respect to the x axis (rad)
laser_polangle = 0.
# Gaussian pulse:
# Laser amplitude at focus
laser_a0 = 1.
# Waist at focus (meters)
laser_w0 = 4.e-6
# Length of the pulse (length from the peak to 1/e of the amplitude ; meters)
laser_ctau = 3.e-6
# Initial position of the centroid (meters)
laser_z0 = -1.5 * laser_ctau
# Focal position
laser_zfoc = 4.5e-05

# Plasma macroparticles
# ---------------------
# Number of macroparticles per cell in each direction
# In Circ, nppcelly is the number of particles along the
# azimuthal direction. Use a multiple of 4*circ_m
plasma_nx = 2
plasma_ny = 4
plasma_nz = 2

# Plasma content and profile
# --------------------------
# Reference plasma density (in number of particles per m^3)
n_plasma = 2.5e25
# Positions between which the plasma is initialized
# (Transversally, the plasma is initialized between -plasma_xmax and
# plasma_xmax, along x, and -plasma_ymax and plasma_ymax along y)
plasma_zmin = 0.e-6
plasma_zmax = 1500.e-6
plasma_xmax = xmax
plasma_ymax = ymax

# Define your own profile and profile parameters below
ramp_length = 20.e-6
def plasma_dens_func( x, y, z ):
    """
    User-defined function: density profile of the plasma

    It should return the relative density with respect to n_plasma,
    at the position x, y, z (i.e. return a number between 0 and 1)

    Parameters
    ----------
    x, y, z: 1darrays of floats
        Arrays with one element per macroparticle
    Returns
    -------
    n : 1d array of floats
        Array of relative density, with one element per macroparticles
    """
    # Allocate relative density
    n = ones_like(z)
    # Make linear ramp
    n = where( z<ramp_length, z/ramp_length, n )
    # Supress density before the ramp
    n = where( z<0, 0., n )

    return(n)

# -----------------------------------------------------------------------------
# Initialization of the simulation (Normal users should not modify this part.)
# -----------------------------------------------------------------------------

# Set some general options for warp
set_diagnostics( interactive=0 )
set_boundary_conditions( f_boundz, f_boundxy, p_boundz, p_boundxy )
set_simulation_box( Nz, Nx, Ny, zmin, zmax, xmax, ymax, dim )
set_moving_window( l_moving_window=1, v_moving_window=clight )

# See smoothing.py
set_smoothing_parameters( use_smooth, dim, npass_smooth,
                         alpha_smooth, stride_smooth )

# Creation of the species
# -----------------------
# Create the plasma species
# Reference weight for plasma species
elec_weight = prepare_weights( n_plasma, plasma_nx, plasma_ny,
                            plasma_nz, dim, circ_m )
elec = Species(type=Electron, weight=elec_weight, name='electrons')

# Set the numerical parameters only now: they affect the newly created species
set_numerics( depos_order, efetch, particle_pusher, dim)

# Setup the field solver object
# -----------------------------
em = initialize_em_solver( stencil, dim,
    npass_smooth, alpha_smooth, stride_smooth,
    circ_m = (dim =="circ")*circ_m )
registersolver(em)

# Introduce the laser
# -------------------
add_laser( em, dim, laser_a0, laser_w0, laser_ctau, laser_z0,
    zf=laser_zfoc, theta_pol=laser_polangle, source_z=laser_source_z )

# Introduce the plasma
# --------------------
# Create an object to store the information about plasma injection
plasma_injector = PlasmaInjector( elec, None, w3d, top, dim,
        plasma_nx, plasma_ny, plasma_nz, plasma_zmin,
        plasma_zmax, plasma_xmax, plasma_ymax, plasma_dens_func )
# Continuously inject the plasma, if the moving window is on
installuserinjection( plasma_injector.continuous_injection )

# Setup the diagnostics
# ---------------------
remove_existing_directory( ['diags'] )
# Particle output
diag1 = FieldDiagnostic( period=diag_period, top=top, w3d=w3d, em=em,
            comm_world=comm_world, lparallel_output=False )
installafterstep( diag1.write )
# Field output
diag2 = ParticleDiagnostic( period=diag_period, top=top, w3d=w3d,
        species={'plasma electrons': elec},
        comm_world=comm_world, lparallel_output=False )
installafterstep( diag2.write )

print('\nInitialization complete\n')

# ---------------
# Simulation loop
# ---------------
# Split loop for printing purposes
for n in range( int(N_steps/10) ):
    step(10)
