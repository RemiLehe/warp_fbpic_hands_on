"""
This is a typical input script that runs a simulation of
laser-wakefield acceleration in a boosted-frame, using Warp in 2D / Circ / 3D.

Usage
-----
- Modify the parameters below to suit your needs
- Type "python -i lpa_boostedframe_script.py" in a terminal
- When the simulation finishes, the python session will *not* quit.
    Therefore the simulation can be continued by running step()
    Otherwise, one can just type exit()
"""
# Import warp-specific packages
from warp.init_tools import *

# -----------------------------------------------------------------------------
# Parameters (Modify the values below to suit your needs)
# -----------------------------------------------------------------------------

# General parameters
# ------------------
# Dimension of simulation ("3d", "circ" or "2d", "1d")
dim = "2d"
# Number of azimuthal modes beyond m=0, for "circ" (not used for "2d" and "3d")
circ_m = 1
# Total number of timesteps in the simulation
N_steps = 301

# Simulation box
# --------------
# Number of grid cells in the longitudinal direction
Nz = 400
# Number of grid cells in transverse direction (represents Nr in "circ")
Nx = 100
# Number of grid cells in the 3rd dimension (not used for "2d" and "circ")
Ny = 100
# Dimension of the box in longitudinal direction (meters)
zmin_lab = -20.e-6
zmax_lab = 0.e-6
# Dimension of the box in transverse direction (box ranges from -xmax to xmax)
xmax = 75.e-6
# Dimension of the box in 3rd direction (not used for "2d" and "circ")
ymax = 75.e-6

# Field boundary conditions (longitudinal and transverse respectively)
f_boundz  = openbc
f_boundxy = openbc
if dim=="circ":
    f_boundxy = dirichlet
# Particles boundary conditions (longitudinal and transverse respectively)
p_boundz  = absorb
p_boundxy = absorb

# Boosted frame
gamma_boost = 10.

# Diagnostics
# -----------
# Period of diagnostics (in number of timesteps)
diag_period = 50
Ntot_snapshot_lab = 20
dt_snapshot_lab = 6*(zmax_lab-zmin_lab)/clight

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
# Laser amplitude at focus
laser_a0 = 4.
# Waist at focus (meters)
laser_w0 = 25.e-6
# Length of the pulse (length from the peak to 1/e of the amplitude ; meters)
laser_ctau = 5.e-6
# Initial position of the centroid (meters)
laser_z0 = -2 * laser_ctau
# Focal position
laser_zfoc = 0.e-6
# Position of the antenna (meters)
laser_source_z = -0.1e-6
# Polarization angle with respect to the x axis (rad)
laser_polangle = 0.
# Wavelength
laser_lambda0 = 0.8e-6

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
n_plasma =  2.e18*1.e6
# The different elements used. (Only used if use_ions is different than 0.)
# relative_density is the density relative to n_plasma.
# q_start is the ionization state of the ions at the beginning of the simulation
# q_max is the maximum ionization state
# If q_start is not equal to q_max, ionization between states will be computed.
ion_states = { 'Hydrogen': {'relative_density':1., 'q_start':1, 'q_max':1} }
# Positions between which the plasma is initialized
# (Transversally, the plasma is initialized between -plasma_xmax and
# plasma_xmax, along x, and -plasma_ymax and plasma_ymax along y)
plasma_zmin = 0.e-6
plasma_zmax = 5.e-3
plasma_xmax = 70.e-6
plasma_ymax = 70.e-6
plasma_ramp = 500.e-6

# Perform a boost of the different quantities
# -------------------------------------------
boost = BoostConverter( gamma_boost )
# Plasma
n_plasma, = boost.static_density([ n_plasma ])
plasma_uz_m, = boost.longitudinal_momentum([ 0. ])
plasma_ramp, = boost.static_length([ plasma_ramp ])
def plasma_dens_func( x, y, z ):
    # Allocate relative density
    n = np.ones_like(z)
    # Make ramp up
    n = np.where( z<ramp_up, z/ramp_up, n )
    n = np.where( z<0, 0, n )
    return(n)
# Simulation box
zmin, zmax = boost.copropag_length([ zmin_lab, zmax_lab ])
# NB: Do not boost the laser quantities: these are boosted in add_laser

# -----------------------------------------------------------------------------
# Initialization of the simulation (Normal users should not modify this part.)
# -----------------------------------------------------------------------------

# Set some general options for warp
set_diagnostics( 0 )
set_boundary_conditions( f_boundz, f_boundxy, p_boundz, p_boundxy )
set_simulation_box( Nz, Nx, Ny, zmin, zmax, xmax, ymax, dim )
set_moving_window( l_moving_window=1, v_moving_window=clight )

# See smoothing.py
set_smoothing_parameters( use_smooth, dim, npass_smooth,
                         alpha_smooth, stride_smooth )

# Creation of the species
# -----------------------

elec = None
ions = None
elec_from_ions = None
# Create the plasma species
# Reference weight for plasma species
plasma_weight = prepare_weights( n_plasma, plasma_nx, plasma_ny,
                            plasma_nz, dim, circ_m )
elec = Species(type=Electron, weight=plasma_weight, name='electrons')
ions, elec_from_ions = initialize_ion_dict( ion_states, plasma_weight )

# Set the numerical parameters only now: they affect the newly created species
set_numerics( depos_order, efetch, particle_pusher, dim)

# Setup the field solver object
# -----------------------------
em = EM3D(
    stencil = stencil,
    npass_smooth = npass_smooth,
    alpha_smooth = alpha_smooth,
    stride_smooth = stride_smooth,
    l_2dxz = (dim in ["2d", "circ"]),
    l_2drz = (dim in ["circ"]),
    l_1dz = (dim=="1d"),
    l_getrho = True,
    circ_m = (dim == "circ")*circ_m,
    l_correct_num_Cherenkov = True,
    type_rz_depose = 1,
    l_setcowancoefs = True )
registersolver(em)

# Introduce the laser
# -------------------
add_laser( em, dim, laser_a0, laser_w0, laser_ctau, laser_z0,
    zf=laser_zfoc, theta_pol=laser_polangle, source_z=laser_source_z,
    lambda0=laser_lambda0, gamma_boost=gamma_boost )

# Introduce the plasma
# --------------------
# Create an object to store the information about plasma injection
plasma_injector = PlasmaInjector( elec, ions, w3d, top, dim,
        plasma_nx, plasma_ny, plasma_nz, plasma_zmin,
        plasma_zmax, plasma_xmax, plasma_ymax, None,
        uz_m = plasma_uz_m )
# Continuously inject the plasma, if the moving window is on
installuserinjection( plasma_injector.continuous_injection )

# Setup the diagnostics
# ---------------------
remove_existing_directory( ['diags', 'lab_diags'] )
diag1 = FieldDiagnostic( period=diag_period, top=top, w3d=w3d, em=em,
            comm_world=comm_world )
installafterstep( diag1.write )
diag0 = BoostedFieldDiagnostic( zmin_lab, zmax_lab, clight,
    dt_snapshot_lab, Ntot_snapshot_lab, gamma_boost,
    period=diag_period, top=top, w3d=w3d, em=em, comm_world=comm_world )
installafterstep( diag0.write )

# ---------------
# Simulation loop
# ---------------
# Split loop for printing purposes
for n in range( int(N_steps/10) ):
    step(10)
