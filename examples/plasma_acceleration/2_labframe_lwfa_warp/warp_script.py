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
stencil = 0
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
npass_smooth = array([[ 0 , 0 ], [ 0 , 0 ], [ 1 , 1 ]])
# Smoothing coefficients in each direction (x, y, z)
alpha_smooth = array([[ 0.5, 3.], [ 0.5, 3.], [0.5, 3./2]])
# Stride in each direction (x, y, z)
stride_smooth = array([[ 1 , 1 ], [ 1 , 1 ], [ 1 , 1 ]])

# Laser parameters
# ----------------
# Initialize laser (0:off, 1:on)
use_laser = 1
# Position of the antenna (meters)
laser_source_z = 0.e-6
# Polarization angle with respect to the x axis (rad)
laser_polangle = pi/2

# Laser file:
# When using a laser profile that was experimentally
# measured, provide a string with the path to an HDF5 laser file,
# otherwise provide None and a Gaussian pulse will be initialized
laser_file = None
laser_file_energy = 2. # When using a laser file, energy in Joule of the pulse

# Gaussian pulse:
# Laser amplitude at focus
laser_a0 = 1.
# Waist at focus (meters)
laser_w0 = 4.e-6
# Length of the pulse (length from the peak to 1/e of the amplitude ; meters)
laser_ctau = 3.e-6
# Initial position of the centroid (meters)
laser_z0 = -2 * laser_ctau
# Focal position
laser_zfoc = 4.5e-05

# Plasma macroparticles
# ---------------------
# Initialize some preexisting plasmas electrons (0:off, 1:on)
# (Can be used in order to neutralize pre-ionized ions, if any,
# or in order to simulate a plasma without having to initialize ions)
use_preexisting_electrons = 1
# Initialize plasma ions (0:off, 1:on)
use_ions = 1
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
# Relative density of the preexisting electrons (relative to n_plasma)
rel_dens_preexisting_electrons = 1.
# The different elements used. (Only used if use_ions is different than 0.)
# relative_density is the density relative to n_plasma.
# q_start is the ionization state of the ions at the beginning of the simulation
# q_max is the maximum ionization state
# If q_start is not equal to q_max, ionization between states will be computed.
ion_states = { 'Hydrogen': {'relative_density':1., 'q_start':1, 'q_max':1 },
                 'Helium': {'relative_density':0.25, 'q_start':0, 'q_max':2 } }
# Positions between which the plasma is initialized
# (Transversally, the plasma is initialized between -plasma_xmax and
# plasma_xmax, along x, and -plasma_ymax and plasma_ymax along y)
plasma_zmin = 1.e-6
plasma_zmax = 1500.e-6
plasma_xmax = xmax
plasma_ymax = ymax

# Define your own profile and profile parameters below
ramp_start = 0.e-6
ramp_length = 20.e-6
ramp_plateau = 20.e-6
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
    n = where( z<ramp_start+ramp_length, (z-ramp_start)/ramp_length, n )
    # Supress density before the ramp
    n = where( z<ramp_start, 0., n )
    # Reduce density by half after the ramp
    n = where( z> ramp_start+ramp_length+ramp_plateau, 0.5*n, n )
    # Put the density to 0 later
    n = where( z> ramp_start+ramp_length+2*ramp_plateau, 0., n )

    return(n)

# Relativistic beam
# -----------------
# Initialize beam electrons (0:off, 1:on)
# (Please be aware that initializing a beam in 2D geometry makes very little
# physical sense, because of the long range of its space-charge fields)
use_beam = 0
# Longitudinal momentum of the beam
beam_uz = 100.
# Beam density
n_beam = 1.e26
# Number of macroparticles per cell in each direction
beam_nx = 2*plasma_nx
beam_ny = 2*plasma_ny
beam_nz = 2*plasma_nz
# Positions between which the beam is initialized
# (Transversally, the plasma is initialized between -plasma_xmax and
# plasma_xmax, along x, and -plasma_ymax and plasma_ymax along y)
beam_zmin = -12.e-6
beam_zmax = -10.e-6
beam_xmax = 3.e-6
beam_ymax = 3.e-6

# Define your own profile and profile parameters below
beam_rmax = beam_xmax
def beam_dens_func(x, y, z):
    """
    User-defined function: density profile of the beam

    It should return the relative density with respect to n_beam,
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
    # Longitudinal profile: parabolic
    n = n*(z - beam_zmin)*(beam_zmax - z) * 4/(beam_zmax - beam_zmin)**2
    # Transverse profile: parabolic
    r = sqrt( x**2 + y**2)
    n = n*(1 - (r/beam_rmax)**2 )
    # Put the density above rmax to 0
    n[r > beam_rmax] = 0.

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

elec = None
ions = None
elec_from_ions = None
beam = None
# Create the plasma species
# Reference weight for plasma species
plasma_weight = prepare_weights( n_plasma, plasma_nx, plasma_ny,
                            plasma_nz, dim, circ_m )
if use_preexisting_electrons:
    elec_weight = rel_dens_preexisting_electrons * plasma_weight
    elec = Species(type=Electron, weight=elec_weight, name='electrons')
if use_ions:
    ions, elec_from_ions = initialize_ion_dict( ion_states, plasma_weight,
                                                group_elec_by_element=True )
# Create the beam
if use_beam:
    beam_weight = prepare_weights( n_beam, beam_nx, beam_ny,
                                   beam_nz, dim, circ_m )
    beam = Species(type=Electron, weight=beam_weight, name='beam')
# Set the numerical parameters only now: they affect the newly created species
top.ssnpid = nextpid()
set_numerics( depos_order, efetch, particle_pusher, dim)

# Setup the field solver object
# -----------------------------
em = initialize_em_solver( stencil, dim,
    npass_smooth, alpha_smooth, stride_smooth,
    circ_m = (dim =="circ")*circ_m )
registersolver(em)

# Introduce the laser
# -------------------
if use_laser==1:
    add_laser( em, dim, laser_a0, laser_w0, laser_ctau, laser_z0,
        zf=laser_zfoc, theta_pol=laser_polangle, source_z=laser_source_z,
        laser_file=laser_file, laser_file_energy=laser_file_energy )

# Introduce the beam
# ------------------
# Load the beam
if use_beam:
    PlasmaInjector( beam, None, w3d, top, dim, beam_nx, beam_ny, beam_nz,
                beam_zmin, beam_zmax, beam_xmax, beam_ymax,
                dens_func = beam_dens_func, uz_m=beam_uz )
    initialize_beam_fields( em, dim, beam, w3d, top )

# Introduce the plasma
# --------------------
# Create an object to store the information about plasma injection
plasma_injector = PlasmaInjector( elec, ions, w3d, top, dim,
        plasma_nx, plasma_ny, plasma_nz, plasma_zmin,
        plasma_zmax, plasma_xmax, plasma_ymax, plasma_dens_func )
# Continuously inject the plasma, if the moving window is on
installuserinjection( plasma_injector.continuous_injection )

# Setup the diagnostics
# ---------------------
remove_existing_directory( ['diags'] )
# Whether to write the diagnostics in parallel
parallel_output = False
# Particle output
diag1 = FieldDiagnostic( period=diag_period, top=top, w3d=w3d, em=em,
            comm_world=comm_world, lparallel_output=parallel_output )
installafterstep( diag1.write )
# Field output
diag2 = ParticleDiagnostic( period=diag_period, top=top, w3d=w3d,
        species={species.name : species for species in listofallspecies},
        particle_data={"position","momentum","weighting","id"},
        comm_world=comm_world, lparallel_output=parallel_output )
installafterstep( diag2.write )

print('\nInitialization complete\n')

# ---------------
# Simulation loop
# ---------------
# Split loop for printing purposes
for n in range( int(N_steps/10) ):
    step(10)
