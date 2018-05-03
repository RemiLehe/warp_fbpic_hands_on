"""
Python input script for a nonrelativistic Warp xy slice simulation of a 
K+ ion beam with intense space-charge focused by a hard-edge magnetic 
quadrupole doublet focusing lattice. A multi-grid field solver for the 
beam self-fields is employed.  This script is heavily commented to aid users 
in making modifications of this script for their specific problems.  More 
outputs and diagnostics can be easily added by uncommenting/modifying various 
lines.  For help and/or additional information contact:

     Dave Grote     dpgrote@lbl.gov   (510) 495-2961
     Steve Lund     smlund@lbl.gov    (510) 486-6937
     Jean-Luc Vay   jlvay@lbl.gov     (510) 486-4934

For more extensive documentation see the Warp web-site:

     https://warp.lbl.gov

All code inputs are mks units with the exception of the particle kinetic 
energy (eV).   

To run this Warp script:
  Interactive (return to interpreter after executing):

    % python -i xy-quad-mag-mg.py 

  Non-interactive (exit interpreter after executing):

    % python xy-quad-mag-mg.py
"""

# Load Warp and various script packages 
from warp        import *               # Warp code 
from warp.utils.errorcheck  import checksymmetry   # Check for input errors

# --- main run parameters and switches
e_kin = 100.*kV       # ion kinetic energy [eV] 
Q     = 1.e-4         # dimensionless perveance [1] 
emit  = 10.e-7        # rms edge emittance [m-rad]

l_automatch = False   # automatically matches beam initial conditions to lattice if True

n_grid = 200          # number grid cells (no symmetries) 

# Set informational labels included on all output cgm plots.   
top.pline2   = "xy Slice Simulation: Magnetic Quadrupole Doublet" 
top.pline1   = " "   # Add more info, if desired.  

# Invoke setup routine for graphics and output files (THIS IS MANDATORY)
setup()

# The run id will be the name of the input file. Other files created during
# a run, such as plot and dump files will have that name as a prefix.
# Plot files will also have a three digit number as a suffix. This is the run
# number and is incremented for each run, starting at 000. This suffix is
# stored in setup.pnumb. It is recommended that this suffix be added to any
# user created output files to keep the files connected.

# Set runmaker - it is included in informational labels on output plots
top.runmaker = "User Name"

# Beam parameters for simulation
#   Other than numerical parameters, the physical simulation is specified
#   by the numbers input immediately below.  

# --- Define species
beam = Species(type=Potassium, charge_state = +1, name = "Ion")

# --- Set initial beam envelope conditions.  Numbers given below are for an 
#     intial rms envelope matched tor the quadrupole focusing structure 
#     later defined. If the beam and/or focusing lattice is changed, 
#     matching conditions can be iterated by hand in (rapid) envelope code 
#     runs or a matching package based on a paper by:
#         Lund, Chilton, and Lee, PRSTAB 9, 064201 (2006)   
#     can be executed (requires sciPy installed) as follows: 
#
#   1) After setting up the focusing lattice and the beam, read in 
#      the envelope matching package: 
#        >>> from envmatch_KVinvariant import *    
#   2) Execute matching package and adjust values initial beam envelope 
#      radii r_x, r_y and initial beam envelope angles r_xp and r_yp 
#      consistent with returned values in beam.a0, beam.b0, etc.   
#        >>> Match()
#        >>> r_x  = beam.a0;  r_y  = beam.b0 
#        >>> r_xp = beam.ap0  r_yp = beam.bp0 
#
#     This procedure should be robust so long as the lattice phase advances 
#     in the absence of space-charge are less than 180 degrees/period. 

r_x  = 3.943*mm      # initial beam envelope in lattice period: r_x = a [m] 
r_y  = r_x           # initial beam envelope in lattice period: r_y = b [m]
r_xp = 12.86e-3      # initial beam envelope in lattice period: r_x' = a' [rad]
r_yp = -r_xp         # intiial beam envelope in lattice period: r_y' = b' [rad]

# --- factors to setup code inputs:
#     v_b  = axial velocity of nonrelativistic beam [m/sec] 
v_b = sqrt( 2.*jperev*e_kin/beam.mass )

##########################################################################################
# Input Warp parameters describing the beam.
##########################################################################################
#
# The beam kinetic energy (ekin) and axial velocity (vbeam) should not both
# be set unless done so consistently.  If one is zero the code sets from
# the other on generation:
#    vbeam = 0    => set vbeam from ekin
#    ekin  = 0    => set ekin  from vbeam
#
top.lrelativ  = 0            # turn off relativity 
beam.ekin      = e_kin       # kinetic energy of beam particle [eV]
beam.vbeam     = 0.          # beam axial velocity [m/sec]
beam.ibeam     = 2.*pi*eps0*beam.mass*v_b**3*Q/beam.charge
                             # beam current [amps]
beam.emitx     = emit        # beam x-emittance, rms edge [m-rad] 
beam.emity     = emit        # beam y-emittance, rms edge [m-rad]
beam.vthz      = 0.          # axial velocity spread [m/sec] 

# This routine will calculate vbeam and other quantities.
derivqty()

# Beam centroid and rms envelope initial conditions at s=0      
#
# --- centroid 
beam.x0  = 0.   # initial x-centroid xc = <x> [m]
beam.y0  = 0.   # initial y-centroid yc = <y> [m]
beam.xp0 = 0.   # initial x-centroid angle xc' = <x'> = d<x>/ds [rad]
beam.yp0 = 0.   # initial y-centroid angle yc' = <y'> = d<y>/ds [rad]
# --- envelope (must be set consistently for matched beam) 
beam.a0   =  r_x   # initial x-envelope edge a = 2*sqrt(<(x-xc)^2>) [m]
beam.b0   =  r_y   # initial y-envelope edge b = 2*sqrt(<(y-yc)^2>) [m]
beam.ap0  =  r_xp  # initial x-envelope angle ap = a' = d a/ds [rad]
beam.bp0  =  r_yp  # initial y-envelope angle bp = b' = d b/ds [rad]


##########################################################################################
# Setup linear applied focusing lattice for a piecewise constant 
##########################################################################################
# syncopated quadrupole doublet in the range s = z = [0,L_p].  The lattice is 
# extended periodically and z = 0 is taken (arbitrary choice) to be the 
# midpoint of the 1st drift length.  For notation employed, see:
#
#   USPAS Course Notes:  hifweb.lbl.gov/USPAS_2011
#   Lund and Bukh, PRSTAB 7, 024801 (2004) 
#   Lund, Chilton, and Lee, PRSTAB 9, 064201 (2006)   
#
#   L_p   = Lattice Period 
#   eta   = occupancy in (0,1]   
#             0 => Thin Lens,   1 => Full Occupancy  
#   alpha = syncopation parameter in [0,1] 
#             0   => Focus touches DeFocus to one   side 
#             1/2 => Symmetric FODO with equal drifts 
#             1   => Focus touches DeFocus to other side 
#   kappa_max = max value of focus function 
#                 set for desired phase advance (magnet strength)  
#
#   kappa_x(z) = - kappa_y(z)    Quadrupole symmetry 
#
#  If a specific value of lattice focusing strength is desired, this 
#  can be accoumplished by iterating the quadrupole gradient db by 
#  hand in successive (rapid) envelope code runs.  Or a special lattice 
#  rescaling package (requires sciPy installed) can be executed as follows:
#
#   1) After setting up lattice with a trial value of db in 
#      the quadrupole element load, read-in lattice rescale package 
#        >>> from lattice_rescale import *   
#   2) To obtain phase advance sigma0x = x-plane phase advance (degrees/period) 
#      run the rescale function below and then scale db by returned value.  
#        >>> rescalefunc(sigma0target=sigma0x,plane="x") 

L_p        = 0.5    # lattice period [m] 
eta        = 0.5    # Quadrupole occupancy [1] 
alpha      = 0.5    # Quadrupole syncopation factor [1] 
r_p        = 15.*mm # clear bore aperture (round beam pipe) [m]
kappa_max  = 50.49  # max kappa value of lattice [m^-2] 
                    #   (alternatively just set dbxdy field gradient below)  

# --- drift lengths between quadrupoles 
d_1 = alpha*(1-eta)*L_p       # 1st drift [m] (place 1/2 before 1st quad)
d_2 = (1-alpha)*(1-eta)*L_p   # 2nd drift [m]

# --- z starts of quadrupoles 
zs_q1 = d_1/2.                      # quad 1 [m]
zs_q2 = zs_q1 + (eta/2.)*L_p + d_2  # quad 2 [m]

# --- quadrupole field gradient 
dbxdy = beam.mass*beam.vbeam*kappa_max/beam.charge    # dbxdy = d B_y/dx [Tesla/m^2]

# --- Add Focus (db > 0) and DeFocus (db < 0) quadrupole elements.
addnewquad(zs = zs_q1, ze = zs_q1 + (eta/2.)*L_p, db =  dbxdy)
addnewquad(zs = zs_q2, ze = zs_q2 + (eta/2.)*L_p, db = -dbxdy) 

# --- Add grounded conducting pipe of radius r_p to grid 
aperture = ZCylinderOut(radius=r_p,length=largepos,voltage=0.,
                        zcent=L_p/2.,condid="next")

# --- Add a circular aperture particle scraper at pipe radius aperture. 
#       Could also use ParticleScraper(conductors=aperture) but 
#       setting prwall is faster for a simple cylinder. 
top.prwall = r_p 

# --- Lattice periodicity  
top.zlatperi  = L_p     # periodicity length [m]  
top.zlatstrt  = 0.      # z of lattice start; added to element z's [m] 
                        #   (can use to change lattice phase) 

##########################################################################################
# Define transverse simulation grid
##########################################################################################

# --- Symmetries.  Set for increased statistical efficiency.  These should
#     only be used in cases where lattice symmetries and initial beam
#     conditions and loads allow.  For no symmetry, set both options to false.
w3d.l2symtry = false     # 2-fold symmetry
w3d.l4symtry = true      # 4-fold symmetry

# -- Grid increments
#      First choose number grid cells without symmetry and reset 
#      consistent with symmetry options

sym_x = 1
sym_y = 1
if w3d.l4symtry:
  sym_x = 0.5
  sym_y = 0.5
elif w3d.l2symtry:
  sym_x = 0.5 

w3d.nx = int(sym_x*n_grid) 
w3d.ny = int(sym_y*n_grid)

# ---- Grid bounds 
#      Some bounds will be reset to zero by code on generation
#      if symmetry options are set.
l_grid = 2.*r_p               # length edge of simulation grid [m]      
w3d.xmmax =  l_grid/2.        # x-grid max limit [m] 
w3d.xmmin = -l_grid/2.        # x-grid min limit [m] 
w3d.ymmax =  l_grid/2.        # y-grid max limit [m] 
w3d.ymmin = -l_grid/2.        # y-grid min limit [m] 

# --- grid increments to use before code generation in setup
dx = l_grid/float(n_grid)


##########################################################################################
# Setup particle loading
##########################################################################################
#
# Set simulation macro-particle number (top.npmax) by specifying the
#  number of macro-particles to load per xy grid cell (nppg) using an
# rms equivalent uniform density beam measure.  This number is set
# consistently with the symmetry options.  
# 

nppg = 100    # number of particles per grid cell
top.npmax = int(nppg*pi*(beam.a0*beam.b0)/dx**2*sym_x*sym_y) # max initial particles loaded

# Distribution loads type 
# rms equivalent beam loaded with the specified distribution form 
#     KV => KV Distribution 
#
#     SG => semi-Gaussian distribution 
#             (KV density and local Gaussian angle spread about KV flutter) 
#
#     TE => Pseudoequilibrium with Thermal  Equilibrium form 
#     WB => Pseudoequilibrium with Waterbag Equilibrium form
#             The Pseudoequilibrium distributions use continuous focused 
#             equilibrium forms which are canoically transformed to AG 
#             symmetry of the lattice. 
#
#     For more info on loads, see review paper:
#       Lund, Kikuchi, and Davidson, PRSTAB 12, 114801 (2009) 

#w3d.distrbtn = "KV"          # initial KV distribution
#w3d.distrbtn = "TE"          # initial thermal distribution
#w3d.distrbtn = "WB"          # initial waterbag distribution
w3d.distrbtn = "SG"          # initial semi-Gaussian distribution 

# --- random number options to use in loading 
w3d.xrandom  = "digitrev"    # load x,y,z  with digitreverse random numbers 
w3d.vtrandom = "digitrev"    # load vx, vy with digitreverse random numbers
w3d.vzrandom = "digitrev"    # load vz     with digitreverse random numbers 
w3d.cylinder = true          # load a cylinder

##########################################################################################
# Setup particle pushing parameters
##########################################################################################
#
# Set increments and moving options by varying the number of steps to take
# per period (nstep) and the total number of periods to advance (nadvance).
#  Note:  Hard edge quadrupoles have
#         residence corrections allowing larger axial step-sizes.   
#

nstep    = 100     # number s-steps per lattice period
nadvance = 10.     # number lattice periods to advance (can be non-int)
ndiag    = 25      # number of moment accumulations per lattice period 
                   #   nstep/ndiag must be integer to accumulate at
                   #   the same phase in each period.  

top.lrelativ   =  false    # turn off relativistic kinematics
top.relativity = 0         # turn off relativistic self-field correction
                           #   to account for approx diamagnetic B-field of beam

wxy.ds = L_p/float(nstep)  # ds for part adv [m] 
wxy.lvzchang = true        # Use iterative stepping, which is needed if
                           # the vz of the particles changes.
                           #  ... must change even in linear lattice 
                           #          for high-order energy conservation 
top.ibpush   = 2           # magnetic field particle push, 
                           #   0 - off, 1 - fast, 2 - accurate 


# Setup field solver using 2d multigrid field solver. 

w3d.boundxy = 0              # Neuman boundary conditions on edge of grid.
w3d.solvergeom = w3d.XYgeom  # fieldsolve type to 2d multigrid 

# --- Uncomment to turn off space-charge deposition for simulation of particles 
#     moving in applied field  
#top.depos = "none"


# Turn on x-window plots, if desired; use winkill() to close interactively.  
#winon()


##########################################################################################
# Envelope calculation, beam matching and plotting
##########################################################################################

# Set envelope solver parameters  
env.zl      = 0.           # starting z for envelope calculation [m]
env.zu      = L_p          # ending   z for envelope calculation [m] 
env.dzenv   = L_p/1000.    # step size in envelope calculation [m] 

top.tunelen   = L_p     # length for calculation of sigma, sigma_0 [m]  

# Run envelope solver.  First the envelope package must be generated.  In
# the case of the envelope solver, a single call to step() advances the
# envelope from zl to zu.  The envelope solver can be run multiple times
# within a single session without problems (useful for iterating match
# conditions).
package("env");generate(); step()

# --- matches beam parameters to lattice using 15 iterations
if l_automatch:
    from warp.envelope.env_match import match2
    match2(15)

# --- generate envelope solver and run
generate(); step()

# Diagnostic plots of centroid and envelope evolution

# --- Plot beam centroid
plg( env.xenv/mm,env.zenv ) 
plg( env.yenv/mm,env.zenv, color="red")
ptitles("Env. Code: x [b] and y [r] Centroids",
        "s [m]","<x> and <y> [mm]",)
fma()

# --- Plot envelope as a function of z - superimpose quad focus strengths
plg( env.aenv/mm,env.zenv ) 
plg( env.benv/mm,env.zenv, color="red")
env_max = max(maxnd(env.aenv), maxnd(env.benv))/mm
env_min = min(minnd(env.aenv), minnd(env.benv))/mm
plg( env_min - 0.5*(env_max-env_min) +
     0.25*(env_max-env_min)*env.fqxenv/maxnd(env.fqxenv), env.zenv,
     color = "green")
ptitles("Env Code: RMS x [b] and y [r] Envelopes + Foc Func [g]",
        "s [m]","r_x_, r_y_  [mm]",)
fma()

# --- Plot envelope as a function of z - without quad focus strengths   
plg( env.aenv/mm,env.zenv )
plg( env.benv/mm,env.zenv, color = "red") 
ptitles("Env code: RMS x [b] and y [r] Envelopes",
        "s [m]","r_x_, r_y_  [mm]",)
fma() 

# Printout envelope quantities at the end of lattice
# (to check that match is consistent)  
print("a(0) - a(L_p) = ",beam.a0-env.aenv[env.nenv]," m") 
print("b(0) - b(L_p) = ",beam.b0-env.benv[env.nenv]," m") 
print("a'(0) - a'(L_p) = ",beam.ap0-env.apenv[env.nenv]," rad") 
print("b'(0) - b'(L_p) = ",beam.bp0-env.bpenv[env.nenv]," rad")


##########################################################################################
# Genrate 2-D XY package, install conductors and perform initial field solve
##########################################################################################

# Generate the xy PIC code.  In the generate, particles are allocated and
# loaded consistent with initial conditions and load parameters
# set previously.  Particles are advanced with the step() command later
# after various diagnostics are setup.
package("wxy"); generate()

# Install conducting aperture on mesh 
installconductors(aperture,dfill=largepos)

# Carry out explicit fieldsolve after generate to include conducing pipe 
# with initial beam 
fieldsolve()

# Check that inputs are consistent with symmetries (errorcheck package function)
checksymmetry()

##########################################################################################
# Setup diagnostics
##########################################################################################
# Diagnostics are grouped into several classes:
#   - Particle:  Snapshot plots of distribution function projections 
#   - Field:     Snapshot plots of self fields 
#   - History:   History plots on the evolution of moments and particle counts 
#                   accumulated as the simulation advances.   

# --- set max simulation step for diagnostic setup 
max_diag_step = int(1e5*nstep)

# --- set history diagnostic and moment accumulations 
top.nhist = nstep/ndiag                     # step interval for histories 
top.itmomnts[0:3] = [0,max_diag_step,top.nhist]   # do loop ranges for moments 
                                            #   and status writes to tty

# --- Plot limits for particle phase space plots. If lframe = true (default
#     false) diagnostics such as ppxxp for x-x' particle phase space will
#     use these ranges.  
#      max/min x,y   plot coordinates (m) 
#      max/min x',y' plot coordinates (rad)
l_diag = r_p
top.xplmax =  l_diag  
top.xplmin = -l_diag
top.yplmax =  l_diag
top.yplmin = -l_diag         
top.xpplmax = 2.5*max(beam.emitx/beam.a0,beam.emity/beam.b0) 
top.xpplmin = -top.xpplmax    
top.ypplmax =  top.xpplmax 
top.ypplmin = -top.xpplmax

# --- Color palette for phase-space plots (comment for default)
#     Search for .gp suffix files in the Warp scripts directory for possible
#     choices.  Some useful ones include:
#       earth.gp   (default)        heat.gp     (heat) 
#       gray.gp    (gray scale)     rainbow.gp  (rainbow) 
#palette("heat.gp")

# --- Set a chop factor for particle phase space plots to avoid plotting
#     too many particles (large storage and features will obscure).  Set
#     for approx 10 K particles per plot.  
chop_fraction = 10.e3/float(top.npmax) 

# --- Particle phase space diagnostics.
#     The list diag_step_part contains all steps where diagnostics in
#     diag_part() are made.  The list can contain repeated elements
#     and need not be ordered: take special" cases 1st period and then 
#     every 5 periods  
diag_step_part = ([0, int(nstep/4), int(nstep/2), int(3*nstep/4)] +
                  list(range(0,max_diag_step,5*nstep)))

def diag_part(plt_xy=False,plt_xxp=False,plt_yyp=False,plt_xpyp=False,
              plt_trace=False, plt_dens=False, out_dist=False):
  print("Making particle diagnostic plots")
  scale_fac = 1.e3/clight  # scale factor for canonical momentum [mrad]
  #
  # --- x-y projection
  if plt_xy:
    ppxy(lframe=true,chopped=chop_fraction,color='density',ncolor=25,
         titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("x-y Phase Space: Periods = %5.2f"%(top.zbeam/L_p),
            "x [mm]","y [mm]", )
    fma()
  # --- x-x' projection
  if plt_xxp: 
    ppxvx(lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
          titles=false,yscale=scale_fac,xscale=1./mm)
    ptitles("x-x' Phase Space: Periods = %5.2f"%(top.zbeam/L_p),
            "x [mm]","x' [mrad]", )
    fma()
  # --- y-y' projection
  if plt_yyp:
    ppyvy(lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
          titles=false,yscale=scale_fac,xscale=1./mm)
    ptitles("y-y' Phase Space: Periods = %5.2f"%(top.zbeam/L_p),
            "y [mm]","y' [mrad]", )
    fma()
  # --- x'-y' projection
  if plt_xpyp:
    ppvxvy(lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25,
           titles=false,yscale=scale_fac,xscale=scale_fac)
    ptitles("x'-y' Phase Space: Periods = %5.2f"%(top.zbeam/L_p),
            "x' [mrad]","y' [mrad]", )
    fma()
  # --- x-y, x-px, y-py, px-py projections, 4 to a page (trace-space)
  if plt_trace:
    pptrace(lframe=true,chopped=chop_fraction,slope='auto',color='density',ncolor=25)
    fma()
  # --- norm. number density on x and y axis
  if plt_dens:
    # --- density along principal axes and center, works for all symmetries. 
    ix_cen = sum(where(w3d.xmesh < 0.,1,0))
    iy_cen = sum(where(w3d.ymesh < 0.,1,0))
    rho_x = getrho(iy=iy_cen)
    rho_y = getrho(ix=ix_cen) 
    rho_max = maxnd(rho_x,rho_y) 
    if rho_max < beam.charge/(w3d.dx)**3: rho_max = 1 
    # 
    plg(rho_x/rho_max,1000*w3d.xmesh)
    if w3d.l4symtry: plg(rho_x/rho_max,-1000*w3d.xmesh) 
    plg(rho_y/rho_max,1000*w3d.ymesh,color="red")
    if w3d.l4symtry or w3d.l2symtry: 
      plg(rho_y/rho_max,-1000*w3d.ymesh,color="red")
    ptitles("Density on x [b] and y [r] Axes: Periods = %5.2f"%(top.zbeam/L_p),
            "x,y [mm]","Normalized Density [1]", )
    fma()      


# --- Field diagnostics.  
#     The list diag_step_field containins all steps where
#     diagnostics in diag_field() are made. The list can contain repeated
#     elements and need not be ordered.   
diag_step_field = ([0, int(nstep/4), int(nstep/2), int(3*nstep/4)] +
                   list(range(0,max_diag_step,5*nstep)))

def diag_field(plt_pc=False,plt_pc_xy=False):
  print("Making field diagnostic plots")
  #
  # --- self-field electrostatic potential
  if plt_pc:
    pfxy(cond=true,titles=false,yscale=1./mm,xscale=1./mm,iz = 0)
    ptitles("Self-Field Potential: Periods = %5.2f"%(top.zbeam/L_p),
            "x [mm]","y [mm]", )
    fma()
  # --- self-field electrostatic potential and particles together
  if plt_pc_xy:
    pfxy(cond=true,titles=false,yscale=1./mm,xscale=1./mm)
    ppxy(lframe=true,chopped=chop_fraction,color='density',ncolor=25,
         titles=false,yscale=1./mm,xscale=1./mm)
    ptitles("Self-Field Potential: Periods = %5.2f"%(top.zbeam/L_p),
            "x [mm]","y [mm]", )
    fma()

# --- History diagnostics.  These can be made at intermediate stages of the
#     run as well as at the end.  The list diag_step_hist contains all
#     steps where diagnostics in diag_hsit() are made. The list can
#     contain repeated elements and need not be ordered.
#     Notes:
#      * Many additional history diagnostics can be added by looking for
#        relevant moments accumulated in the Warp (see the variable group
#        "Hist" in top.v for an extensive list of variables that can be
#         used) and using gist commands to make relevant plots
#      * An example of this method can be found in the plt_emit_g option. 
diag_step_hist = list(range(10*nstep,max_diag_step,10*nstep))

def diag_hist(plt_cen=False,plt_env=False,plt_emit=False,plt_emit_g=False,
              plt_emitn=False,plt_pnum=False,plt_temp=False,plt_temp_flow=False):
  print("Making history diagnostic plots")
  #
  # --- centroid
  if plt_cen:
    hpxbar(titles=true,yscale=1./mm,xscale=top.vbeam/L_p)
    hpybar(titles=true,yscale=1./mm,xscale=top.vbeam/L_p,color="red")
    ptitles("History: Beam Centroid: x [b], y [r]",
            "s/L_p_, Lattice Periods","<x>, <y> Centroids [mm]", )
    fma()
  # --- envelope
  if plt_env:
    hpenvx(titles=false,yscale=1./mm,xscale=top.vbeam/L_p)    
    hpenvy(titles=false,yscale=1./mm,xscale=top.vbeam/L_p,color="red")
    ptitles("History: RMS Envelope: x [b], y [r]",
            "s/L_p_, Lattice Periods","RMS Edge Radii [mm]", )
    fma()
  # --- emittance
  scale_fac = 1.e6*top.vbeam/clight
  if plt_emit:
    hpepsx(titles=false,yscale=scale_fac,xscale=top.vbeam/L_p)
    hpepsy(titles=false,yscale=scale_fac,xscale=top.vbeam/L_p,color="red")
    ptitles("RMS Edge Emittance History: x [b], y [r]",
            "s/L_p_, Lattice Periods","RMS Edge Emittance [mm-mrad]", )
    fma()
  # --- emittance growth
  if plt_emit_g:
    if top.hepsx[0,0,0] > 0.: 
      plg(top.hepsx[0,0:top.jhist+1,0]/top.hepsx[0,0,0],
          top.hzbeam[0:top.jhist+1]/L_p )
    if top.hepsy[0,0,0] > 0.:
      plg(top.hepsy[0,0:top.jhist+1,0]/top.hepsy[0,0,0],
          top.hzbeam[0:top.jhist+1]/L_p,color="red")
      ptitles("History: Emittance Growth: x [b], y [r]",
              "s/L_p_, Lattice Periods","RMS Emittance Growth", )
    fma()      
  # --- particle number (to check for lost particles)
  if plt_pnum:
    hppnum(titles=false,xscale=top.vbeam/L_p)
    ptitles("History: Live Particle Number", "s/L_p_, Lattice Periods", 
            "Particle Number (simulation)", )
    fma()
  # --- temperature using <x'2> and corelation of x and x' via <x*x'> 
  if plt_temp and top.jhist > 0:
    xrms  = top.hxrms[0,1:top.jhist+1,0]
    emitx = top.hepsx[0,1:top.jhist+1,0]
    plg(beam.mass*top.vbeam**2*emitx**2/(jperev*16.*xrms**2),
        top.hzbeam[1:top.jhist+1]/L_p )
    yrms  = top.hyrms[0,1:top.jhist+1,0]
    emity = top.hepsy[0,1:top.jhist+1,0]
    plg(beam.mass*top.vbeam**2*emity**2/(jperev*16.*yrms**2),
        top.hzbeam[1:top.jhist+1]/L_p,color="red")
    ptitles("History: Spatial Avg Temp: x [b], y [r]",
            "s/L_p_, Lattice Periods","Spatial Avg Temp (eV)", )
    fma()

#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 
#  -- Install diagnostics at appropriate intervals after steps
#       Add options to generate plots desired 

# Function to call diagnostics at a timestep in step control lists 
def diag_calls():
  if top.it in diag_step_part:
    diag_part(plt_xy=true,plt_xxp=true,plt_yyp=true,plt_xpyp=false,
              plt_dens=true, plt_trace=false)
  if top.it in diag_step_field:
    diag_field(plt_pc_xy=true)
  if top.it in diag_step_hist:
    diag_hist(plt_env=true,plt_emit=true,plt_emit_g=false,plt_pnum=true,
              plt_temp=true)

# Install diagnostic calls after simulation step
installafterstep(diag_calls)

# Step 0 diagnostics (if any) of the initial distribution loaded 
diag_calls() 

# Advance simulation specified steps 
step(nstep*int(nadvance))

# Make additional history plots for final run if not already called 
if not(top.it in diag_step_hist):
  diag_hist(plt_env=true,plt_emit=true,plt_pnum=true,plt_temp=true)

# Save restart dump of run.  By default the name of the dump file is
# top.runid (or script name if this is not set) with the step number (iii)
# and ".dump" appended to the name:
#       runidiii.pdb 
# To restart:
#   % python
#     >>> from warp import *
#     >>> restart("runidiii.dump") 
#
#dump() 

# Print out timing statistics of run 
printtimers() 

# Make sure that last plot is flushed from buffer
fma() 

