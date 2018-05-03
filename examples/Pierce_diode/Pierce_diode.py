"""
Example Pierce diode calculation.
Hot plate source emitting singly ionized potassium
"""
from warp import *
from warp.run_modes.egun_like import gun
from warp.utils.constantcurrentinjection import SpecifiedCurrentRiseTime

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Pierce diode example"
top.pline1   = "Injected beam. Semi-Gaus."
top.runmaker = "DPG, JLV"

# --- Invoke setup routine for the plotting
setup()

# ------------------------------------------------------------------------------
# main simulation parameters
# ------------------------------------------------------------------------------

# --- Set the dimensionality
w3d.solvergeom = w3d.RZgeom
#w3d.solvergeom = w3d.XYZgeom

# --- Sets method of running
# ---   Steady state gun mode
# ---   Time dependent simulation (when False)
steady_state_gun = False

# --- Automatically adjust diode voltage to produce constant current if True
l_constant_current = False

# --- if true, inject one particle at each grid node (random injection otherwise)
w3d.l_inj_regular = False   

# --- Approximate number of particles injected each step
top.npinject    = 150                     

# --- type of injection
top.inject      = 1    # 1=constant current, 2=space-charge limited injection
                       # type doc('inject') for more options
     
# ------------------------------------------------------------------------------
# setup source
# ------------------------------------------------------------------------------

# --- Basic parameters
channel_radius = 15.*cm
diode_voltage = 93.*kV

# --- Setup source plate
source_radius = 5.5*cm
source_temperature = 1.e-6 # in eV
source_curvature_radius = 30.*cm # --- radius of curvature of emitting surface
pierce_angle = 67.

# --- Setup diode aperture plate
zplate = 8.*cm       # --- plate location
rplate = 5.5*cm      # --- aperture radius
plate_width = 2.5*cm # --- thickness of aperture plate

# --- Create aperture plate
rround = plate_width/2.
plate = ZRoundedCylinderOut(radius=rplate, length=plate_width, radius2=rround, voltage=0., zcent=zplate)

piercezlen = (channel_radius - source_radius)*tan((90.-pierce_angle)*pi/180.)
piercezlen = 0.04

# --- Create source conductors

# --- Outer radius of Pierce cone
rpierce = source_radius + piercezlen*tan(pierce_angle*pi/180.)

# --- Depth of curved emitting surface
sourcezlen = (source_radius**2/(source_curvature_radius + sqrt(source_curvature_radius**2 - source_radius**2)))

# --- the rsrf and zsrf specify the line in RZ describing the shape of the source and Pierce cone.
# --- The first segment is an arc, the curved emitting surface.
source = ZSrfrv(rsrf = [0.,                      source_radius, channel_radius,          channel_radius],  # surface of revolution points R positions
                zsrf = [0.,                      sourcezlen,    sourcezlen + piercezlen, 0.            ],  # surface of revolution points Z positions
                zc   = [source_curvature_radius, None,          None,                    None          ],  # surface of revolution arc center R position
                rc   = [0.,                      None,          None,                    None          ],  # surface of revolution arc center Z position
                voltage=diode_voltage)

# --- create collector
collector = Plane(z0=12.5*cm,name='collector')

# --- Specify injection of the particles
top.rinject     = source_curvature_radius # Source radius of curvature
top.vinject     = diode_voltage           # Voltage on the injection source
w3d.inj_nz = 0
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are 
                                          # computed analytically rather than interpolated
w3d.l_inj_rz = (w3d.solvergeom == w3d.RZgeom) # If using the RZ geometry, set so injection uses the same geometry
w3d.l_inj_addtempz_abs = True             # Sets initial velocity spread only in the forward direction

# ------------------------------------------------------------------------------
# setup beam
# ------------------------------------------------------------------------------

# --- Setup simulation species
beam = Species(type=Potassium, charge_state=+1, name='beam')

# --- Child-Langmuir current between parallel plates
j = 4./9.*eps0*sqrt(2.*echarge*beam.charge_state/beam.mass)*diode_voltage**1.5/zplate**2
diode_current = pi*source_radius**2*j

print(("Child-Langmuir current density = ", j))
print(("Child-Langmuir current = ", diode_current))

# --- Set basic beam parameters
beam.a0       = source_radius
beam.b0       = source_radius
beam.ap0      = .0e0
beam.bp0      = .0e0
beam.ibeam    = diode_current
beam.vthz     = sqrt(source_temperature*jperev/beam.mass)
beam.vthperp  = sqrt(source_temperature*jperev/beam.mass)
derivqty()

# ------------------------------------------------------------------------------
# setup grid and boundary conditions
# ------------------------------------------------------------------------------

# --- Length of simulation box
runlen = zplate + 5.*cm

# --- Variables to set symmetry, when using 3D
#w3d.l4symtry = true
#w3d.l2symtry = false

# --- Set boundary conditions
# ---   for field solve
w3d.bound0  = dirichlet
w3d.boundnz = dirichlet
w3d.boundxy = neumann
# ---   for particles
top.pbound0  = absorb
top.pboundnz = absorb
top.prwall   = channel_radius

# --- Set field grid size
w3d.xmmin = w3d.ymmin = -channel_radius
w3d.xmmax = w3d.ymmax = +channel_radius
w3d.zmmin = 0.
w3d.zmmax = runlen

# --- Field grid dimensions - note that nx and ny must be even.
w3d.nx = w3d.ny = 64
if w3d.solvergeom == w3d.RZgeom:w3d.nx/=2
w3d.nz = 32

# --- Set the time step size. This needs to be small enough to satisfy the Courant limit.
dz = (w3d.zmmax - w3d.zmmin)/w3d.nz
vzfinal = sqrt(2.*diode_voltage*jperev/beam.mass)
top.dt = 0.4*(dz/vzfinal)

# ------------------------------------------------------------------------------
# setup field solver, install conductors & scraper
# ------------------------------------------------------------------------------
f3d.mgtol = 1.e-1 # Multigrid solver convergence tolerance, in volts

if w3d.solvergeom == w3d.XYZgeom:
    solver = MultiGrid3D()
else:
    solver = MultiGrid2D()
registersolver(solver)

# --- install conductors+collector
installconductor(source+plate+collector, dfill=largepos)

# --- Setup the particle scraper
scraper = ParticleScraper([plate,collector],lcollectlpdata=True)

# ------------------------------------------------------------------------------
# load main package and initialize 
# ------------------------------------------------------------------------------

# --- Set pline1 to include appropriate parameters
if w3d.solvergeom == w3d.RZgeom:
    top.pline1 = ("Injected beam. Semi-Gaus. %dx%d. npinject=%d, dt=%d"%
                  (w3d.nx, w3d.nz, top.npinject, top.dt))
else:
    top.pline1 = ("Injected beam. Semi-Gaus. %dx%dx%d. npinject=%d, dt=%d"%
                  (w3d.nx, w3d.ny, w3d.nz, top.npinject, top.dt))

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")
generate()

# ------------------------------------------------------------------------------
# setup plots
# ------------------------------------------------------------------------------

# --- Open up plotting windows
if not steady_state_gun:
    window(1, hcp='current.cgm',dump=1)
winon()

def beamplots(l_plotcurrent=True):
    window(0)
    fma()
    if l_plotcurrent:
        plsys(10)
        pzcurr()
        limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., diode_current*1.5)
        view=9
        xlim=channel_radius/2
        pfzr(plotsg=0, cond=0, titles=False, filled=1, contours=100,view=view)
    else:
        view=1
        xlim=channel_radius
    if w3d.l_inj_regular:
        msize=2
    else:
        msize=1
    if w3d.solvergeom == w3d.XYZgeom:
        pfzx(plotsg=0, cond=0, titles=False)
        source.draw(filled=150)
        plate.draw(filled=100)
        ppzx(titles=False,msize=msize, color=red, view=view)
        limits(w3d.zmminglobal, w3d.zmmaxglobal, -xlim, xlim)
        ptitles('Hot plate source', 'Z (m)', 'X (m)')
    else:
        pfzr(plotsg=0, cond=0, titles=False,view=view)
        source.draw(filled=150, fullplane=False)
        plate.draw(filled=100, fullplane=False)
        ppzr(titles=False,msize=msize, color=red, view=view)
        limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., xlim)
        ptitles('Hot plate source', 'Z (m)', 'R (m)')
    refresh()

# --- Turn on optional automatic adjustment of diode voltage for specific current history
if l_constant_current and not steady_state_gun:
    def emittedcurrentprofile(time):
        return j
    TDC = SpecifiedCurrentRiseTime(sourceid=source,
                                  currentdensityfunc=emittedcurrentprofile,
                                  sourcevolt=diode_voltage,
                                  otherids=[],
                                  othervolts=[],
                                  endplatevolt=0.)#

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

if steady_state_gun:
    # --- Steady-state operation
    # --- This does steady-state gun iterations, plotting the z versus r
    # --- after each iteration.
    gun(10, ipstep=1, insertafteriter=beamplots, l_savepart_always=True)

else:

    # --- Call beamplots after every 20 steps
    @callfromafterstep
    def makebeamplots():
        if top.it%20 == 0:
            beamplots(True)

    step(700)

# ------------------------------------------------------------------------------
# Final diagnostics
# ------------------------------------------------------------------------------

if not steady_state_gun:
    window(1)
    plsys(9)
    if l_constant_current:
        pla(TDC.hsourcevolt[...]*1.e-3,TDC.htime[...]*1.e6,color=red,width=3)
    else:
        pldj([0.],[diode_voltage*1.e-3],[top.time*1.e6],[diode_voltage*1.e-3],color=red,width=3)
    ptitles('Voltage history on emitter','Time [microsec.]','Emitter voltage [kV]')
    limits(0.,top.time*1.e6,0.,150.)
    plsys(10)
    collector.plot_current_history()

# --- Make sure that last plot frames get sent to the cgm file
window(0)
hcp()
if not steady_state_gun:
    window(1)
    hcp()
    
    