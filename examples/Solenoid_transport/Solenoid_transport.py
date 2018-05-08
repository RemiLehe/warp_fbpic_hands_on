"""
Example Pierce diode with subsequent solenoid transport.
Hot plate source emitting singly ionized potassium.
"""
from warp import *

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Pierce diode with solenoid transport example"
top.pline1   = "Injected beam. Semi-Gaus."
top.runmaker = "DPG, JLV"

# --- Invoke setup routine for the plotting
setup()

# ------------------------------------------------------------------------------
# main simulation parameters
# ------------------------------------------------------------------------------

# --- turns Solenoid on/off
l_solenoid = False

# --- Set the dimensionality
w3d.solvergeom = w3d.RZgeom
#w3d.solvergeom = w3d.XYZgeom

# --- if true, inject one particle at each grid node (random injection otherwise)
w3d.l_inj_regular = False

# --- Approximate number of particles injected each step
top.npinject    = 45                     

# ------------------------------------------------------------------------------
# setup source
# ------------------------------------------------------------------------------

# --- Basic parameters
channel_radius = 15.*cm
diode_voltage = 93.*kV

# --- Setup source plate
source_radius = 5.5*cm
source_temperature = 0.1 # in eV
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

conductors = source+plate

# --- Specify injection of the particles
top.inject      = 2                       # 2 means space-charge limited injection
top.rinject     = source_curvature_radius # Source radius of curvature

top.vinject     = diode_voltage           # Voltage on the injection source
w3d.l_inj_exact = True                    # if true, position and angle of injected particle are 
                                          # computed analytically rather than interpolated
w3d.l_inj_rz = (w3d.solvergeom == w3d.RZgeom) # If using the RZ geometry, set so injection uses the same geometry

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
if l_solenoid:
   runlen = 1.
else:
   runlen = zplate + 5.*cm
   
# --- Variables to set symmetry, when using 3D
#w3d.l4symtry = true
#w3d.l2symtry = false

# --- Set boundary conditions
# ---   for field solve
w3d.bound0  = dirichlet
w3d.boundnz = neumann
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
w3d.nx = w3d.ny = 38
if l_solenoid:
    w3d.nz = 250
else:
    w3d.nz = 64

# --- Set the time step size. This needs to be small enough to satisfy the Courant limit.
dz = (w3d.zmmax - w3d.zmmin)/w3d.nz
vzfinal = sqrt(2.*diode_voltage*jperev/beam.mass)
top.dt = 0.4*(dz/vzfinal)

# ------------------------------------------------------------------------------
# Set up solenoid lattice
# ------------------------------------------------------------------------------
if l_solenoid:
    match_solenoid_length = 0.07
    match_solenoid_Bpeak = 2.8
    solenoid_radius = 8.*cm
    solenoid_length = 0.3
    solenoid_period = 0.35
    solendoid_halfgap = (solenoid_period - solenoid_length)/2.
    solenoid_Bpeak = 2.

    addnewsolenoid(zi=zplate+solendoid_halfgap,
                   zf=zplate+solendoid_halfgap+match_solenoid_length,
                   ri=solenoid_radius,
                   maxbz=match_solenoid_Bpeak)

    nsolenoids = 5
    z0 = zplate + solendoid_halfgap + match_solenoid_length + solendoid_halfgap
    for i in range(nsolenoids):
        addnewsolenoid(zi=z0 + solendoid_halfgap+i*solenoid_period,
                       zf=z0 + solendoid_halfgap+i*solenoid_period + solenoid_length,
                       ri=solenoid_radius,
                       maxbz=solenoid_Bpeak)

    def plotsolenoids():
        cc = array([10],dtype=ubyte)
        rr = [solenoid_radius, solenoid_radius, solenoid_radius+1.*cm, solenoid_radius+1.*cm]
        zz = zplate+solendoid_halfgap + array([0., match_solenoid_length, match_solenoid_length, 0.])
        plfp(cc, rr, zz, [4])

        z0 = zplate + solendoid_halfgap + match_solenoid_length + solendoid_halfgap
        zz = array([0., solenoid_length, solenoid_length, 0.])
        for i in range(nsolenoids):
            plfp(cc, rr, z0 + solendoid_halfgap+i*solenoid_period + zz, [4])

    # --- Pipe in the solenoid transport
    pipe = ZCylinderOut(radius=solenoid_radius, zlower=zplate+plate_width/2., zupper=zplate+nsolenoids*solenoid_period)
    conductors+=pipe

# ------------------------------------------------------------------------------
# setup field solver, install conductors & scraper
# ------------------------------------------------------------------------------
f3d.mgtol = 1.e-1 # Multigrid solver convergence tolerance, in volts

solver = MultiGrid2D()
registersolver(solver)

installconductor(conductors)

# --- Setup the particle scraper
#scraper = ParticleScraper([source, plate])#, pipe])
scraper = ParticleScraper(conductors)

# ------------------------------------------------------------------------------
# load main package and initialize 
# ------------------------------------------------------------------------------

# --- Set pline1 to include appropriate parameters
top.pline1 = ("Injected beam. Semi-Gaus. %dx%d. npinject=%d, dt=%d"%(w3d.nx, w3d.nz, top.npinject, top.dt))

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")
generate()

# ------------------------------------------------------------------------------
# setup plots
# ------------------------------------------------------------------------------

# --- Open up plotting windows
winon(1, suffix='current')
winon()

def rzplot(view=1):
    plsys(view)
    pfzr(plotsg=0, cond=0, titles=False, view=view)
    source.draw(filled=150, fullplane=False)
    plate.draw(filled=100, fullplane=False)
    if l_solenoid:
        pipe.draw(filled=100, fullplane=False)
        plotsolenoids()
    ppzr(titles=False, view=view)
    limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., channel_radius)
    ptitles('', 'Z (m)', 'R (m)')
    
def beamplots():
    window(0)
    fma()
    rzplot()
    plsys(1)
    ptitles('Hot plate source into solenoid transport', '','')
    refresh()

    window(1)
    fma()
    pzcurr()
    limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., diode_current*1.5)
    refresh()

# --- Call beamplots after every 20 steps
@callfromafterstep
def makebeamplots():
    if top.it%20 == 0:
        beamplots()

step(1000)

# --- Make sure that last plot frames get sent to the cgm file
window(0)
hcp()
window(1)
hcp()
