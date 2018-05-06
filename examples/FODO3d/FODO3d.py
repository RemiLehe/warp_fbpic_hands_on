"""Basic 3D simulation of an ion beam in a periodic FODO lattice.
This input file sets up a periodic FODO lattice and creates a beam
that is matched to the lattice. The beam is propagated one lattice period.
"""
# --- This imports the Warp code into python, giving access to all
# --- of the Warp data and routines. This is typically the first command
# --- of a Warp input file.
from warp import *
from warp.init_tools import *
from warp.data_dumping.openpmd_diag import ParticleDiagnostic
from warp.data_dumping.openpmd_diag import ElectrostaticFieldDiagnostic

l_movieplot = False
l_movieplot3d = False

# --- Setup the description text which will be included at the bottom
# --- of every plot frame. This is for user convenience, documenting
# --- what the simulation is on the graphical output.
top.pline2   = "Example 3D beam in a FODO lattice"
top.pline1   = "Semi-Gaussian cigar beam. 32x32x128"
top.runmaker = "David P. Grote"

# --- Invoke plotting setup routine - it is needed to create a cgm output file for plots.
# --- The plot file will have a name with the format FODO3D.###.cgm. The prefix is the
# --- same as the input file name, and the ### is a number, increasing each time the
# --- simulation is carried out.
setup()

# --- Create the beam species. This instance of the Species class
# --- is used to configure aspects related to the beam and particles.
beam = Species(type=Potassium,charge_state=+1,name="Beam species")

# --- Set input parameters describing the beam, with a tune depression of 72 to 20 degrees.
# --- These numbers were generated by hand, by using the env package and adjusting the
# --- parameters to get the desired results.
# --- Note the units multipliers, e.g. mm. Almost all variables are MKS units
# --- (except ekin). The multipliers provide a nice way of converting to MKS
# --- while providing documentation about the units that are being used.
beam.a0    = 8.760439903086566*mm
beam.b0    = 15.599886448447793*mm
beam.emit  = 6.247186343204832e-05
beam.ap0   = 0.
beam.bp0   = 0.
beam.ibeam = 2.*mA
beam.vbeam = 0.
beam.ekin  = 80.*kV

# --- This call does some further processing on the input parameters.
# --- For example, in the above, the beam energy is specified. derivqty
# --- will calculate the beam velocity from the energy. It is not necessary
# --- to call derivqty (since it will be called during the generate), but it
# --- is called here to calculate beam.vbeam which is used in the calculation
# --- of the longitudinal thermal velocity spread below.
derivqty()

# --- Specify the longitudinal thermal velocity spread. In this case, it
# --- is the same as the transverse thermal velocity spread, as set by
# --- the emittance.
beam.vthz = 0.5*beam.vbeam*beam.emit/sqrt(beam.a0*beam.b0)

# --- Setup the FODO lattice
# --- These are user created python variables describing the lattice.
hlp     = 36.*cm   # half lattice period length
piperad = 3.445*cm # pipe radius
quadlen = 11.*cm   # quadrupole length

# --- Magnetic quadrupole field gradient - calculated to give sigma0 = 72 degrees.
dbdx    = 0.93230106124518164/quadlen

# --- Set up the quadrupoles. Only one lattice period is defined.
# --- This period is repeated to fill all space.
# --- The lattice consists of two quadrupoles, one defocusing, one focusing.
addnewquad(zs= 0.5*hlp - quadlen/2.,
           ze= 0.5*hlp + quadlen/2.,
           db=+dbdx)
addnewquad(zs=1.5*hlp - quadlen/2.,
           ze=1.5*hlp + quadlen/2.,
           db=-dbdx)

# --- zlatstrt is the start of the periodicity, relative to the quadrupoles position.
top.zlatstrt  = -4.*hlp

# --- zlatperi is the length of the lattice period, the length of the periodic repeat.
top.zlatperi  = 2.0*hlp

# ------------------------------------------------------------------------
# --- The next section sets up and run the envelope equation solver.
# --- Given the initial conditions specified above (a0, b0 etc.),
# --- the envelope package solves the KV envelope equations.
# --- The envelope solution will be used to specify the transverse
# --- shape of the beam where simulation particles will be loaded.

# --- The lattice period length, used to calculate phase advances.
top.tunelen = 2.*hlp

# --- The start and end of the envelope calculation. The initial conditions
# --- are the values at env.zl. Note that zl and zu must cover
# --- the longitudinal extent where the beam particles will be loaded.
# --- dzenv is the step size used in the envelope solver.
env.zl = -2.5*hlp # z-lower
env.zu = -env.zl  # z-upper

env.dzenv = top.tunelen/100.

# --- Select the envelope solver, do any initialization, and solve the equations.
package("env")
generate()
step()

# --- Make a plot of the resulting envelope solution.
penv()
fma()

# ------------------------------------------------------------------------
# --- Now, set up the parameters describing the 3D simulation.

# --- Specify the time step size. In this case, it is set so that
# --- it takes the specified number of time steps per lattice period.
steps_p_perd = 50
top.dt = (top.tunelen/steps_p_perd)/beam.vbeam

# --- Specify the number of grid cells in each dimension.
w3d.nx = 32
w3d.ny = 32
w3d.nz = 128

# --- Specify the extent of the field solve grid.
w3d.xmmin = -piperad
w3d.xmmax =  piperad
w3d.ymmin = -piperad
w3d.ymmax =  piperad
w3d.zmmin = -hlp*2
w3d.zmmax = +hlp*2

# --- Specify the boundary conditions on the outer sides of the grid.
# --- Possible values are dirichlet, periodic, and neumann.
w3d.bound0 = dirichlet # at iz == 0
w3d.boundnz = dirichlet # at iz == nz
w3d.boundxy = dirichlet # at all transverse sides

# --- Set the particle boundary conditions at the outer sides of the grid.
# --- Possible values are absorb, periodic, and reflect.
top.pbound0 = absorb
top.pboundnz = absorb
top.pboundxy = absorb

# --- Set the beam pulse length.
# --- Here, it is set to 80% of the grid length.
beam.zimin = w3d.zmmin*0.8
beam.zimax = w3d.zmmax*0.8

# --- Setup the parameters describing how the beam is created.
# --- npmax is the number of simulation particles to create.
top.npmax = 200000

# --- The distribution of the beam.
# --- There are a number of possible values, including "semigauss", "KV", and "WB".
w3d.distrbtn = "semigaus"
#w3d.distrbtn = "KV"

# --- The longitudinal velocity distribution of the beam. Possible values are "gaussian", "neuffer".
w3d.distr_l = "gaussian"
#w3d.distr_l = "neuffer"

# --- Turn on the "cigar" loading option This imposes a parabolic taper in the line-charge
# --- at the ends of the beam, adjusting the beam envelope to stay matched.
# --- beam.straight specifies the fraction of the beam that is in the middle, without the tapering.
# --- The length of each end will be (1 - beam.straight)/2.
w3d.cigarld = true
beam.straight = 0.5

# --- Set up field solver.
# --- fstype == 0 species the FFT solver.
top.fstype = 0

# --- Optional symmetries can be imposed on the solver.
# --- If l4symtry is true, the fields are calculated in only in transverse
# --- quadrant, and are replicated in the other quadrants.
# --- Note that the particles still occupy all of transverse space.
# --- When the charge is deposited, it would be mapped into the one quadrant.
w3d.l4symtry = false

# --- Setup various diagnostics and plots.
# --- By default, Warp calculates all 1st and 2nd order moments of the particles
# --- as a function of z position.

# --- Warp can save histories of the values of these moments at a selected number
# --- of z-locations relative to the beam frame. These locations are specified
# --- by zwindows. Note that the moments at the center of the window are saved.
# --- The zwindows are given a finite extent since that can also be used to
# --- select particles within the range, for plotting for example.
# --- Note that top.zwindows[:,0] always includes the whole longitudinal extent
# --- and should not be changed.
top.zwindows[:,1] = [-0.35, -0.3]
top.zwindows[:,2] = [-0.25, 0.25]
top.zwindows[:,3] = [0.3, 0.35]

# --- Since it can use a significant amount of memory, only time histories of the
# --- line-charge and vzbar are saved by default. These lines turn on the saving
# --- of time histories of other quantities.
top.lhxrmsz = true
top.lhyrmsz = true
top.lhepsnxz = true
top.lhepsnyz = true
top.lhcurrz = true

# --- nhist specifies the period, in time steps, of saving histories of
# --- the particle moments.
top.nhist = 1

# --- Define some plots to make and the frequency.
# --- zzplalways defines how often certain plots are generated.
# --- zzplalways = [zstart, zend, zperiod, extra_z_values, ...]
top.zzplalways[0:4] = [0.,100000.,2*hlp,0.]

# --- These specify that the plots ppzxy and ppzvz will be called as specified by zzplalways.
top.ipzxy[-2] = always
top.ipzvz[-2] = always

# --- User defined functions can be called from various points within the time step loop.
# --- This "@callfromafterstep" is a python decorator that says that this function will
# --- be called after every time step.
@callfromafterstep
def runtimeplots(nsteps=steps_p_perd):
    "Make user defined plots, every steps_p_perd steps"
    if top.it%nsteps != 0:
        return
    # --- Create overlaid plots in subframes of the plot window.
    plsys(9)
    pfzx(cellarray=1, contours=0, centering='cell')
    pzxedges(color=red, titles=False)
    plsys(10)
    pfzy(cellarray=1, contours=0, centering='cell')
    pzyedges(color=red, titles=False)
    fma()

    # --- Make plots of the transverse distribution in two zwindows.
    plsys(3)
    ppxy(iw=1)
    limits(-0.02, +0.02, -0.02, +0.02)
    plsys(4)
    ppxxp(iw=1)
    limits(-0.02, +0.02, -0.04, +0.04)
    plsys(5)
    ppxy(iw=3)
    limits(-0.02, +0.02, -0.02, +0.02)
    plsys(6)
    ppxxp(iw=3)
    limits(-0.02, +0.02, -0.04, +0.04)
    fma()


# --- Switch to the w3d package, which runs the 3D PIC model.
# --- The generate command does the initialization, including creating
# --- the particles, doing the initial Poisson solve, and calculating
# --- initial diagnostics and moments.
solverE=MultiGrid3D()
registersolver(solverE)
package("w3d")
generate()

# --- Define diagnostics
remove_existing_directory( ['diags'] )

particleperiod = 10
particle_diagnostic_0 = ParticleDiagnostic(period = particleperiod, top = top, w3d = w3d,
                                          species = {species.name: species for species in listofallspecies},
                                          comm_world=comm_world, lparallel_output=False)
fieldperiod = 10
efield_diagnostic_0 = ElectrostaticFieldDiagnostic(
    solver=solverE, top=top, w3d=w3d, comm_world = comm_world,
    period=fieldperiod, write_dir='diags/hdf5')

installafterstep(particle_diagnostic_0.write)
installafterstep(efield_diagnostic_0.write)

# --- Directly call the user defined function, producing plots of the initial conditions.
runtimeplots()

if l_movieplot:
    window(1,hcp='movie.cgm',dump=1)
    @callfromafterstep
    def movieplot():
        window(1)
        pfzx(view=9,titles=0,filled=0)
        pfzy(view=10,titles=0,filled=0)
        ppzxy(color=black)
        limits()
        fma();
        window(0)
        
if l_movieplot3d:
    try:
        from Opyndx import *
        iframe = 0
        os.system('rm -fr movie3d')
        os.system('mkdir movie3d')
        @callfromafterstep
        def movieplot3d():
            global iframe
            if top.it%3!=0:return
            pp,clm=viewparticles(getx()[::10],
                                 gety()[::10],
                                 getz()[::10]-ave(getz()),
                                 getr()[::10]*100,
                                 size=3.e-4,
                                 ratio=1.,
                                 color='auto',
                                 display=0,
                                 opacity=1.,
                                 cmin=0.,cmax=1.7)
            c = DXCamera(lookto=[0.,0.,0.],lookfrom=[0.00025484, 0.201298, -0.163537],width=100., \
                         resolution=640,aspect=0.75,up=[0.873401, 0.333893, 0.354521],perspective=1,angle=30., \
                         background='black')
            #DXImage(pp[0],scale=[1.,1.,0.1],camera=c)
            cc=[pp]
            cc.append(DXColorBar(clm,labelscale=1.,label='Radius [cm]',position=[0.01,0.95],min=0.,max=1.7))
            ccc = DXCollect(cc)
            im = DXScale(ccc,[1.,1.,0.1])
            DXWriteImage('movie3d/img%g'%(iframe+10000),im,c,format='tiff')
            iframe+=1
    except:
        pass


# --- Run for 50 time steps.
# --- Note that after each time step, the routine runtimeplots will be automatically called.
step(50)

# --- Make various post processing diagnostic plots.
ptitles('','','', 'Envelope is 2*Xrms, 2*Yrms')
ppgeneric(gridt=2.*top.hxrmsz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax,view=9)
limits()
ptitles('Beam X envelope history, in beam frame', 'Lattice periods', 'Beam frame location (m)')
ppgeneric(gridt=2.*top.hyrmsz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax,view=10)
limits()
ptitles('Beam Y envelope history, in beam frame', 'Lattice periods', 'Beam frame location (m)')
fma()

ppgeneric(gridt=top.hepsnxz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax,view=9)
limits()
ptitles('Beam slices X normalized emittance history, in beam frame', 'Lattice periods', 'Beam frame location (m)')
ppgeneric(gridt=top.hepsnyz[:,:top.jhist,0], xmin=0., xmax=top.zbeam/(2.*hlp), ymin=w3d.zmmin, ymax=w3d.zmmax,view=10)
limits()
ptitles('Beam slices Y normalized emittance history, in beam frame', 'Lattice periods', 'Beam frame location (m)')
fma()

hpepsnx(titles=0)
hpepsny(titles=0,color=red)
ptitles('X (black) and Y (red) emittances for whole beam','time [s]','!p-mm-mrad')
fma()

if l_movieplot:
    os.system("python cgmtomp4 movie.cgm 50")
    
if l_movieplot3d and iframe>0:
    step(150)
    os.system('cd movie3d;ffmpeg -y -r 12 -i img1%04d.tiff -strict -1 -qscale 1.0 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" FODO3D.mp4;mv FODO3D.mp4 ../.')
#    os.system('')
#    os.system('')
#    os.system('cd ..')
    
    
