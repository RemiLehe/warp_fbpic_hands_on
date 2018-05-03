"""
Launches laser at in vacuum or with dielectric medium.
"""
import warpoptions    # --- first import warpoptions module
# --- add arguments
warpoptions.parser.add_argument('-ll','--lambda',  dest='lambda_laser',
                                type=float,default=1.e-6,help='laser wavelength')
warpoptions.parser.add_argument('-er','--er',      dest='er',
                                type=float,default=1.,   help='relative permittivity')
warpoptions.parser.add_argument('-rf','--rfact',   dest='rfact',
                                type=float,default=1.,   help='resolution scaling factor')
 
warpoptions.parser.add_argument('-lp','--laserpos',dest='laser_position',
                                type=str,default='[0.,0.,-5.e-6]',help='center of laser emitter')
warpoptions.parser.add_argument('-lv','--laservec',dest='laser_vector',
                                type=str,default='[0.,0.,1.]',help='laser emission vector')
warpoptions.parser.add_argument('-lo','--laserpol',dest='laser_polvector',
                                type=str,default='[0.,1.,0.]',help='laser polarization vector')

warpoptions.parser.add_argument('-bp','--boxposx',  dest='boxposx',
                                type=str,default='[2.5e-6,7.5e-6]',help='box min & max in x')

from warp import *    # --- import Warp

list_options = ['lambda_laser','er','rfact','laser_position','laser_vector']
list_options += ['laser_polvector','boxposx']
for o in list_options:
    ol = o
    exec(o+' = warpoptions.options.'+o)
    exec('t = type('+ol+')')
    if t is str:
        exec('newval = '+o)
        exec('newval = '+newval)
        exec(o+' = asarray(newval)')

from warp.field_solvers.laser import *
from warp.field_solvers.laser.laser_profiles import GaussianProfile

dim = "2d"

l_test = 1   # --- open window on screen if true, save on disk in cgm file otherwise
ncells = 256 # --- nb cells in x

# --- grid dimensions and nb cells
w3d.xmmin = -10.e-6             # --- min x of simulation box
w3d.xmmax =  10.e-6             # --- max x of simulation box
w3d.zmmin = -10.e-6             # --- min z of simulation box
w3d.zmmax = w3d.zmmin+50.e-6    # --- max z of simulation box
w3d.nx    = nint(256*rfact)     # --- nb cells in x
w3d.nz    = nint(640*rfact)     # --- nb cells in z
if dim=='2d':
  # --- sets y min/max and ny so that dy=1 in 2D
  w3d.ymmin = -1.
  w3d.ymmax = 1.
  w3d.ny = 2
else:
  # --- sets y min/max and ny to x values in 3D
  w3d.ymmin = w3d.xmmin
  w3d.ymmax = w3d.xmmax
  w3d.ny = ncells
w3d.dx = dx = (w3d.xmmax-w3d.xmmin)/w3d.nx
w3d.dy = dy = (w3d.ymmax-w3d.ymmin)/w3d.ny
w3d.dz = dz = (w3d.zmmax-w3d.zmmin)/w3d.nz

# --- field boundary conditions
bounds = zeros(6)
bounds[:] = openbc

k0                 = 2.*pi/lambda_laser
w0                 = k0*clight
a0                 = 3.                         # normalized potential vector (amplitude)
laser_waist        = 0.1*(w3d.xmmax-w3d.xmmin)  # laser waist
laser_duration     = 12./w0
laser_source_v     = 0.
W0                 = laser_waist        # Radius of beam at waist
Zr                 = (pi*laser_waist**2)/lambda_laser # Rayleigh length  
Zf                 = 1.*Zr                  # position of laser focus (relative to antenna)
laser_emax         = a0*w0*emass*clight/echarge

laser = GaussianProfile(k0, 
                        laser_waist, 
                        laser_duration, 
                        4*laser_duration, 
                        a0, 
                        dim,
                        focal_length=-Zf, 
                        temporal_order=2, 
                        boost=None,
                        source_v=0)

# --- add dielectric
box_xmin = boxposx[0]
box_xmax = boxposx[1]
dielectric_box = Box(zsize=w3d.zmmax-w3d.zmmin, \
                     ysize=w3d.ymmax-w3d.ymmin, \
                     xsize=box_xmax-box_xmin, \
                     xcent=0.5*(box_xmin+box_xmax), \
                     zcent=0.5*(w3d.zmmin+w3d.zmmax), \
                     permittivity=er)

# --- open graphics window (on screen if l_test=1, on disk otherwise)
if l_test:
  window(0,dpi=100)
else:
  setup()
  
palette('bluewhitered.gp')

ions = Species(type=Proton) # --- unused but necessary to declare at least one particle species
top.fstype=-1               # --- turns off electrostatic solver
generate()                  # initializes internal arrays

#-------------------------------------------------------------------------------
# --- initialization of electromagnetic solver
#-------------------------------------------------------------------------------
em = EM3D(bounds=bounds,
          l_2dxz=dim=='2d',
          stencil=1)
  
registersolver(em)
installconductors(dielectric_box)
em.finalize()

# --- add one pass bilinear smoothing (of currents) in x and z
em.npass_smooth = array([[1],[0],[1]])

# --- add laser antenna
laser_antenna = LaserAntenna(laser, 
                             laser_vector,
                             laser_polvector, 
                             laser_position,
                             laser_emax, 
                             None, 
                             laser_source_v,
                             polangle=None, 
                             w3d=w3d, 
                             dim=dim, 
                             circ_m=em.circ_m)
em.laser_antenna.append( laser_antenna )

def rotate(x,c,angle):
    cr = cos(angle)
    sr = sin(angle)
    return c + (x-c)*array([cr,cr])+(x-c)[::-1]*array([-sr,sr])

#-------------------------------------------------------------------------------
# --- definition and installation of plotting routine
#-------------------------------------------------------------------------------
def mkplots():
    print em.get_tot_energy()
    if top.it%10!=0: return # --- plots every 10 time steps
    window(0);
    # --- 2D plot of Ey at top of frame
    if laser_polvector[0]==0.:
        pf = em.pfey
    else:
        pf = em.pfex
    fma();pf(direction=1,l_transpose=1,view=9,cmin=-laser_emax,cmax=laser_emax);
    
    # --- trace emission plane
    angle = arctan2(laser_vector[0],laser_vector[2])    
    ptc = array([laser_position[2],laser_position[0]])
    pt1 = array([laser_position[2],w3d.xmmin*10])
    pt2 = array([laser_position[2],w3d.xmmax*10])
    
    pt1 = rotate(pt1,ptc,angle)
    pt2 = rotate(pt2,ptc,angle)
    
    pldj([pt1[0]],[pt1[1]],[pt2[0]],[pt2[1]],color=green)

    # --- trace focus plane
    ptc = array([laser_position[2]+Zf,laser_position[0]])
    pt1 = array([laser_position[2]+Zf,w3d.xmmin*10])
    pt2 = array([laser_position[2]+Zf,w3d.xmmax*10])
    
    pt1 = rotate(pt1,ptc,angle)
    pt2 = rotate(pt2,ptc,angle)
    
    pldj([pt1[0]],[pt1[1]],[pt2[0]],[pt2[1]],color=magenta)

    # --- if dielectric, trace dielectric box contours
    if er>1.:
        dielectric_box.draw()
        
    limits(w3d.zmmin,w3d.zmmax,w3d.xmmin,w3d.xmmax)

    # --- refresh window
    refresh()
installafterstep(mkplots)
 
# inject laser for some time
step(nint(15/2*laser_duration/top.dt))  

if er==1. and all(laser_vector==array([0.,0.,1.])): 
    # --- if no dielectric and propagating along z, then compare to paraxial solution

    # --- now plot electric field from analytical formula under paraxial approximation 

    x0,z0 = getmesh2d(w3d.xmmin,w3d.dx,w3d.nx-1,
                      w3d.zmmin,w3d.dz,w3d.nz-1)

    Zc = laser_position[2]

    z0 -= Zc

    zd = z0.copy()-Zf

    z0 += 4*laser_duration*clight
    z0 -= top.time*clight

    diffract_factor = 1. + 1j*(zd)/Zr

    # --- Calculate the argument of the complex exponential
    exp_argument = 1j * k0*z0  \
                 - (x0**2) / (diffract_factor*laser_waist**2) \
                 - (z0 / (laser_duration*clight))**2

    if dim=='2d':geom_factor=0.5
    if dim=='3d':geom_factor=1.
    profile = exp(exp_argument) / ( diffract_factor**geom_factor )

    E = laser_emax * profile.imag

    ppg(transpose(E),view=10,cmin=-laser_emax,cmax=laser_emax,xmin=w3d.zmmin,xmax=w3d.zmmax,ymin=w3d.xmmin,ymax=w3d.xmmax)
    pldj([Zc],[w3d.xmmin],[Zc],[w3d.xmmax],color=green)
    pldj([Zc+Zf],[w3d.xmmin],[Zc+Zf],[w3d.xmmax],color=magenta)
    ptitles('E_y ^^analytic','z','x')
    refresh()

