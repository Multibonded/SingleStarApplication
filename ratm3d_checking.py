import numpy as np
from itertools import islice
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def read_mesh(filename):
    """Read in the mesh files using a generator and spit out one value at a time.

    Parameters
    ----------
    filename : str
        The name of the mesh file.

    Returns
    -------
    val : str
        Yields the next value from the file per call.
    """

    with open(filename, 'r') as f:
        for line in f:
            for val in line.split():
                yield val

def ratm3d(atm3d_in, mesh_in = None, swap = False, rnh = False, rvturb = False, ismhd = False):
    """Read in the atmosphere in the atm3d and mesh files.

    Parameters
    ----------
    atm3d_in : str
        The path to the atm3d file.
    mesh_in : str
        The path to the mesh file.
    swap : Bool
        Swap the endian. Default is small.
    rnh : Bool
        Read the hydrogen level populations, n = 1-5 + HII.
    rvturb : Bool
        Read the microturbulence.
    ismhd : Bool
        Read magnitude magnetic field.

    Returns
    -------
    atm : dict
       The read in atmosphere.
    """

    if mesh_in is None:
        mesh = 'mesh.' + atm3d_in
        atm3d = 'atm3d.' + atm3d_in
    else:
        mesh = mesh_in
        atm3d = atm3d_in

    # mesh
    mesh_reader = read_mesh(mesh)
    nx = int(next(mesh_reader))
    x = np.array(list(map(float, (islice(mesh_reader, nx)))), dtype = 'float64')
    ny = int(next(mesh_reader))
    y = np.array(list(map(float, (islice(mesh_reader, ny)))), dtype = 'float64')
    nz = int(next(mesh_reader))
    z = np.array(list(map(float, (islice(mesh_reader, nz)))), dtype = 'float64')

    endian = '<'
    if swap:
        endian = '>'

    # snapshot
    with open(atm3d, 'rb') as f:
        count = nx*ny*nz
        nne = np.fromfile(f, dtype = endian + 'f', count = count)
        nne.shape = (nz, ny, nx)
        nne = nne.astype('float64')
        temp = np.fromfile(f, dtype = endian + 'f', count = count)
        temp.shape = (nz, ny, nx)
        temp = temp.astype('float64')
        vx = np.fromfile(f, dtype = endian + 'f', count = count)
        vx.shape = (nz, ny, nx)
        vx = vx.astype('float64')
        vy = np.fromfile(f, dtype = endian + 'f', count = count)
        vy.shape = (nz, ny, nx)
        vy = vy.astype('float64')
        vz = np.fromfile(f, dtype = endian + 'f', count = count)
        vz.shape = (nz, ny, nx)
        vz = vz.astype('float64')
        rho = np.fromfile(f, dtype = endian + 'f', count = count)
        rho.shape = (nz, ny, nx)
        rho = rho.astype('float64')
        if rnh:
            nh = np.fromfile(f, dtype = endian + 'f', count = count*6)
            nh.shape = (6, nz, ny, nx)
            nh = nh.astype('float64')
        else:
            nh = None
        if rvturb:
            vturb = np.fromfile(f, dtype = endian + 'f', count = count)
            vturb.shape = (nz, ny, nx)
            vturb = vturb.astype('float64')
        else:
            vturb = None
        if ismhd:
            bscal = np.fromfile(f, dtype = endian + 'f', count = count)
            bscal.shape = (nz, ny, nx)
            bscal = bscal.astype('float64')
        else:
            bscal = None

    atm3d = {'nx':nx, 'x':x, 'ny':ny, 'y':y, 'nz':nz, 'z':z,
             'nne':nne, 'temp':temp,
             'vx':vx, 'vy':vy, 'vz':vz,
             'rho':rho, 'nh':nh, 'vturb':vturb,
             'bscal':bscal}

    return atm3d

################################################################################
################################################################################
def XYcrossection(var, atm3d, layer):  #layer = 0 - 100, with 0 the lower boundary
    cp =plt.contourf(atm3d['x'],atm3d['y'], atm3d['temp'][layer],100,cmap = 'inferno')
    plt.colorbar(cp, label = var) #, ticks = [-0.5,0,0.5,1,2]
    plt.xlabel('$x/10^6$ ', fontsize = 30)
    plt.ylabel('$z/10^6$ ', fontsize = 30)
    plt.show()
    return

def Zcrossection(var, atm3d):
    Zsection = [0]*101
    for i in range(0,101):
        Zsection[i] = atm3d[var][i][20]

    plt.figure()
    cp = plt.contourf(atm3d['x']/10**6, atm3d['z']/10**6, Zsection ,100, cmap = 'inferno') #atm3d['rho'][0]
    plt.colorbar(cp, label = var) #, ticks = [-0.5,0,0.5,1,2]
    plt.xlabel('$x/10^6$ ', fontsize = 30)
    plt.ylabel('$z/10^6$ ', fontsize = 30)
    plt.show()
    return


def radialvar(N, atm3d, var):
    average = np.array([0.1]*101)
    nzz = np.array([[0.1]*N]*3600)  #3600rows*101columns
    for j in range(0,N):
        tempAvrg = 0
        for count, item in enumerate(nzz, start=0):
            k, u = (count//60, count-60*(count//60))
            nzz[count][j] = atm3d[var][j][k][u]
            tempAvrg = tempAvrg + nzz[count][j]
        average[j] = tempAvrg/3600
    print('start plotting')
    for b in range(0,N):
        plt.plot(-1*atm3d['z']/10**6,nzz[b],'k', linewidth = 0.5)
    plt.plot(-1*atm3d['z']/10**6, average, 'r--', linewidth = 3)
    return average

def radialvarMARCS(atm3d, var):
    nzz = np.array([0.1]*220)  #220 grid points
    for count, item in enumerate(atm3d[var], start=0):
        nzz[count] = item
    print('start plotting')
    plt.plot(-1*atm3d['z']/10**6,nzz,'c', linewidth = 2)
    return
#################################################################################

meshfileS = r"D:\PhD\TiFeAtoms\Jack\balder\atmospheres/mesh5792g+3.65z-2.36.sun"
atmfileS =  r"D:\PhD\TiFeAtoms\Jack\balder\atmospheres/atm3d5792g+3.65z-2.36.sun"


atmfileM =  r"D:\PhD\TiFeAtoms\Jack\balder\atmospheres/atm3d.hd122563" #220 grid points
meshfileM =  r"D:\PhD\TiFeAtoms\Jack\balder\atmospheres/mesh.hd122563"
atm3dS = ratm3d(atmfileS, meshfileS, swap = False, rnh = False, rvturb = False, ismhd = False)

atm3dM = ratm3d(atmfileM, meshfileM, swap = False, rnh = False, rvturb = False, ismhd = False)



RadialVar = ['temp', '$T$ [K]'] #radial variable to plot and label
XYVar = ['temp', '$T$ [K]'] #horizontal crosssection variable to plot and label
ZVar = ['vy', '$T$ [K]'] #vertical crosssection variable to plot and label

#################################################################################
plt.figure()
custom_legend = [mlines.Line2D([],[],color = 'k', linestyle = '-', linewidth = 1, label = 'Stagger'),
    mlines.Line2D([],[],color = 'r', linestyle = '--', linewidth = 1, label = 'Stagger avrg'),
    mlines.Line2D([],[],color = 'c', linestyle = '-', linewidth = 1, label = 'Marcs')]

radialvar(101,atm3dS, RadialVar[0])
radialvarMARCS(atm3dM, RadialVar[0])


plt.xlabel('radial coordinate [cm?]')
plt.ylabel(RadialVar[1])
plt.legend(handles = custom_legend)
plt.show()

XYcrossection(XYVar[0],atm3dS, 50)
Zcrossection(ZVar[0], atm3dS)
