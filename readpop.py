import numpy as np
import numpy as np
import struct
import numpy as np
import matplotlib.pyplot as plt
import const
import sys
from readsp import readsp
from readdet import readdet
from os import listdir
from ratom import ratom
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Radiative_Run\output"
files = listdir(PathToOutput)
for f in files:
	if "stagger" in f:
		hydro = "stagger"
	elif "marcs" in f:
		hydro = "marcs"
filename_det = PathToOutput + "/detailed_"+hydro+"_sun"
filename_sp = PathToOutput +"/sp"
file_iec = PathToOutput +"/iec_"+hydro+"_sun"
file_iecl = PathToOutput +"/iecl_"+hydro+"_sun"
file_iet = PathToOutput +"/iet_"+hydro+"_sun"
file_ietl = PathToOutput +"/ietl_"+hydro+"_sun"
atom=ratom(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\atoms"+'/atom.ti_Nq140k_k0001')
sp = readsp(filename_sp)
out = readdet(filename_det, sp)

from astropy.io import fits
atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"

atom11 = fits.open(atomname)
atom1 = atom11[1].data
atom2 = atom11[2].data

def readpop(filename, nx, ny, nz, nlevel):
    """Reads population.

    Parameters
    ----------
    filename : str
        Path to binary file.
    nx : int
        x dimension. From readdet.
    ny : int
        y dimension. From readdet.
    nz : int
        z dimension. From readdet.
    nlevel : int
        number of levels. From ratom.

    Returns
    -------
    pop : dict
        The population data.
    """
    with open(filename, 'rb') as f:
        b = np.fromfile(f, dtype = 'float32', count = nx*ny*nz*nlevel)
        b.shape = (nlevel, nz, ny, nx)
        nstar = np.fromfile(f, dtype = 'float32', count = nx*ny*nz*nlevel)
        nstar.shape = (nlevel, nz, ny, nx)
        totn = np.fromfile(f, dtype = 'float32', count = nx*ny*nz)
        totn.shape = (nz, ny, nx)

        b = b.astype('float64')
        nstar = 10**nstar.astype('float64')
        totn = 10**totn.astype('float64')
        n = nstar*b

        pop = {'n': n, 'nstar': nstar, 'totn': totn, 'b': b}
        return pop

"""filename = PathToOutput+"/pop_ti_marcs_sun"
dep = readpop(filename, out['nx'], out['ny'], out['nz'], atom['nlevel'])
print(dep.keys())

print(len(dep['n'][0]))
print(len(dep['n']))
print((out['nx']))
print((out['ny']))
print((out['nz']))"""

