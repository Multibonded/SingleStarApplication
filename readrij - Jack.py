import numpy as np



import struct
import numpy as np
import matplotlib.pyplot as plt
import const
import sys
from readsp import readsp
from readdet import readdet
from readieJack import readie
from os import listdir
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Radiative_Run\output\\"
#PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Sun\TiRun_Nq140k_k0001_490Abund\output\\"
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

sp = readsp(filename_sp)
out = readdet(filename_det, sp)

def readrij(filename, out = None):
    """Read the radiative brackets.

    Parameters
    ----------
    filename : str
        Path to binary file. 

    Returns
    -------
    rij : 4d array
        The radiative brackets. 
    """

    with open(filename, 'r') as f:
        rij = np.fromfile(f, dtype = 'float32', count = out['nx']*out['ny']*out['nz']*(out['nline']+out['ncont']))
        rij.shape = (out['nline']+out['ncont'], out['nz'], out['ny'], out['nx'])
        rij = np.float64(rij)
        return rij


filename = PathToOutput + "rij_marcs_sun"
with open(filename, 'r') as f:
    rij = np.fromfile(f, dtype = 'float32', count = out['nx']*out['ny']*out['nz']*(out['nline']+out['ncont']))
    rij.shape = (out['nline']+out['ncont'], out['nz'], out['ny'], out['nx'])
    print(rij.shape          )

    rij = np.float64(rij)







d = {}
print(rij.shape)
print(rij)
import pickle




atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"
from astropy.io import fits

# Loading the BB and BF transitions to get the indexes and rates
atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data
k16 = "K16"












def findlowest():
    for x in range(len(rij)):
        print(max(rij[x]))
        d[x] = max(rij[x])


    # Get value for 90% of length of dict
    bottom90 = int(0.9*len(d.keys()))
    print(bottom90)
    print(int(bottom90))
    for x in range(bottom90):

        # index of lowest value
        mini = (min(d, key=d.get))
        # Remove lowest value
        d.pop(mini)
    print(len(d))
    print(d.keys())
    pickle.dump(d, open("radiative brackets/Old_lowestrad", "wb"))

print(rij)
print(rij.shape)