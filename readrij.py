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
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Sun\TiRun_Nq140k_k0001_490Abund\output\\"
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

