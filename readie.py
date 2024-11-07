import numpy as np
import struct
import numpy as np
import matplotlib.pyplot as plt
import const
import sys
from readsp import readsp
from readdet import readdet
from os import listdir
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\TiRun_Nq80k_k003_full\output"
files = listdir(PathToOutput)
for f in files:
	if "stagger" in f:
		hydro = "stagger"
	elif "marcs" in f:
		hydro = "marcs"
filename_det = PathToOutput + "/detailed_"+hydro+"_sun"
atom2check = True
if atom2check:
    filename_sp = PathToOutput +"/sp2"
else:
    filename_sp = PathToOutput + "/]#[]sp"

file_iec = PathToOutput +"/iec_"+hydro+"_sun"
file_iecl = PathToOutput +"/iecl_"+hydro+"_sun"
file_iet = PathToOutput +"/iet_"+hydro+"_sun"
file_ietl = PathToOutput +"/ietl_"+hydro+"_sun"

sp = readsp(filename_sp)
out = readdet(filename_det, sp)

def readie(filename, out, mode = 'avrg_xy'):
    """Read in the intensities.

    Parameters
    ----------
    filename : str
        Path to binary file.
    out : dict
        The detailed output parameters dictionary produced with readdet.
    mode : str
        The type of ie file. Possible modes: 'resolved', 'avrg_xy', 'avrg_xy_phi'.

    Returns
    -------
    ie : ndarray
        The line spectrum across all angles.
    """

    with open(filename, 'rb') as f:
        ie = np.fromfile(f, dtype = 'float64')
        if mode == 'resolved':
            ie.shape = (out['nnu'], out['nangle'], out['ny'], out['nx'])
        elif mode == 'avrg_xy':
            ie.shape = (out['nnu'], out['nangle'])
        elif mode == 'avrg_xy_phi':
            ie.shape = (out['nnu'], out['nmuav'])
        else:
            raise ValueError('No match found for mode.')
    return ie.T
filename = PathToOutput+r"\ietl_marcs_sun"

ie = readie(filename, out, mode = 'avrg_xy')
print(ie.shape)
