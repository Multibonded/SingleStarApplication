import numpy as np
import struct
import numpy as np
import matplotlib.pyplot as plt
import const
import sys
from readsp import readsp
from readdet import readdet
from os import listdir

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

"""output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Radiative_Run\output\\"
sp = readsp(output + 'sp')

out = readdet(output + 'detailed_marcs_sun', sp=sp)
filename = output+'change_ie'
with open(filename, 'r') as f:
    rij = np.fromfile(f, dtype='float32')
    rij.shape = (out['nline'] + out['ncont'], out['nz'], out['ny'], out['nx'])

print(len(rij))
"""