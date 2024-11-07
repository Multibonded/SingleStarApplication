from __future__ import division
import struct
import numpy as np
import const

def readsp(filename):
    """Read in the statistical equilibrium setup.

    Parameters
    ----------
    filename : str
       Path to binary file.

    Returns
    -------
    sp : dict
        The statistical equilibrium setup.
    """

    with open(filename, "rb") as f:
        nnu, maxac, maxal = struct.unpack('3i', f.read(4*3))
        nu = np.fromfile(f, dtype = 'float64', count = nnu)
        wnu = np.fromfile(f, dtype = 'float64', count = nnu)

        if maxac > 0:
            ac = np.fromfile(f, dtype = 'int32', count = maxac*nnu)
            ac.shape = (nnu, maxac)
            ac = ac.T
        else:
            ac = None

        if maxal > 0:
            al = np.fromfile(f, dtype = 'int32', count = maxal*nnu)
            al.shape = (nnu, maxal)
            al = al.T
        else:
            al = None

        nac = np.fromfile(f, dtype = 'int32', count = nnu)
        nal = np.fromfile(f, dtype = 'int32', count = nnu)

        wl = const.cc/nu

        sp = {'nnu':nnu, 'maxac':maxac, 'maxal':maxal, 'nu':nu, 'wnu':wnu,
              'ac':ac, 'al':al, 'nac':nac, 'nal':nal, 'wl':wl}

        return sp
