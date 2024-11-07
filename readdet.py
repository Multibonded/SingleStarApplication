import struct
import numpy as np
import const
from readsp import readsp

def readdet(filename, sp):
    """Read in the detailed output parameters.

    Parameters
    ----------
    filename : str
        Path to binary file to be read in.
    sp : dict
        The statistical equilibrium setup dictionary produced with readsp.

    Returns
    -------
    out : dict
        The detailed output parameters.
    """

    with open(filename, 'rb') as f:
        # read frequency indices
        nnu, = struct.unpack('i', f.read(4)) # struct.unpack returns a tuple even for a single element, we unpack that with the ,
        off = np.fromfile(f, dtype = 'int32', count = nnu)
        if len(sp) > 0:
            print("off", off, len(sp['wl']))
            wl = sp['wl'][off-1]
        else:
            wl = None
        freq = {'nnu':nnu, 'off':off, 'wl':wl}

        # read lines
        nline, = struct.unpack('i', f.read(4))
        if nline > 0:
            iatoml, wl0 = np.zeros((2, nline), dtype = 'int32')
            for iline in list(range(0, nline)):
                jatoml, = struct.unpack('i', f.read(4))
                iatoml[iline] = jatoml
                lam, = struct.unpack('d', f.read(8))
                wl0[iline] = lam*1e-2 # m
        else:
            iatoml, wl0 = None, None

        # read continua
        ncont, = struct.unpack('i', f.read(4))
        if ncont > 0:
            iatomc, nu0 = np.zeros((2, ncont), dtype = 'int32')
            for icont in list(range(0, ncont)):
                jatomc, = struct.unpack('i', f.read(4))
                iatomc[icont] = jatomc
                lam, = struct.unpack('d', f.read(8))
                nu0[icont] = lam/10**12 # s^-1
        else:
            iatomc, nu0 = None, None

        # read angles
        nangle, = struct.unpack('i', f.read(4))
        if nangle > 0:
            mus = np.fromfile(f, dtype = 'float64', count = 4*nangle)
            mux, muy, muz, wmu = np.reshape(mus, (4, nangle))
        else:
            mux, muy, muz, wmu = None, None, None, None

        nmuav, = struct.unpack('i', f.read(4))
        if nmuav > 0:
            muavs = np.fromfile(f, dtype = 'float64', count = 2*nmuav)
            muav, wmuav = np.reshape(muavs, (2, nmuav))
            nphiav = np.fromfile(f, dtype = 'int32', count = nmuav)
        else:
            muav, wmuav, nphiav = None, None, None

        # read geo
        nx, ny, nz, ngeo = struct.unpack('4i', f.read(4*4))
        if ngeo > 0:
            geos = np.fromfile(f, dtype = 'float64', count = 4*ngeo)
            geo_mux, geo_muy, geo_muz, geo_wmu = np.reshape(geos, (4, ngeo))
        else:
            geo_mux, geo_muy, geo_muz, geo_wmu = None, None, None, None

        out = {'nnu':nnu, 'off':off, 'wl':wl, 'nline':nline, 'iatoml':iatoml,
               'wl0':wl0, 'ncont':ncont, 'iatomc':iatomc, 'nu0':nu0,
               'nangle':nangle, 'mux':mux, 'muy':muy, 'muz':muz, 'wmu':wmu,
               'nmuav':nmuav, 'muav':muav, 'wmuav':wmuav, 'nphiav':nphiav,
               'nx':nx, 'ny':ny, 'nz':nz, 'ngeo':ngeo, 'geo_mux':geo_mux,
               'geo_muy':geo_muy, 'geo_muz':geo_muz, 'geo_wmu':geo_wmu}

    return out












