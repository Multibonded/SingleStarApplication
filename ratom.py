from __future__ import division
import const
import numpy as np
from units import cmm1_to_ev

def nextline(openfile):
    """Reads a line. Discards the line if it begins with '*'. Picks up new line if discarded.

    Parameters
    ----------
    openfile : file
        Open file object.

    Returns
    -------
    line : str
        The next line of the file that is not a comment.
    """

    line = openfile.readline()
    while line[0] == '*':
        line = openfile.readline()
    return line

def ratom(filename):
    """Reads in the m3d model atom. Takes filename of atom without 'atom.'. Returns atom as dictionary.

    Parameters
    ----------
    filename : str
        Name of the atom file without 'atom.'.

    Returns
    -------
    atom : dict
        The atom details.
    """

    max_nq = 2000
    #filename = 'atom.' + filename

    with open(filename, 'r') as f:
        species = nextline(f).split()[0] # gets rid of the '\n'
        abund, awgt = map(float, nextline(f).split())
        nlevel, nline, ncont, nfix = map(int, nextline(f).split())

        #*   E [cm^-1]       g  label
        level = []
        for ilevel in range(nlevel):
            level_dict = {'ev': None, 'g': None, 'label': None} # enforce ordering of inputs
            E_g, label, ion = nextline(f).split("'")
            level_dict['label'] = label
            ev, g = map(float, E_g.split())
            level_dict['ev'], level_dict['g'] = cmm1_to_ev(ev), g
            ion = ion.split()
            j = None
            if len(ion) == 1:
                ion = ion[0]
            elif len(ion) == 2:
                ion, _ = ion
            elif len(ion) == 3:
                ion, _, j = ion
                j = float(j)
            level_dict['ion'] = int(ion)
            level_dict['j'] = j
            level.append(level_dict)

        #* upp  low          f   nq     qmax    q0 iw         Grad           GW           GQ       wl [A]
        line = []
        for iline in range(0, nline):
            line_dict = {'ul': None, 'll': None, 'f': None, 'nq': None,
                         'qmax': None, 'q0': None, 'ncomp': None, 'ga': None,
                         'gw': None, 'gs': None, 'wl': None}
            ul, ll, line_dict['f'], nq, line_dict['qmax'], line_dict['q0'], ncomp, line_dict['ga'], line_dict['gw'], line_dict['gs'], *_ = map(float, nextline(f).split()) # the last input is not always present
            ul, ll, nq, ncomp = int(ul), int(ll), int(nq), int(ncomp)
            line_dict['ul'] = ul
            line_dict['ll'] = ll
            line_dict['nq'] = nq
            line_dict['ncomp'] = ncomp
            line_dict['wl'] = const.hh*const.cc/(const.ee*(level[ul-1]['ev'] - level[ll-1]['ev']))
            #*      CEN     CWG          F             GA        GV         GQ
            if ncomp > 1:
                line_dict['denergy'] = np.zeros(ncomp)
                line_dict['weight'] = np.zeros(ncomp)
                for icomp in range(0, ncomp):
                    denergy, weight, *_ = map(float, nextline(f).split())
                    line_dict['denergy'][icomp] = denergy
                    line_dict['weight'][icomp] = weight
            line.append(line_dict)

        #* bound-free transitions
        if ncont > 0:
            cont = []
            for icont in range(0, ncont):
                cont_dict = {'ul': None, 'll': None, 'f': None, 'nq': None,
                             'qmax': None}
                ul, ll, cont_dict['f'], nq, qmax = map(float, nextline(f).split())
                ul, ll, nq = int(ul), int(ll), int(nq)
                cont_dict['ul'] = ul
                cont_dict['ll'] = ll
                cont_dict['nq'] = nq
                cont_dict['qmax'] = qmax
                if qmax < 0:
                    cont_dict['wavelen'] = np.zeros(nq, dtype = 'float64')
                    cont_dict['alpha'] = np.zeros(nq, dtype = 'float64')
                    for iq in range(0, nq):
                        wl, al = nextline(f).split()
                        cont_dict['wavelen'][iq] = wl
                        cont_dict['alpha'][iq] = al
                cont.append(cont_dict)
        else:
            cont = None

        colroutine = nextline(f).split()[0] # NEWCOL or GENCOL

        coldata = []
        row = nextline(f)
        while not row.split()[0] == 'END':
            coldata_dict = {'il': None, 'fl': None, 'type': None, 'nT': None,
                            'T': None, 'rates': None}
            row = row.split(maxsplit = 1) # if there are comments
            if row[0] == 'TEMP':
                temps = nextline(f).split()
                nT = int(temps[0])
                T = list(map(float, temps[1:]))
            else:
                coldata_dict['type'] = row[0]
                if len(row) > 1:
                    coldata_dict['comment'] = row[1][:-1] # the -1 gets rid of the '\n' in the comment
                il, fl, *rates = nextline(f).split()
                il, fl = int(il), int(fl)
                rates = list(map(float, rates))
                coldata_dict['il'] = il
                coldata_dict['fl'] = fl
                coldata_dict['nT'] = nT
                coldata_dict['T'] = T
                coldata_dict['rates'] = rates
                coldata.append(coldata_dict)
            row = nextline(f)

        ncoldata = len(coldata)

    # atom struct
    atom={'species':species, 'abund':abund, 'awgt':awgt,
          'nlevel':nlevel, 'nline':nline, 'ncont':ncont,
          'nfix':nfix, 'level':level, 'line':line, 'cont':cont,
          'colroutine':colroutine, 'ncoldata':ncoldata,
          'coldata':coldata}

    return atom
