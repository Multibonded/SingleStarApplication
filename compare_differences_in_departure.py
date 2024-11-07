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
import os


def compare():
    depart_dict = {}
    # We wanna look at each folder with the output at different k values
    for subdir, dirs, files in os.walk(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs"):

        for dir in dirs:
            if "80k" not in dir or "Checkup" in dir:
                continue
            name = dir.split("_")
            kvalue = name[-1]
            print(kvalue)
            depart_dict[kvalue] = {}

            # The k value we used to weight the lower equivalent widths when dividing up frequncy points
            k_string = kvalue[1:]
            if "." in k_string:
                k_string = str(k_string).replace(".", "")
                print(k_string)
            PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+dir+"\output"
            # removing the output
            atom=ratom(PathToOutput[:-6]+'/atom.ti_Nq80k_k'+k_string)

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

            filename = PathToOutput+"/pop_ti_marcs_sun"
            dep = readpop(filename, out['nx'], out['ny'], out['nz'], atom['nlevel'])
            print(dep.keys())

            print(len(dep['b'][0]))
            print(len(dep['b']))
            print((out['nx']))
            print((out['ny']))
            print((out['nz']))

            # The indexes of the same levels that berggeman tested: They represent a3F E
            # 0
            # a3H E
            # 12
            # t3F O
            # 80
            # x3H O
            # 86
            # a4F E
            # 459
            # a2P E
            # 465
            # z4D O
            # 479
            indexes = [0, 12, 80, 86, 459, 465, 479]

            xpoints = []
            ypoints = []
            ypoints_total = []
            for level in indexes:
                for depth in range(out['nz']):

                    xpoints.append(depth)
                    ypoints.append(dep['b'][level][depth][0][0])
                    ypoints_total.append(dep['b'][level][depth][0][0])
                if level >= 459:
                    plt.plot(xpoints, ypoints, "--", label= atom1['term'][level])
                else:
                    plt.plot(xpoints, ypoints,  label= atom1['term'][level])
                depart_dict[kvalue][level] = ypoints
                plt.legend()
                xpoints = []
                ypoints = []

            plt.title("Departure coefficients of levels using nq=80k and k="+ str(k_string))
            plt.xlabel("Geometric depth")
            plt.ylabel("Departure Coefficient")
            plt.yscale('log')
            plt.xscale('log')

            plt.figtext(0.2,0.14,"Total departure coefficients are: "+str(sum(ypoints_total)))
            #plt.show()
    pickle.dump(depart_dict, open("depart_dict", "wb"))

import pickle
verbose = False
#
differences = []
def examine():
    indexes = [0, 12, 80, 86, 459, 465, 479]

    dd = pickle.load(open("depart_dict", "rb"))

    # Temp list to check if a departure coeff ever changes (hint, yes)

    for index in indexes:
        departure_coeffs = []
        kvals = []
        for kval in dd:

            # For each k value, check the first index value ([0]) at the first geometric depth
            kvals.append(kval)
            if dd[kval][index][0] in departure_coeffs:
                print("EXISTS!", departure_coeffs)
                exit()
            else:
                departure_coeffs.append(dd[kval][index][0])


        # We want to rearrange them to make sure highest k is first (more nq for important levels)

        sorted_kvals, sorted_departure_coeffs = (list(t) for t in zip(*sorted(zip(kvals, departure_coeffs), reverse=True)))
        print("In decending order of k, the difference in departure coeff. from highest k is:", "\n")
        if verbose:
            for k_index in range(len(sorted_kvals)):
                print(sorted_kvals[0], "->", sorted_kvals[k_index])
                print(sorted_departure_coeffs[0], "->", sorted_departure_coeffs[k_index])
                print("Difference:", sorted_departure_coeffs[k_index]-sorted_departure_coeffs[0])
                print("----\n")
        else:
            temp_diff = []
            for k_index in range(len(sorted_kvals)):

                temp_diff.append(sorted_departure_coeffs[k_index] - sorted_departure_coeffs[0])
            differences.append(min(differences))
        exit()

examine()