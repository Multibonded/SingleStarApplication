# Plot departure coefficients of a run, either atom2 or 1 depending on the Path chosen and filenames.

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
# The k value we used to weight the lower equivalent widths when dividing up frequncy points
atom2check = True
Nq = 140


optical_depth_allinfo = """k    log(tau)  T    log(Pe)   log(Pg)    rhox
1  -4.9153  4063.00  -1.6770   2.4269    0.545404E-02
2  -4.7893  4084.51  -1.6066   2.4988    0.719670E-02
3  -4.6654  4110.15  -1.5356   2.5701    0.924411E-02
4  -4.5431  4138.08  -1.4643   2.6409    0.116486E-01
5  -4.4219  4166.44  -1.3930   2.7116    0.144528E-01
6  -4.3013  4194.94  -1.3218   2.7822    0.177448E-01
7  -4.1810  4223.68  -1.2508   2.8526    0.216080E-01
8  -4.0608  4252.53  -1.1801   2.9229    0.261400E-01
9  -3.9404  4281.39  -1.1095   2.9931    0.314558E-01
10  -3.8199  4310.21  -1.0391   3.0630    0.376903E-01
11  -3.6992  4339.00  -0.9691   3.1328    0.450010E-01
12  -3.5782  4367.73  -0.8992   3.2025    0.535717E-01
13  -3.4569  4396.37  -0.8295   3.2722    0.636166E-01
14  -3.3353  4424.99  -0.7600   3.3417    0.753873E-01
15  -3.2135  4453.43  -0.6909   3.4110    0.891824E-01
16  -3.0913  4481.34  -0.6219   3.4804    0.105350E+00
17  -2.9689  4509.58  -0.5531   3.5496    0.124287E+00
18  -2.8464  4537.99  -0.4843   3.6187    0.146459E+00
19  -2.7236  4566.07  -0.4159   3.6878    0.172421E+00
20  -2.6007  4594.10  -0.3477   3.7568    0.202818E+00
21  -2.4775  4622.20  -0.2795   3.8257    0.238409E+00
22  -2.3544  4650.55  -0.2114   3.8945    0.280061E+00
23  -2.2311  4679.40  -0.1433   3.9632    0.328812E+00
24  -2.1078  4708.95  -0.0752   4.0318    0.385859E+00
25  -1.9844  4739.82  -0.0067   4.1007    0.452602E+00
26  -1.8610  4772.46   0.0621   4.1690    0.530664E+00
27  -1.7377  4807.62   0.1319   4.2375    0.621944E+00
28  -1.6145  4846.15   0.2024   4.3058    0.728621E+00
29  -1.4913  4888.24   0.2743   4.3741    0.853252E+00
30  -1.3684  4937.62   0.3483   4.4420    0.998741E+00
31  -1.2456  4991.77   0.4240   4.5099    0.116841E+01
32  -1.1231  5054.81   0.5027   4.5776    0.136610E+01
33  -1.0006  5127.63   0.5853   4.6450    0.159609E+01
34  -0.8784  5213.00   0.6733   4.7119    0.186270E+01
35  -0.7563  5313.33   0.7690   4.7782    0.216983E+01
36  -0.6341  5431.80   0.8763   4.8429    0.251960E+01
37  -0.5117  5571.29   1.0003   4.9055    0.290989E+01
38  -0.3889  5736.02   1.1478   4.9641    0.333083E+01
39  -0.2659  5930.44   1.3251   5.0172    0.376340E+01
40  -0.1425  6158.57   1.5341   5.0631    0.418257E+01
41  -0.0188  6426.71   1.7745   5.1014    0.456504E+01
42   0.1047  6738.81   2.0410   5.1320    0.489469E+01
43   0.2271  7127.58   2.3464   5.1555    0.516185E+01
44   0.3473  7533.54   2.6373   5.1725    0.536484E+01
45   0.4653  7870.54   2.8592   5.1862    0.552910E+01
46   0.5817  8146.43   3.0293   5.1982    0.568057E+01
47   0.6976  8380.52   3.1662   5.2096    0.583030E+01
48   0.8136  8586.84   3.2819   5.2210    0.598484E+01
49   0.9301  8773.23   3.3829   5.2327    0.614916E+01
50   1.0474  8947.45   3.4745   5.2452    0.632776E+01
51   1.1655  9116.06   3.5606   5.2586    0.652442E+01
52   1.2845  9280.76   3.6426   5.2729    0.674264E+01
53   1.4045  9443.13   3.7214   5.2882    0.698650E+01
54   1.5256  9604.79   3.7983   5.3049    0.726045E+01
55   1.6477  9766.94   3.8736   5.3230    0.756949E+01
56   1.7709  9930.74   3.9481   5.3425    0.792324E+01"""

optical_depth = []
for x in optical_depth_allinfo.split("\n"):
    if "log" in x.split()[1]:
        continue
    optical_depth.append(float(x.split()[1]))

print(optical_depth)
star_number = 0
star_list = ['sun', 'Arcturus', '122563', '140283', '84937']
star = star_list[star_number]
k = "0001"
if k != "0005":
    pass
else:
    k=""
k_string = "0001"

# using our abundance calc
if star_number == 3:
    folder_name = "TiRun_Nq140k_k0001_270Abund"
elif star_number == 0:
    folder_name = "TiRun_Nq140k_k0001_490Abund"
elif star_number == 4:
    folder_name = "TiRun_Nq140k_k0001_330Abund"
elif star_number == 2:
    folder_name = "TiRun_Nq140k_k0001_230Abund"
elif star_number == 1:
    folder_name ="TiRun_Nq140k_k0001_470Abund"
#folder_name = "TiRun_Nq"+str(Nq)+"k_k"+k_string+"_full_halved"

PathNoOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\"+str(folder_name)+"\\"
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\"+str(folder_name)+"\output\\"
PathToAtom = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\atoms\\"
# removing the output from the string. Not atom2 as indexing gets messed up, and fine structures have identical DCs
atom=ratom(PathToAtom+'atom.ti_Nq'+str(Nq)+'k_k'+k_string)

files = listdir(PathToOutput)
for f in files:
	if "stagger" in f:
		hydro = "stagger"
	elif "marcs" in f:
		hydro = "marcs"
filename_det = PathToOutput + "/detailed_"+hydro+"_sun"
if atom2check:
    filename_sp = PathToOutput +"/sp2"
else:
    filename_sp = PathToOutput + "/sp"

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


# The indexes of the same levels that berggeman tested:
atom111 = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
atom22 = r"D:\PhD\TiFeAtoms\Jack\balder\atom2\atom.Ti782.fits"
atom111 = fits.open(atom111)
atom22 = fits.open(atom22)
atom111 = atom111[1].data
atom22 = atom22[1].data

# terms wanted: z4Do, a2P, a4F, x3Ho, t3Fo, a3H, a3F
indexes = [0, 12, 80, 86, 459, 465, 479]

indexes2 = []
# The indexes are for the atom1, but we need the departure coefficients
print("here start", atom1['term'][indexes])

for i in range(0, len(indexes)):
    atom2i = np.where(np.logical_and(atom22['term'] == atom1['term'][indexes[i]], atom22['par'] == atom1['par'][indexes[i]]))
    print(atom2i, atom22['par'][atom2i])
    print(atom22['term'][atom2i])
    indexes2.append(atom2i[0][0])
indexes = indexes2
print("indexes", indexes)
print("terms used", atom22['term'][indexes])




if atom2check:
    filename = PathToOutput+"/pop2_ti_marcs_sun"
else:
    filename = PathToOutput+"/pop_ti_marcs_sun"


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
dep = readpop(filename, out['nx'], out['ny'], out['nz'], atom['nlevel'])
print(dep.keys())

print(len(dep['b'][0]))
print(len(dep['b']))
print((out['nx']))
print((out['ny']))
print((out['nz']))
exit()




xpoints = []
ypoints = []
ypoints_total = []
for level in indexes:
    for depth in range(out['nz']):

        xpoints.append(depth)
        ypoints.append(np.log10(dep['b'][level][depth][0][0]))
        ypoints_total.append(dep['b'][level][depth][0][0])
    if level >= 459:
        plt.plot(xpoints, ypoints, "--", label= atom22['term'][level])
    else:
        plt.plot(xpoints, ypoints,  label= atom22['term'][level])

    plt.legend()
    xpoints = []
    ypoints = []

#plt.title("Departure coefficients of levels using nq="+str(Nq)+"k and k="+ str(k))
plt.xlabel("Geometric depth")
plt.ylabel("$Log_{10}$(b)")
plt.xscale('log')
#plt.figtext(0.2,0.14,"DC sum: "+str(sum(ypoints_total)))

if atom2check:
    balder_run = "atom2TiRun_Nq"+str(Nq)+"k_k"+k_string
else:
    balder_run = "TiRun_Nq"+str(Nq)+"k_k"+k_string

print(PathNoOutput+balder_run+"_departures.png")
plt.savefig(PathNoOutput+balder_run+"_departures.png")

plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\\departures\\"+balder_run+"_departures.png")


plt.show()
