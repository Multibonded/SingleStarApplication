import pickle
import numpy as np
import  matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"
from os import listdir
from readsp import readsp
from readdet import readdet
from ratom import ratom

atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data
k16 = "K16"
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Radiative_Run\output\\"
read_atom=ratom(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\atoms"+'/atom.ti_Nq140k_k0001')


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



indexes = pickle.load(open("radiative brackets/atom1_indexes", "rb"))
print(indexes)
rad = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\PYtools\radiative brackets\Radiative_RunAtom1.pkl", "rb"))
print(rad[0].shape)
print(len(atom2))











radiative = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\Final_outputs\full_radiative_xsection_nahar+regemorter.pkl", "rb"))

initial = []
final = []
bf_wl = []
for bf_rad_transition in range(len(radiative)):
    # Loading interpolated cross sections, in order of original atom1 transitions. Each file is a different level to level transition, consisting of many different wavelengths on a set grid.
    xsections = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\interpolation_wavelengths\interpolated_cross_sections/cross_section_interpolated_transition"+str(bf_rad_transition)+".pkl", "rb"))
    # Loading the wavelengths that gave rise to the interpolated cross sections
    wavelengths = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\interpolation_wavelengths\interpolated_grids\wavelength_grid_interpolated_transition"+str(bf_rad_transition)+".pkl", "rb"))
    # Is still in the same order as the original radiative, so all other information is correct (Target, e.g)
    # Priting transition + 1 as we don't start at 0 in balder.
    # Mistake iwth radiative code, should be going to the 2nd state in atom2 from singlets.
    if int(radiative[bf_rad_transition]['target']) == 460:
        (radiative[bf_rad_transition]['target']) = 461
    initial.append(bf_rad_transition+1)
    final.append(int(radiative[bf_rad_transition]['target'])+1)
    bf_wl.append(radiative[bf_rad_transition]['wavelength'][0])
    # So for each photon wavelength/energy, we have an individual cross section.


from readpop import readpop

pop = readpop(PathToOutput+"/pop_ti_marcs_sun", out['nx'], out['ny'], out['nz'], read_atom['nlevel'])
print(pop.keys())
print(len(pop['n']))
n = pop['n']



ti_i = []
ti_f = []
ti1_indexes = []
ti2_indexes = []
initial_ec = []
y = []
wl = []
for trans in range(len(rad)):
    rate = rad[trans][0][0][0]
    if trans < 4784:
        radini = (atom2[indexes[trans]]['irad'])
        radfin = atom2[indexes[trans]]['jrad']
        wl.append(atom2[indexes[trans]]['alamb'])

    else:
        radini = initial[trans-4784]
        radfin = final[trans-4784]
        wl.append(bf_wl[trans-4784])
    ediff = atom1[radfin-1]['ec'] - atom1[radini-1]['ec']
    if radini <= 459:
        ti1_indexes.append(trans)
    else:
        ti2_indexes.append(trans)

    initial_ec.append(atom1[radini-1]['ec'])
    ti_i.append(radini-1)
    ti_f.append(radfin-1)
    y.append(rate)

nlist = []
# Only the population for thje first depth of the atmosphere (first [0])
for x in range(len(n)):
    nlist.append(n[x][0][0][0])
nlist = np.asarray(nlist)
# Number population of all initial pops for each transition to be applied to rij directly (as it represnts all transitions)
ni = nlist[ti_i]
nj = nlist[ti_f]


filename = PathToOutput + "rij_marcs_sun"
filename2 = PathToOutput + "rji_marcs_sun"

with open(filename, 'r') as f:
    rij = np.fromfile(f, dtype = 'float32', count = out['nx']*out['ny']*out['nz']*(out['nline']+out['ncont']))
    rij.shape = (out['nline']+out['ncont'], out['nz'], out['ny'], out['nx'])
    #rij = np.float64(rij)



filename2 = PathToOutput + "rji_marcs_sun"
with open(filename2, 'r') as f2:
    rji = np.fromfile(f2, dtype = 'float32', count = out['nx']*out['ny']*out['nz']*(out['nline']+out['ncont']))
    rji.shape = (out['nline']+out['ncont'], out['nz'], out['ny'], out['nx'])

    #rji = np.float64(rji)

rad_bracket = []
for r in range(len(rij)):
    bracket_ij = ni[r]*rij[r][20][0][0] - nj[r]*rji[r][20][0][0]
    print(bracket_ij)
    rad_bracket.append(bracket_ij)
print(rad_bracket)

initial_ec = np.asarray(initial_ec)
rad_bracket = np.asarray(rad_bracket)
print(initial_ec.shape)
print(rad_bracket.shape)
plt.scatter(initial_ec[ti1_indexes],rad_bracket[ti1_indexes],color= "blue", s=2 )
plt.scatter(initial_ec[ti2_indexes],rad_bracket[ti2_indexes], color="red", s=2 )



dictionary = {"initial_ec": initial_ec, "rad_bracket":rad_bracket, "ti1_indexes":ti1_indexes, "ti2_indexes":ti2_indexes, "wl": wl}

#plt.xscale('log')
plt.yscale('symlog')


pickle.dump(dictionary, open("radiative brackets/radiative_dict", "wb"))
plt.show()