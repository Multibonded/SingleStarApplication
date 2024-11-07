import struct
import numpy as np
import matplotlib.pyplot as plt
import const
import sys
from readsp import readsp
from readdet import readdet
from readie import readie
from os import listdir
from PyAstronomy import pyasl as pa

#outputfile =str(sys.argv[1])
Mode = "avrg_xy"


k = 0.005
k_string = str(k).replace(".", "")
print(k_string)
atom2check = True

Nq = 140
folder_name = "TiRun_Nq"+str(Nq)+"k_k"+k_string
#folder_name = "TiRun_Nq"+str(Nq)+"k_k"+k_string+"_545Abund"




PathNoOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+str(folder_name)+""
PathToOutput = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+str(folder_name)+"\output"
PathToAtom = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\atoms"
# removing the output from the string. Not atom2 as indexing gets messed up, and fine structures have identical DCs

files = listdir(PathToOutput)
for f in files:
    print("f", f)
    if "stagger" in f:
        hydro = "stagger"
    elif "marcs" in f:
        hydro = "marcs"
# needs to be in AA, not nm
def Vac2Air(wavelength):
    # Vac2Air. Changes wavelengths.

    vac = wavelength
    s = 10 ** 4 / vac
    res = vac / (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2))
    air_wave = res





    return air_wave

#############################################################
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
iec = readie(file_iec, out, mode = Mode)
iecl = readie(file_iecl, out, mode = Mode)
iet = readie(file_iet, out, mode = Mode)
ietl = readie(file_ietl, out, mode = Mode)

length = len(range(0,out['nnu']))
fluxc = np.array([0]*length, dtype=np.float64)
fluxcl = np.array([0]*length, dtype=np.float64)
fluxt = np.array([0]*length, dtype=np.float64)
fluxtl = np.array([0]*length, dtype=np.float64)

wl = out['wl']*10**10

wl2 = np.sort(wl)
print(wl2)
intherightrange = np.where(np.logical_and(wl2>612, wl2<613))
print(wl[intherightrange])
print(wl2[intherightrange])


for i in range(0,out['nnu']):
    fluxc[i] = sum(iec[:,i]*out['muz']*out['wmu'])
    fluxt[i] = sum(iet[:,i]*out['muz']*out['wmu'])

"""    # useless?
    fluxcl[i] = sum(iecl[:,i]*out['muz']*out['wmu'])
    # useless?
    fluxtl[i] = sum(ietl[:,i]*out['muz']*out['wmu'])"""





plt.figure()


#plt.plot(np.sort(wl),fluxt/fluxc)
plt.plot(Vac2Air(wl),iet[-1,:]/iec[-1,:])

plt.title("Spectra for "+folder_name)

lim = ""

lim = 6120
# ***   4641:  526.59639309 nm
#     150        22         0.036191560328006744       240   48.0 0   1   229086592.0   1.5 2.8840318009315524e-06
#    15220.393              11.0        'a3G:      3d3.(2G).4s'                 1    22     0
#     34204.985              9.0         'v3F:      3d2.(1G).4s.4p.(3P)'         1   150     0
#***   4641:  526.59639309 nm
    #150        22         0.036191560328006744       240   48.0 0   1   229086592.0   1.5 2.8840318009315524e-06


#atom1
#     14947.571708994708     27.0        'a3G:      3d3.(2G).4s'                 1     8      0
#     33886.93214285714      21.0        'v3F:      3d2.(1G).4s.4p.(3P)'         1    56      0
#***   10069:  528.0009340822619 nm
    #56         8          0.03229428827762604        27    27 0   1   229432896.0   1.5   0.0
plt.xlim([5265.5,5266.5])


plt.ylabel("Flux")
plt.xlabel("wavelength (AA)")
#-6 removes the output part of the name.
print(PathNoOutput+folder_name+"_spectra.png")
plt.savefig(PathNoOutput+folder_name+"_spectra.png")
plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\graphs\nq"+str(Nq)+"k\spectra\\"+folder_name+"_spectra_"+str(lim)+".png")

plt.show()
