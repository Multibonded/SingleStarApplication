
from readsp import readsp
from readdet import readdet
from readieJack import readie
import specutils
import numpy as np
from specutils import analysis
from specutils import Spectrum1D
import matplotlib.pyplot as plt
import astropy.units as u
from scipy.interpolate import UnivariateSpline
import pickle
Mode = "avrg_xy"

abundances_string = ["265", "270", "275", "280", "290"]
star = "chris"
k = "0001"
ew_dict = {}
wavelength_dict = {}


def Vac2Air(wavelength):
    # Vac2Air. Changes wavelengths.

    vac = wavelength
    s = 10 ** 4 / vac
    res = vac / (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2))
    air_wave = res

    return air_wave

lte = False
for abundance in abundances_string:
    count = 0
    wavelength_dict[abundance] = []
    ew_dict[abundance] = []
    # new list for equiv widths for each abundance
    ew_list = []
    if k == "":
        folder = star + "\TiRun_Nq140k_k0005_" + abundance + "Abund"

    else:
        folder = star + "\TiRun_Nq140k_k" + k + "_" + abundance + "Abund"
    # folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
    # folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
    if star == "chris":
        folder = star + "\\Abund" + abundance
    output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + folder + "\\output\\"
    sp = readsp(output + 'sp2')
    out = readdet(output + 'detailed_marcs_sun', sp=sp)
    if lte:
        iet = readie(output + 'ietl_marcs_sun', out=out, mode=Mode)
        iec = readie(output + 'iecl_marcs_sun', out=out, mode=Mode)
    else:
        iet = readie(output + 'iet_marcs_sun', out=out, mode=Mode)
        iec = readie(output + 'iec_marcs_sun', out=out, mode=Mode)

    fluxlistt = []
    fluxlist2c = []

    iet2 = np.copy(iet)
    iec2 = np.copy(iec)
    # HAcky way to replace disk centre intensity with the flux. Literally replaces that column if not for the sun
    for x in range((out['nnu'])):
        fluxlistt.append(sum(iet[:, x] * out['muz'] * out['wmu']))
        fluxlist2c.append(sum(iec[:, x] * out['muz'] * out['wmu']))
    iet[-1] = np.asarray(fluxlistt)
    iec[-1] = np.asarray(fluxlist2c)
    print("c", sum(iec[-1]), "\n")
    print("t", sum(iet[-1]), "\n")
    print("ratio of flux", sum(iet2[-1]) / sum(iec2[-1]))
    print(out['muz'])
    # exit()
    # Changing wl from vac to air, but from nm to AA first
    out['wl'] = out['wl'] * 1E10
    air_wave = Vac2Air(out['wl'])

    limy = (iet[-1, :] / iec[-1, :])
    with open("Chris/ChrisAbundance"+abundance+".txt", "w") as text_file:
        text_file.write("Wavelength (AA)             Flux\n")
        for x in range(len(limy)):
            text_file.write(f"{np.round(out['wl'][x], 4):<28} {limy[x]:<25}")
            text_file.write("\n")
