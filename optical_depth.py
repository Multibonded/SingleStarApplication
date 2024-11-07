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
wavelength_dict = {}
ew_dict = {}
Mode = "avrg_xy"
rtc = 1.0973731569E5  # Conversion factor between Rydberg and cm^(-1)
rte = 13.605698  # Conversion factor between Rydberg and eV

# run for each abundance folder we ran through balder.
# new list for equiv widths for each abundance
ew_list = []
# folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
# folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\Radiative_Run\output\\"
sp = readsp(output + 'sp')
out = readdet(output + 'detailed_marcs_sun', sp=sp)
"""if lte:
    iet = readie(output + 'ietl_marcs_sun', out=out, mode=Mode)
    iec = readie(output + 'iecl_marcs_sun', out=out, mode=Mode)
else:"""
iet = readie(output + 'iet_marcs_sun', out=out, mode=Mode)
iec = readie(output + 'iec_marcs_sun', out=out, mode=Mode)
print(iet[-1].shape)
exit()
fluxlistt = []
fluxlist2c = []
import numpy as np

print("c", sum(iec[-1]), "\n")
print("t", sum(iet[-1]), "\n")
print("ratio of intensity", sum(iet[-1]) / sum(iec[-1]), "\n")

iet2 = np.copy(iet)
iec2 = np.copy(iec)
# HAcky way to replace disk centre intensity with the flux. Literally replaces that column if not for the sun
for x in range((out['nnu'])):
    fluxlistt.append(sum(iet[:, x] * out['muz'] * out['wmu']))
    fluxlist2c.append(sum(iec[:, x] * out['muz'] * out['wmu']))
iet[-1] = np.asarray(fluxlistt)
iec[-1] = np.asarray(fluxlist2c)
limy = (iet[-1, :] / iec[-1, :])
print(iet[-1])
