
from readsp import readsp
from readdet import readdet
from readie import readie
import specutils
import numpy as np
from specutils import analysis
Mode = "avrg_xy"
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

# this is old vturb=1.5 results
abundances = ["445", "455", "460", "470", "480", "485", "490", "495", "502", "520", "545"]
#for vturb=1
abundances=["445", "455", "460",  "465", "475", "485", "495", "502"]
#abundances = ["485"]
# saving all abundances to run in a different program
abund_dict = {}
ti2 = False
lte = True


# needs to be in AA, not nm
def Vac2Air(wavelength):
    # Vac2Air. Changes wavelengths.

    vac = wavelength
    s = 10 ** 4 / vac
    res = vac / (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2))
    air_wave = res





    return air_wave


# run for each abundance folder we ran through balder.
for abundance in abundances:

    # new list for equiv widths for each abundance
    ew_list = []
    folder = "TiRun_Nq140k_k0005_"+abundance+"Abund"
    #folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
    #folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"

    output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + folder + "\\output\\"

    sp = readsp(output + 'sp2')
    out = readdet(output + 'detailed_marcs_sun', sp=sp)
    if lte:
        iet = readie(output + 'ietl_marcs_sun', out=out, mode=Mode)
        iec = readie(output + 'iecl_marcs_sun', out=out, mode=Mode)
    else:
        iet = readie(output + 'iet_marcs_sun', out=out, mode=Mode)
        iec = readie(output + 'iec_marcs_sun', out=out, mode=Mode)

    # Changing wl from vac to air
    out['wl'] = out['wl'] * 1E10
    air_wave = Vac2Air(out['wl'])

    loop_count = 0
    plt.plot(air_wave, iet[-1,:]/iec[-1,:])
    plt.xlim(-100, 30000)
    plt.show()