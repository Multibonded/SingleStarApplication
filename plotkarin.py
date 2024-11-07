from readsp import readsp
from readdet import readdet
from readie import readie
import specutils
import numpy as np
from specutils import analysis
import pickle

Mode = "avrg_xy"


# needs to be in AA, not nm
def Vac2Air(wavelength):
    # Vac2Air. Changes wavelengths.

    vac = wavelength
    s = 10 ** 4 / vac
    res = vac / (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2))
    air_wave = res

    return air_wave


abundances = ["445", "455", "460", "470", "480", "490", "495", "502", "520", "545"]
abundances = ["485"]
# saving all abundances to run in a different program
abund_dict = {}
karin = pickle.load(open("EquivWidths_Abundances_Karin.pkl", "rb"))


# run for each abundance folder we ran through balder.
for abundance in abundances:

    # new list for equiv widths for each abundance
    ew_list = []
    folder = "TiRun_Nq140k_k0005_" + abundance + "Abund"
    # folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"
    # folder = "TiRun_Nq140k_k0005_TenthQmaxAtom2"

    output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + folder + "\\output\\"
    output = r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\TiRun_Nq140k_k0005_485Abund_KarinsMods\\output\\"
    sp = readsp(output + 'sp2')
    out = readdet(output + 'detailed_marcs_sun', sp=sp)
    iet = readie(output + 'iet_marcs_sun', out=out, mode=Mode)
    iec = readie(output + 'iec_marcs_sun', out=out, mode=Mode)

    # Changing wl from vac to air
    out['wl'] = out['wl'] * 1E10
    air_wave = Vac2Air(out['wl'])

    # Lines from scott. 8000 isn't well produced for us... Maybe due to lack of lab data.
    lines = [4281.363, 4465.805, 4758.118, 5022.866, 5490.147, 6092.789, 7357.726, 8426.504]
    # After 7357 is ti II lines
    lines = [4281.363, 4465.805, 4758.118, 5022.866, 5490.147, 6092.789, 7357.726, 4409.52, 4657.212, 4719.533]
    lines = [4281.363, 4465.805, 4758.118, 5022.866, 5490.147, 6092.789, 7357.726]

    # lines = [2325.77, 8426.504]
    lower_val = [x - 5.3 for x in lines]
    lower_val = [x - 0.1 for x in lines]
    upper_val = [x + 5.3 for x in lines]
    upper_val = [x + 0.1 for x in lines]

    for lower_val, upper_val in zip(lower_val, upper_val):
        print("Working on", lower_val, "-", upper_val)
        # Limits the wavelengths to our selection (vacuum). L stands for Limit, so limx is limited x value
        lindex = np.where(np.logical_and(air_wave < upper_val, air_wave > lower_val))
        limx = air_wave[lindex]
        limy = (iet[-1, :] / iec[-1, :])[lindex]

        from specutils import Spectrum1D
        import matplotlib.pyplot as plt
        import astropy.units as u

        """vac_lower = 4281
        vac_upper = 4285
        # Now we have the AirWavelength, we want to constrain the wavelengths further as we can be more accurate now.
        lindex = np.where(np.logical_and(limx < vac_upper, limx > vac_lower))
        limx = limx[lindex]
        limy = limy[lindex]


        limx = Vac2Air(limx)
        limy = Vac2Air(limy)"""

        spec1d = Spectrum1D(spectral_axis=limx * u.AA, flux=limy * u.Jy)
        ew = analysis.equivalent_width(spec1d)
        print(ew)
        ew_list.append(ew.value)
        plt.figtext(0.2, 0.14, "EW: " + str(np.round(ew, 8)))

        plt.title("karin stuff")
        plt.plot(limx, limy,marker='o' )
        plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\spectra\485Abund_Karin\line"+str(int(np.round(lower_val+0.1)))+".png")

        plt.show()

    abund_dict[abundance] = ew_list


pickle.dump(abund_dict, open("EquivWidths_Abundances_Karin.pkl", "wb"))