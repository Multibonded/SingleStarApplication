import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp1d
import numpy as np
from astropy.io import fits

atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"

atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data

# Which atom/lte version do we run? And for which atmosphere?
ti2= True
lte = False
# which stellar atmos are we looking at
star_number = 19
            #  0        1          2          3         4        5                 6               7     8(10% electron rates)     9    10                      11                  12              13
star_list = ['sun', 'Arcturus', '122563', '140283', '84937', "140283Karin", "arcturusKarin", "122563_3D", "122563_10th", "122563_100th", "122563_0kaulakys", "Arcturus_100th", "Arcturus_10th", "Arcturus_0electrons"
    , "Sun_0kaulakys100electrons", "Arcturus_0kaulakys100electrons", "140283_0kaulakys100electrons",     "84397_0kaulakys100electrons", "122563_0kaulakys100electrons", "chris"]
        # 14                             15                              16                                      17                          18
star = star_list[star_number]

#abundances = [["445", "455", "460", "470", "480", "485", "490", "495", "502", "520", "545"]]
# k0001
k = "0001"
# WE did more abundances for the sun at k=0.005
# abunds for sun, arcturus, 122563, 140283, and 84937
abundances_all = [["450", "460", "480", "490", "500"], ["460", "470", "490", "500"],
                  ["200", "220", "230", "250", "270"],
                  ["250", "270", "280", "300"], ["290", "310", "320", "330"],
                  ["250", "260", "270", "280", "290", "300"], ["460", "470", "490", "500"],
                  ["200", "220", "230", "240", "250", "260"], ["210", "230", "250", "270"],
                  ["210", "230", "250", "270"], ["210", "230", "250", "270"], ["460", "480", "500"],
                  ["460", "480", "500"], ["460", "480", "500"],
                  ["450", "460", "480", "490", "510"], ["460", "470", "490", "500"],
                  ["250", "270", "280", "300"], ["280", "290", "300", "320", "330"],
                  ["210", "230", "250", "270"], ["265", "270", "275", "280", "290"]]
abundances_string = abundances_all[star_number]


def Vac2Air(wavelength):
    # Vac2Air. Changes wavelengths.

    vac = wavelength
    s = 10 ** 4 / vac
    res = vac / (1 + 0.0000834254 + 0.02406147 / (130 - s ** 2) + 0.00015998 / (38.9 - s ** 2))
    air_wave = res

    return air_wave


# Turning them into floats.
abundances = [float(abund)/100 for abund in abundances_string]

# Dict in form d[abundance] = [x,y,z...] where x,y,z is Equiv Width for lines 4281 to 8426 as above
# Made in jack_spectra to get equiv widths of all lines in literature
if lte:
    if ti2:
        equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"Ti2_EquivWidths_AbundancesLTE.pkl", "rb"))
    else:
        equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"EquivWidths_AbundancesLTE.pkl", "rb"))
else:
    if ti2:
        equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"Ti2_EquivWidths_AbundancesNLTE.pkl", "rb"))
    else:
        equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"EquivWidths_AbundancesNLTE.pkl", "rb"))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def ti1get_chris_data():
    # don'#t forget it takes indentation as text space!
    scott_data = """Ti I    3354.63   0.02   0.11   16.1  -5.32   2.73
Ti I    3371.45   0.05   0.23   18.2  -5.27   2.70
Ti I    3635.46   0.00   0.10   16.8  -5.33   2.71
Ti I    3653.49   0.05   0.28   22.0  -5.22   2.73
Ti I    3741.06   0.02  -0.15   12.1  -5.49   2.72
Ti I    4533.24   0.85   0.54   12.5  -5.56   2.75
Ti I    4534.78   0.84   0.35   11.3  -5.60   2.88
Ti I    5210.38   0.05  -0.82    4.0  -6.12   2.77"""

    # removed Ti II   3348.84   0.12  -1.18   55.4  -4.78   2.79
    # Ti II   3349.03   0.61   0.46   86.1  -4.59   2.64
    # Ti II   3349.40   0.05   0.54  124.0  -4.43   2.81
    # Ti II   3477.18   0.12  -0.95   65.7  -4.72   2.86
    # Ti II   4394.06   1.22  -1.77    9.7  -5.66   2.89
    ewlist = []
    threeDlist = []
    wllist = []
    nlte3dlist = []
    abund_list = []
    ltelist = []
    energylist = []
    gflist = []
    for line in (scott_data.split("\n")):
        scott_split = line.split()
        print(scott_split)

        wl = float(scott_split[2])
        wllist.append(wl)
        rtc = 1.0973731569E5  # Conversion factor between Rydberg and cm^(-1)
        rte = 13.605698  # Conversion factor between Rydberg and eV

        energylist.append((4.135667E-15 * 299792458 / (wl*1E-10)))
        gflist.append(float(scott_split[4]))
        if ti2:
            ew = float(scott_split[-3])
        elif not ti2:
            ew = float(scott_split[-3])
        # To mAA from AA
        ew = np.log10(ew / 1000)
        ewlist.append(ew)

        nlte3d = float(scott_split[-1])
        nlte3dlist.append(nlte3d)
        abund_list.append(nlte3d)
        lte = float(scott_split[-1])
        ltelist.append(lte)

    extra_dict = {}
    extra_dict['ew'] = ewlist
    extra_dict['loggf'] = gflist
    extra_dict['energy'] = energylist
    return wllist, ewlist, abund_list, abund_list, extra_dict


def ti2get_chris_data():
        scott_data = """Ti II   3106.23   1.24  -0.07   52.2  -4.78   2.72
Ti II   3161.20   0.11  -0.74   65.9  -4.68   2.77
Ti II   3202.53   1.08   0.07   58.8  -4.74   2.61
Ti II   3222.84   0.01  -0.40   85.5  -4.58   2.98
Ti II   3224.24   1.58   0.05   45.6  -4.85   2.68
Ti II   3228.61   1.08  -0.21   49.3  -4.82   2.56
Ti II   3229.42   1.13  -0.12   49.2  -4.82   2.52
Ti II   3232.28   1.12  -0.22   48.9  -4.82   2.59
Ti II   3236.57   0.03   0.24  105.8  -4.49   2.85
Ti II   3239.04   0.01   0.04   98.2  -4.52   2.87
Ti II   3241.98   0.00  -0.03   96.2  -4.53   2.88
Ti II   3263.68   1.17  -1.14   17.9  -5.26   2.72
Ti II   3276.77   1.18  -0.89   22.9  -5.16   2.63
Ti II   3278.29   1.23  -0.26   45.8  -4.86   2.65
Ti II   3278.92   1.08  -0.22   49.3  -4.82   2.57
Ti II   3282.33   1.22  -0.34   40.6  -4.91   2.58
Ti II   3287.65   1.89   0.45   45.0  -4.86   2.54
Ti II   3288.14   0.14  -1.99   23.0  -5.16   2.73
Ti II   3288.43   1.24  -1.25   14.2  -5.37   2.77
Ti II   3315.32   1.22  -0.64   31.8  -5.02   2.65
Ti II   3321.70   1.23  -0.34   39.6  -4.92   2.56
Ti II   3329.45   0.14  -0.26   84.9  -4.59   2.88
Ti II   3340.34   0.11  -0.53   75.0  -4.65   2.80
Ti II   3341.88   0.57   0.34   84.9  -4.59   2.68
Ti II   3372.79   0.01   0.28  111.9  -4.48   2.84
Ti II   3380.28   0.05  -0.54   80.9  -4.62   2.94
Ti II   3387.83   0.03  -0.41   83.6  -4.61   2.87
Ti II   3388.75   1.24  -1.02   18.5  -5.26   2.67
Ti II   3409.81   0.03  -1.91   31.3  -5.04   2.74
Ti II   3443.37   2.05  -1.21    3.2  -6.03   2.75
Ti II   3456.38   2.06  -0.11   19.4  -5.25   2.56
Ti II   3491.05   0.11  -1.10   62.5  -4.75   2.88
Ti II   3500.33   0.12  -2.04   23.4  -5.17   2.75
Ti II   3504.89   1.89   0.38   42.2  -4.92   2.50
Ti II   3535.41   2.06   0.01   22.8  -5.19   2.54
Ti II   3561.58   0.57  -2.04    9.6  -5.57   2.69
Ti II   3596.05   0.61  -1.07   40.6  -4.95   2.67
Ti II   3659.76   1.58  -0.54   23.6  -5.19   2.55
Ti II   3662.23   1.57  -0.54   22.7  -5.21   2.51
Ti II   3741.64   1.58  -0.07   44.0  -4.93   2.55
Ti II   3757.68   1.57  -0.44   31.2  -5.08   2.61
Ti II   3759.29   0.61   0.28   93.4  -4.61   2.68
Ti II   3761.32   0.57   0.18   93.0  -4.61   2.73
Ti II   3761.87   2.59  -0.42    5.8  -5.81   2.63
Ti II   3776.05   1.58  -1.24    7.7  -5.69   2.64
Ti II   3813.39   0.61  -1.89   17.2  -5.35   2.77
Ti II   3900.54   1.13  -0.29   61.7  -4.80   2.76
Ti II   4012.38   0.57  -1.78   24.4  -5.22   2.81
Ti II   4025.13   0.61  -2.11   10.3  -5.59   2.70
Ti II   4028.34   1.89  -0.92    9.4  -5.63   2.68
Ti II   4163.64   2.59  -0.13   13.2  -5.50   2.71
Ti II   4320.95   1.17  -1.88    6.6  -5.82   2.76
Ti II   4330.70   1.18  -2.09    3.8  -6.06   2.73
Ti II   4395.03   1.08  -0.54   59.1  -4.87   2.82
Ti II   4395.84   1.24  -1.93    5.5  -5.90   2.80
Ti II   4399.77   1.24  -1.20   21.5  -5.31   2.75
Ti II   4418.33   1.24  -1.99    4.7  -5.97   2.77
Ti II   4443.80   1.08  -0.71   50.4  -4.95   2.77
Ti II   4444.55   1.12  -2.20    3.9  -6.06   2.78
Ti II   4468.49   1.13  -0.63   52.8  -4.93   2.78
Ti II   4501.27   1.12  -0.77   46.6  -4.99   2.77
Ti II   4571.97   1.57  -0.31   45.1  -5.01   2.71
Ti II   5129.16   1.89  -1.34    5.0  -6.01   2.73
Ti II   5336.79   1.58  -1.60    6.0  -5.95   2.77"""
        ewlist = []
        threeDlist = []
        wllist = []
        nlte3dlist = []
        abund_list = []
        ltelist = []
        energylist = []
        gflist = []
        for line in (scott_data.split("\n")):
            scott_split = line.split()
            print(scott_split)

            wl = float(scott_split[2])
            wllist.append(wl)
            rtc = 1.0973731569E5  # Conversion factor between Rydberg and cm^(-1)
            rte = 13.605698  # Conversion factor between Rydberg and eV

            energylist.append((4.135667E-15 * 299792458 / (wl*1E-10) ))
            gflist.append(float(scott_split[4]))
            if ti2:
                ew = float(scott_split[-3])
            elif not ti2:
                ew = float(scott_split[-3])
            # To mAA from AA
            ew = np.log10(ew / 1000)
            ewlist.append(ew)

            nlte3d = float(scott_split[-1])
            nlte3dlist.append(nlte3d)
            abund_list.append(nlte3d)
            lte = float(scott_split[-1])
            ltelist.append(lte)

        extra_dict = {}
        extra_dict['ew'] = ewlist
        extra_dict['loggf'] = gflist
        extra_dict['energy'] = energylist
        return wllist, ewlist, abund_list, abund_list, extra_dict



lines1, equiv_width1, lteabund1, nlteabund1, extra_dict1 = ti1get_chris_data()
lines, equiv_width2, lteabund2, nlteabund2, extra_dict2 = ti2get_chris_data()
energy2 = extra_dict2['energy']
energy1 = extra_dict1['energy']

# Changing the abundances of the lines that literature found to either lte or nlte depending on our run type. ONly useful
# for the pat scott data which has both





ti2lte = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "LTEti2_abundance_comparisons", "rb"))
ti2nlte = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "NLTEti2_abundance_comparisons",
         "rb"))
ti1lte = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "LTE_abundance_comparisons", "rb"))
ti1nlte = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "NLTE_abundance_comparisons", "rb"))

ti1lines = ti1lte['lines']
ti1lteabund = ti1lte['our_abund']
ti1nlteabund = ti1nlte['our_abund']

ti2lines = ti2lte['lines']
ti2lteabund = ti2lte['our_abund']
ti2nlteabund = ti2nlte['our_abund']

# equiv widths
ti2nlte_equiv_width_dict = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "Ti2_EquivWidths_AbundancesNLTE.pkl",
         "rb"))
ti2lte_equiv_width_dict = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "Ti2_EquivWidths_AbundancesLTE.pkl",
         "rb"))

nlte_equiv_width_dict = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "EquivWidths_AbundancesNLTE.pkl",
         "rb"))
lte_equiv_width_dict = pickle.load(
    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "EquivWidths_AbundancesLTE.pkl",
         "rb"))

# collecting ews to mean them, nlte ti1 anbd ti2
ews = []
# mean ews
mean_ew = []
ti1reg_ew = []
# Get the mean equiv widths for ti1 and 2 over all abundaneces
for r in range(len(nlte_equiv_width_dict[list(nlte_equiv_width_dict)[0]])):
    ews = []
    regews = []

    for abund in nlte_equiv_width_dict:
        # Logged previously so we unlog it now
        ew = 10 ** nlte_equiv_width_dict[abund][r]
        regews.append(ew)
        #reduced ew

        ews.append(np.log(ew / ti1lines[r]))
    mean_ew.append(sum(ews) / len(nlte_equiv_width_dict))
    ti1reg_ew.append(sum(regews)/len(nlte_equiv_width_dict))
# mean ews
ti2reg_ew = []

ti2mean_ew = []
for r in range(len(ti2nlte_equiv_width_dict[list(ti2nlte_equiv_width_dict)[0]])):
    ews = []
    regews = []

    for abund in ti2nlte_equiv_width_dict:
        ew = 10 ** ti2nlte_equiv_width_dict[abund][r]
        regews.append(ew)

        ews.append(np.log10(ew / ti2lines[r]))
    ti2mean_ew.append(sum(ews) / len(ti2nlte_equiv_width_dict))
    ti2reg_ew.append(sum(regews)/len(nlte_equiv_width_dict))




#for lte
# mean ews
ltemean_ew = []
lteti1reg_ew = []
# Get the mean equiv widths for ti1 and 2
for r in range(len(lte_equiv_width_dict[list(ti2lte_equiv_width_dict)[0]])):
    ews = []
    regews = []
    for abund in lte_equiv_width_dict:
        ew = 10**lte_equiv_width_dict[abund][r]
        regews.append(ew)

        ews.append(np.log10(ew/ti1lines[r]))
    ltemean_ew.append(sum(ews) / len(lte_equiv_width_dict))
    lteti1reg_ew.append(sum(regews)/len(nlte_equiv_width_dict))

ews = []
# mean ews
lteti2mean_ew = []
ti2ltereg_ew = []
for r in range(len(ti2lte_equiv_width_dict[list(ti2lte_equiv_width_dict)[0]])):
    ews = []
    regews = []

    for abund in ti2lte_equiv_width_dict:
        ew = 10**ti2lte_equiv_width_dict[abund][r]
        regews.append(ew)

        ews.append(np.log10(ew/ti2lines[r]))
    lteti2mean_ew.append(sum(ews) / len(ti2lte_equiv_width_dict))
    ti2ltereg_ew.append(sum(regews)/len(nlte_equiv_width_dict))


ewrun = False
if ewrun:
    fig = plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("$Wavelength _\AA$ (Ti I)")
    ax.set_ylabel("$Equivalent$ $Width$")
    ax2.set_xlabel("$Wavelength _\AA$ (Ti II)")

    aa = ax.plot(ti1lines, ti1reg_ew, marker='s', label="NLTE Ti I", color="blue", markerfacecolor='none',
                 linestyle="None")

    bb = ax2.plot(ti2lines, ti2reg_ew, marker='o', label="NLTE Ti II", color="red", markerfacecolor='none',
                  linestyle="None")

    aa2 = ax.plot(ti1lines, lteti1reg_ew, marker='s', label="LTE Ti I", color="blue",
                 linestyle="None")

    bb2 = ax2.plot(ti2lines, ti2ltereg_ew, marker='o', label="LTE Ti II", color="red",
                  linestyle="None")
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    lns = aa + bb + aa2 + bb2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.title("Wavelength vs Equivalent Width")
    plt.show()







    fig = plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("$Excitation potential (Ti I)")
    ax.set_ylabel("$Equivalent$ $Width$")
    ax2.set_xlabel("Excitation potential (Ti II)")

    aa = ax.plot(energy1, ti1reg_ew, marker='s', label="NLTE Ti I", color="blue", markerfacecolor='none',
                 linestyle="None")

    bb = ax2.plot(energy2, ti2reg_ew, marker='o', label="NLTE Ti II", color="red", markerfacecolor='none',
                  linestyle="None")

    aa2 = ax.plot(energy1, lteti1reg_ew, marker='s', label="LTE Ti I", color="blue",
                  linestyle="None")

    bb2 = ax2.plot(energy2, ti2ltereg_ew, marker='o', label="LTE Ti II", color="red",
                   linestyle="None")
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    lns = aa + bb + aa2 + bb2
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc="upper left")
    plt.title("Excitation potential vs Equivalent Width")
    plt.show()

else:
    wl_coded = True
    if wl_coded:

        import matplotlib.cm as cm
        from matplotlib import cm
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        import matplotlib.colors as mcolors

        colors = plt.cm.rainbow(np.linspace(0, 1, len(energy2)))
        print(ti2lines)
        print(max(ti2lines) - min(ti2lines))
        print(len(ti2lines))

        # plotting wl heatmap for ti2
        ii = 0
        for i in np.linspace(0,1,len(ti2lines)):
            print(i, ii)
            print(ti2nlteabund[ii],"\n")
            plt.scatter(energy2[ii],ti2nlteabund[ii], color = plt.cm.RdYlBu(i),edgecolors="black")
            plt.scatter(energy2[ii],ti2lteabund[ii], marker="s", color = plt.cm.RdYlBu(i),edgecolors="black")
            ii +=1
        plt.scatter(energy2[0], ti2lteabund[0], marker="s", label = "LTE", color= plt.cm.RdYlBu(0),  edgecolors="black", facecolors='none')
        plt.scatter(energy2[0], ti2nlteabund[0], marker="o", label = "NLTE", color= plt.cm.RdYlBu(0),  edgecolors="black", facecolors='none')
        cmappable = ScalarMappable(norm=Normalize(min(ti2lines), max(ti2lines)), cmap=plt.cm.RdYlBu)
        plt.colorbar(cmappable)
        plt.title("Excitation potential vs Derived abundance Ti II")
        plt.ylabel("Abundance")
        plt.legend()
        plt.xlabel("Excitation potential")
        plt.show()



        # for ti1
        ii = 0
        for i in np.linspace(0,1,len(ti1lines)):
            print(i, ii)
            print(ti1nlteabund[ii],"\n")
            plt.scatter(energy2[ii],ti1nlteabund[ii], color = plt.cm.RdYlBu(i),edgecolors="black")
            plt.scatter(energy2[ii],ti1lteabund[ii], marker="s", color = plt.cm.RdYlBu(i),edgecolors="black")
            ii +=1
        plt.scatter(energy2[0], ti1lteabund[0], marker="s", label = "LTE", color= plt.cm.RdYlBu(0),  edgecolors="black", facecolors='none')
        plt.scatter(energy2[0], ti1nlteabund[0], marker="o", label = "NLTE", color= plt.cm.RdYlBu(0),  edgecolors="black", facecolors='none')
        cmappable = ScalarMappable(norm=Normalize(min(ti1lines), max(ti1lines)), cmap=plt.cm.RdYlBu)
        plt.colorbar(cmappable)
        plt.title("Excitation potential vs Derived abundance Ti I")
        plt.ylabel("Abundance")
        plt.legend()
        plt.xlabel("Excitation potential")
        plt.show()



    else:
        plt.scatter(energy1,ti1lteabund,  color= "blue", facecolors='none', label="Ti I LTE")
        plt.scatter(energy1,ti1nlteabund,  color= "blue", label="Ti I NLTE")
        plt.scatter(energy2,ti2lteabund, color = "red", facecolors='none',label = "Ti II LTE")
        plt.scatter(energy2,ti2nlteabund, color = "red", label = "Ti II NLTE")
        plt.title("Wavelength vs Excitation potential")
        plt.legend()
        plt.ylabel("Abundance")
        plt.xlabel("Excitation potential")
        plt.show()
