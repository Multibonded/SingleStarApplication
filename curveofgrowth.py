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

# which stellar atmos are we looking at
star_number = 2
            #  0        1          2          3         4        5                 6               7     8(10% electron rates)     9    10                      11                  12              13
star_list = ['sun', 'Arcturus', '122563', '140283', '84937', "140283Karin", "arcturusKarin", "122563_3D", "122563_10th", "122563_100th", "122563_0kaulakys", "Arcturus_100th", "Arcturus_10th", "Arcturus_0electrons"
    , "Sun_0kaulakys100electrons", "Arcturus_0kaulakys100electrons", "140283_0kaulakys100electrons",     "84397_0kaulakys100electrons", "122563_0kaulakys100electrons", "chris"]
        # 14                             15                              16                                      17                          18
star = star_list[star_number]
ti2l = [True, False]
ltel = [False, True]
worked = 0
def get_122():
    # blended
    #             3717.4        15.970     0.251
    #
    data144ti1ew = """3598.7         9.777     0.220
3635.5         0.000     0.000
3642.7         0.000     0.000
3725.2         4.023     0.142
3729.8        62.128     0.445
3741.1        60.851     0.428
3752.9        74.038     0.929
3904.8        26.459     0.363
3924.5        29.245     0.373
3929.9         0.000     0.000
3947.8        31.212     0.770
3956.3         0.000     0.000
3958.2        73.019     2.408
3981.8         0.000     0.000
3989.8         0.000     0.000
3998.6        75.584     1.503
4008.9        26.877     0.330
4024.6         0.000     0.000
4287.4        19.025     0.231
4298.7         0.000     0.000
4305.9         0.000     0.000
4427.1         0.000     0.000
4449.1         5.708     0.152
4518.0        18.905     0.202
4534.8        43.721     0.333
4535.6        36.210     0.479
4548.8        17.749     0.190
4555.5        13.090     0.158
4617.3         8.111     0.138
4691.3         0.000     0.000
4981.7        62.468     0.396
5014.3         0.000     0.000
5022.9        17.917     0.195
5036.5        12.377     0.162
5173.7        32.679     0.316
5193.0        37.243     0.333"""
    data144ti2ew = """3500.3        99.922     2.727
3504.9         0.000     0.000
3520.3         0.000     0.000
3533.9        15.915     0.250
3535.4        66.974     1.476
3561.6         0.000     0.000
3573.7         0.000     0.000
3596.0         0.000     0.000
3706.2         0.000     0.000
3741.6         0.000     0.000
3748.0         0.000     0.000
3757.7        83.136     2.103
3759.3       209.168     1.279
3761.3       194.250     4.133
3761.9        34.370     0.552
3774.6        55.363     3.764
3776.1        51.734     0.689
3786.3         0.000     0.000
3813.4       100.697     5.068
3882.3         0.000     0.000
3900.5       133.115     0.853
3913.5         0.000     0.000
3982.0         0.000     0.000
3987.6        48.878     0.636
4012.4       103.308     2.660
4025.1        90.895     0.846
4028.3        54.402     0.386
4053.8        41.757     0.460
4161.5        51.710     0.594
4163.6        45.062     0.480
4171.9         0.000     0.000
4290.2         0.000     0.000
4300.0         0.000     0.000
4301.9         0.000     0.000
4312.9         0.000     0.000
4315.0         0.000     0.000
4316.8        18.131     0.207
4321.0         0.000     0.000
4330.2         0.000     0.000
4330.7        49.565     0.624
4367.7         0.000     0.000
4391.0        46.454     0.470
4395.8        51.120     0.621
4409.5        21.144     0.293
4418.3        51.096     0.348
4444.6        47.123     0.503
4488.3        10.120     0.159
4493.5        18.654     0.202
4518.3        32.997     0.353
4545.1        33.627     0.282
4549.6         0.000     0.000
4552.3         0.000     0.000
4583.4        16.976     0.221
4609.3         5.332     0.116
4657.2        34.988     0.629
4708.7        36.021     0.280
4719.5         6.008     0.130
4762.8         0.000     0.000
4763.9        35.184     0.346
4764.5        18.544     0.198
4798.5        28.559     0.650
4865.6        24.999     0.249
4874.0         4.867     0.123
4911.2         8.161     0.129
5005.2         0.000     0.000
5013.7        22.092     0.212
5069.1         0.000     0.000
5072.3         0.000     0.000
5129.2        40.354     0.823
5185.9        35.635     0.319
5211.5         5.966     0.111
5336.8        53.625     0.370
5381.0        34.910     0.280
5396.2         5.427     0.112
5418.8        26.487     0.208
6680.1         0.000     0.000
9432.2         0.000     0.000"""

    """   1- 11 F11.2  0.1nm    lambda             Wavelength in Angstroms, in air units
  12- 18 A7     ---      Ion                Species identification
  19- 27 F9.3   eV       chi                lower excitation level
  28- 35 F8.2   [-]      log(gf)            Log of the oscillator strength
  36- 43 F8.1   10-4nm   EW-122563          Equivalent width, HD 122563, milli-Angstroms
  44- 51 F8.1   10-4nm   EW-140283          ? equivalent width, HD 140283, milli-Angstroms
  52- 59 F8.2   ---      [X/Fe]122563       relative abundance, HD 122563(1)
  60- 67 F8.2   ---      [X/Fe]140283       ? relative abundance, HD 140283(1)"""
    dataabundti1 = """4449.14   Ti I    1.890    0.47     6.4   0.14
    4450.89   Ti I    1.880    0.32     3.6       1     0.02        
    4453.31   Ti I    1.430   -0.03     9.3       1     0.26        
    4455.32   Ti I    1.440    0.13     8.7      1      0.08        
    4457.43   Ti I    1.460    0.26    15.5       1     0.25        
    4512.73   Ti I    0.830   -0.40    14.0     2.0    0.09
    4518.02   Ti I    0.820   -0.25    19.0     3.3    0.08
    4527.30   Ti I    0.810   -0.45    14.4     2.7    0.12
    4533.24   Ti I    0.850    0.54    51.9    15.5   -0.05 
    4534.78   Ti I    0.830    0.35    42.5    11.1   -0.03 
    4535.57   Ti I    0.820    0.14    34.2     7.6    0.02
    4544.69   Ti I    0.820   -0.45    13.8     2.0    0.10
    4548.76   Ti I    0.820   -0.28    16.5      1      0.03        
    4555.48   Ti I    0.850   -0.39    11.7     1       0.00        
    4759.27   Ti I    2.250    0.59     3.8      1      0.22        
    4870.12   Ti I    2.250    0.44     2.0      1      0.06        
    4885.08   Ti I    1.890    0.41     5.0      1      0.07        
    4913.61   Ti I    1.870    0.22     4.7      1      0.21        
    4991.07   Ti I    0.830    0.45    63.0    15.0    0.14
    4999.50   Ti I    0.820    0.32    50.0      1      0.06        
    5007.21   Ti I    0.820    0.17    51.5    13.5    0.23
    5009.65   Ti I    0.020   -2.20     4.1      1      0.26        
    5016.16   Ti I    0.850   -0.48    12.9      1      0.11        
    5024.84   Ti I    0.820   -0.53    13.5      1      0.14        
    5036.46   Ti I    1.440    0.14     8.6      1      0.03        
    5038.40   Ti I    1.430    0.02     8.4      1      0.13        
    5039.96   Ti I    0.020   -1.08    32.9      1      0.16        
    5064.65   Ti I    0.050   -0.94    35.0     4.7    0.09
    5145.46   Ti I    1.460   -0.54     2.5       1     0.17        
    5147.48   Ti I    0.000   -1.94     6.4      1      0.16        
    5173.74   Ti I    0.000   -1.06    31.6     3.7    0.08
    5192.97   Ti I    0.020   -0.95    36.3     4.5    0.08
    5210.38   Ti I    0.050   -0.82    41.4     6.7    0.06
    5219.70   Ti I    0.020   -2.22     4.0      1      0.25        
    5866.45   Ti I    1.070   -0.73     5.5     1       0.17        
    6258.10   Ti I    1.440   -0.39     4.1      1      0.15        
    6258.71   Ti I    1.460   -0.28     6.9      1      0.29        
    6261.10   Ti I    1.430   -0.53     3.2      1      0.16"""

    dataabundti2 = """4395.84  Ti II    1.240   -1.93    48.6     5.9    0.23  
    4399.77  Ti II    1.240   -1.20    90.0    29.5    0.22
    4409.52  Ti II    1.230   -2.53    21.8     1.8    0.31
    4411.07  Ti II    3.090   -0.65     9.0      1      0.26        
    4411.93  Ti II    1.220   -2.62    19.8     1.6    0.34
    4418.33  Ti II    1.240   -1.99    50.9     6.2    0.32
    4421.94  Ti II    2.060   -1.64    14.3      1      0.22        
    4432.10  Ti II    1.240   -3.08     7.3      1      0.32        
    4443.80  Ti II    1.080   -0.71   127.7    61.9    0.33
    4444.55  Ti II    1.110   -2.20    47.5     5.4    0.32
    4450.48  Ti II    1.080   -1.52    84.1    24.0    0.21
    4468.49  Ti II    1.130   -0.63   128.0    62.2    0.31
    4469.15  Ti II    1.080   -2.55    33.0      1      0.38        
    4488.32  Ti II    3.120   -0.50     8.9      1      0.13        
    4493.52  Ti II    1.080   -2.78    19.0      1      0.29        
    4501.27  Ti II    1.110   -0.77   123.8    57.0    0.31 
    4518.33  Ti II    1.080   -2.56    33.8     3.6    0.39
    4529.48  Ti II    1.570   -1.75    42.6     6.3    0.35
    4545.13  Ti II    1.130   -2.45    31.1     2.2    0.29
    4571.97  Ti II    1.570   -0.31   115.8    52.4    0.25
    4583.41  Ti II    1.160   -2.84    16.2      1      0.37        
    4609.27  Ti II    1.180   -3.32     5.0      1      0.31        
    4708.66  Ti II    1.240   -2.35    33.3     3.5    0.35
    4719.51  Ti II    1.240   -3.32     5.0       1     0.38        
    4764.52  Ti II    1.240   -2.69    18.1       1     0.35        
    4849.17  Ti II    1.130   -2.96    14.0       1     0.35        
    4865.61  Ti II    1.110   -2.70    25.2       1     0.38        
    4874.01  Ti II    3.090   -0.86     4.1       1     0.09        
    4911.19  Ti II    3.120   -0.64     8.3       1     0.22        
    5005.17  Ti II    1.560   -2.73     7.3       1     0.34        
    5013.69  Ti II    1.580   -2.14    21.0       1     0.29        
    5185.90  Ti II    1.890   -1.41    35.8     5.0    0.24
    5211.53  Ti II    2.590   -1.41     5.7       1     0.16        
    5268.61  Ti II    2.600   -1.61     4.4       1     0.25        
    5336.79  Ti II    1.580   -1.60    52.1       1     0.29        
    5381.02  Ti II    1.560   -1.97    34.6     3.6    0.36
    5396.25  Ti II    1.580   -3.18     2.9     3.1    0.38
    5418.77  Ti II    1.580   -2.13    26.3      1      0.37"""

    ewdict = {}

    energylist = []
    gflist = []
    wllist = []
    equivwidthlist = []
    abund_list = []

    if ti2:
        # abund data from the literature
        data144 = dataabundti2
        # ew data
        tiews = data144ti2ew

    else:
        data144 = dataabundti1
        tiews = data144ti1ew
    # Making the dictionary of equiv widths with karin's ews to use its wavelengths to compare to the ones in lit.
    for line in (tiews.split("\n")):
        ewdict[float(line.split()[0])] = float(line.split()[1])
    #  36- 43 F8.1   10-4nm   EW-122563          Equivalent width, HD 122563, milli-Angstroms
    for line in (tiews.split("\n")):
        split = line.split()
        wl = float(split[0])
        # Only use the ones in both Karin's and lit
        gf = 0
        equivwidth = float(split[1])

        if equivwidth == 0.0:
            continue
        gflist.append(gf)
        energylist.append(0)

        wllist.append(wl)
        # Go from mAA to AA

        equivwidthlist.append((equivwidth / 1000))
        # Lte abundance
        lte = float(0)
        abund_list.append(lte)

    return wllist, equivwidthlist, abund_list, gflist, energylist


for ti2 in ti2l:
    worked+=1
    for lte in ltel:
        wllist, equivwidthlist, abund_list, gflist, energylist = get_122()

        print("\n\n")
        worked+=1
#abundances = [["445", "455", "460", "470", "480", "485", "490", "495", "502", "520", "545"]]
        # k0001
        k = "0001"
        # WE did more abundances for the sun at k=0.005
        # abunds for sun, arcturus, 122563, 140283, and 84937
        abundances_all = [["450", "460", "480", "490", "500"], ["460", "470", "490", "500"],
                          ["200", "220", "230", "250", "270"],
                          ["250", "270", "280", "300"], ["290", "310", "320", "330"],
                          ["250", "260", "270", "280", "290", "300"], ["460", "470", "480", "490", "500"],
                          ["200", "220", "230", "240", "250", "260"], ["210", "230", "250", "270"],
                          ["210", "230", "250", "270"], ["210", "230", "250", "270"], ["460", "480", "500"],
                          ["460", "480", "500"], ["460", "480", "500"],
                          ["450", "460", "480", "490", "510"], ["460", "470", "490", "500"],
                          ["250", "270", "280", "300"], ["280", "290", "300", "320", "330"],
                          ["210", "230", "250", "270"], ["265", "270", "275", "280", "290"]]
        abundances_string = abundances_all[star_number]




        # Turning them into floats.
        abundances = [float(abund)/100 for abund in abundances_string]


        if lte:
            if ti2:
                equiv_width_dict = pickle.load(open(
                    r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "Ti2_EquivWidths_AbundancesLTE.pkl",
                    "rb"))
            else:
                equiv_width_dict = pickle.load(
                    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "EquivWidths_AbundancesLTE.pkl",
                         "rb"))
        else:
            if ti2:
                equiv_width_dict = pickle.load(open(
                    r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "Ti2_EquivWidths_AbundancesNLTE.pkl",
                    "rb"))
            else:
                equiv_width_dict = pickle.load(
                    open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "EquivWidths_AbundancesNLTE.pkl",
                         "rb"))

        mdict = []
        ewlist = []
        wl = []
        number =                                                                             14
        for metallicity in equiv_width_dict:
            mdict.append(float(metallicity)/100)
            ewlist.append(equiv_width_dict[metallicity][number])
            wl.append(equiv_width_dict[metallicity][number]/wllist[number])
        print(lte)
        if lte:
            print(1)
            plt.plot(mdict, wl, label = "LTE")
        else:
            print(2)
            plt.plot(mdict, wl, linestyle="dashed", color="red", label = "NLTE")

    plt.title(wllist[number])
    plt.ylabel("Reduced equivalent width")
    plt.xlabel("[Fe/H]")
    plt.yscale('log')
    plt.legend()
    #plt.xscale('log')
    plt.show()
