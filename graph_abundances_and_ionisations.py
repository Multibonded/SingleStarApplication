import matplotlib.pyplot as plt
import numpy as np
import pickle
star_number = 0
            #  0        1          2          3         4        5                 6               7     8(10% electron rates)     9    10                      11                  12              13
star_list = ['sun', 'Arcturus', '122563', '140283', '84937', "140283Karin", "122563_0kaulakys100electrons"]
star_list = ['sun', '122563', '140283', '84937', "122563_0kaulakys100electrons"]
metallicity = [0, -2.5, -2.36, -2.06, -2.5]

k = "0001"

dwarfmetal = []
giantmetal = []
sunmetal = []

dwarfion = []
giantion = []
sunion = []

nltedwarfmetal = []
nltegiantmetal = []
nltesunmetal = []

nltedwarfion = []
nltegiantion = []
nltesunion = []


for starnumb in range(len(star_list)):
    star = star_list[starnumb]
    metal = metallicity[starnumb]

    ti2lte = pickle.load(
        open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "LTEti2_abundance_comparisons",
             "rb"))
    ti2nlte = pickle.load(
        open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "NLTEti2_abundance_comparisons",
             "rb"))
    ti1lte = pickle.load(
        open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "LTE_abundance_comparisons",
             "rb"))
    ti1nlte = pickle.load(
        open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "NLTE_abundance_comparisons",
             "rb"))

    ti1lines = ti1lte['lines']
    ti1lteabund = ti1lte['our_abund']
    ti1nlteabund = ti1nlte['our_abund']

    ti2lines = ti2lte['lines']
    ti2lteabund = ti2lte['our_abund']
    ti2nlteabund = ti2nlte['our_abund']

    if star == "sun":
        weightti1 = ti1lte['weight']
        weightti2 = ti2lte['weight']
    else:
        weightti1 = np.ones(len(ti1lines))
        weightti2 = np.ones(len(ti2lines))


    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    nltemeantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    nltemeantiII1d = sum(weighted_ti2_1d) / sum(weightti2)


    # same but for lte now

    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d) / sum(weightti2)
    if star == "122563_0kaulakys100electrons":
        star = "122563_modified"
    if "122563" in star:
        startype = "giant"
        giantmetal.append(metal)
        giantion.append(meantiI1d-meantiII1d)
        nltegiantmetal.append(metal)
        nltegiantion.append(nltemeantiI1d-nltemeantiII1d)
    else:
        startype = "dwarf"
        dwarfmetal.append(metal)
        dwarfion.append(meantiI1d-meantiII1d)
        nltedwarfmetal.append(metal)
        nltedwarfion.append(nltemeantiI1d - nltemeantiII1d)
    print(meantiI1d-meantiII1d)

print(nltegiantion)
plt.scatter(dwarfmetal, dwarfion, label ="dwarf LTE", color = "red")
plt.scatter(nltedwarfmetal, nltedwarfion, label = "giant LTE", color = "red", facecolor='none')

plt.scatter(giantmetal, giantion, label = "giant LTE", color="blue")
plt.scatter(nltegiantmetal, nltegiantion, label = "giant NLTE", color = "blue", facecolor='none')
plt.legend()
plt.show()