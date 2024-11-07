# here we compare the abundances WE CALCULATED with our lines vs scotts by opening a singlke dict
# with the results per line from nlte and lte. So run it AFTER abundance_comparisons

import matplotlib.pyplot as plt
import numpy as np
import pickle
star_number = 7
            #  0        1          2          3         4        5                 6               7
star_list = ['sun', 'Arcturus', '122563', '140283', '84937', "140283Karin", "arcturusKarin", "122563_3D"]
star = star_list[star_number]
k = "0001"
if k != "0005":
    pass
else:
    k=""

ti2 = False
# Load the LTE and NLTE values of scott/our abundances/lines etc
if ti2:
    # there are noaccurate  nlte values for ti2
    lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTEti2_abundance_comparisons", "rb"))
    nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTEti2_abundance_comparisons", "rb"))

else:
    lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTE_abundance_comparisons", "rb"))
    nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTE_abundance_comparisons", "rb"))

# Should actually be the same I think except for Ti I in the sun...
scottabund = lte['literature_abund']
nltescottabund = nlte['literature_abund']
nlteabund= nlte['our_abund']
lteabund = lte['our_abund']

print(len(scottabund))
print(len(lteabund))

nltediff = [nlte_val - scott_val for scott_val, nlte_val in zip(nltescottabund, nlteabund)]
ltediff = [lte_val - scott_val for scott_val, lte_val in zip(scottabund, lteabund)]


lines = lte['lines']
print(len(lines))
l_ew=lte['ew']
l_energy=lte['energy']
l_gf = lte['loggf']
print(len(l_energy))

nl_ew=nlte['ew']
nl_energy=nlte['energy']
nl_gf = nlte['loggf']
plt.scatter(ltediff, l_ew, label = "EW")
plt.legend()
plt.xlabel("Difference in abundance from literature")
plt.ylabel("Variables")
plt.show()

plt.scatter(ltediff, l_gf, label = "Log gf")
plt.legend()
plt.xlabel("Difference in abundance from literature")
plt.ylabel("Variables")
plt.show()



plt.scatter( ltediff,l_energy, label = "Energy")
plt.legend()
plt.xlabel("Difference in abundance from literature")
plt.ylabel("Variables")
plt.show()
