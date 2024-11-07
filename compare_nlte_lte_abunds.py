# here we compare the abundances WE CALCULATED with our lines vs scotts by opening a singlke dict
# with the results per line from nlte and lte. So run it AFTER abundance_comparisons

import matplotlib.pyplot as plt
import numpy as np
import pickle
star_number = 3
            #  0        1          2          3         4        5                 6               7     8(10% electron rates)     9    10                      11                  12              13
star_list = ['sun', 'Arcturus', '122563', '140283', '84937', "140283Karin", "arcturusKarin", "122563_3D", "122563_10th",  "122563_100th", "122563_0kaulakys", "Arcturus_100th", "Arcturus_10th", "Arcturus_0electrons"
    , "Sun_0kaulakys100electrons", "Arcturus_0kaulakys100electrons", "140283_0kaulakys100electrons",             "84397_0kaulakys100electrons", "122563_0kaulakys100electrons", "chris"]
        # 14                             15                              16                                      17                               18
star = star_list[star_number]


#abundances = [["445", "455", "460", "470", "480", "485", "490", "495", "502", "520", "545"]]
# k0001
k = "0001"





ti2 = False
# Load the LTE and NLTE values of scott/our abundances/lines etc
if ti2:
    # there are noaccurate  nlte values for ti2
    lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTEti2_abundance_comparisons", "rb"))
    nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTEti2_abundance_comparisons", "rb"))

else:
    lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTE_abundance_comparisons", "rb"))
    nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTE_abundance_comparisons", "rb"))

# equiv widths
ti2nlte_equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"Ti2_EquivWidths_AbundancesNLTE.pkl", "rb"))
ti2lte_equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"Ti2_EquivWidths_AbundancesLTE.pkl", "rb"))

nlte_equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"EquivWidths_AbundancesNLTE.pkl", "rb"))
lte_equiv_width_dict = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"EquivWidths_AbundancesLTE.pkl", "rb"))






' lines are the same for both nlte and lte'
lines = lte['lines']
lteabund = lte['our_abund']

# Should actually be the same I think except for Ti I in the sun...
literature_abund = lte['literature_abund']
nlte_literature_abund = lte['literature_abund']

nlteabund= nlte['our_abund']

list1, list2 = zip(*sorted(zip(nlteabund, lines)))



# Only the sun had these weird 3d difference additions
if star =="sun":
    lte_threeDdiff = lte['3d_diff']
    nlte_threeDdiff = nlte['3d_diff']
# position of text on figure
if ti2:
    y = 0.34
else:
    y = 0.64

def compare_ews():


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
    # collecting ews to mean them
    ews = []
    # mean ews
    mean_ew = []
    # Get the mean equiv widths for ti1 and 2
    for r in range(len(nlte_equiv_width_dict[list(nlte_equiv_width_dict)[0]])):
        ews = []
        for abund in nlte_equiv_width_dict:
            ew = 10**nlte_equiv_width_dict[abund][r]
            ews.append(np.log(ew/ti1lines[r]))
        mean_ew.append(sum(ews) / len(nlte_equiv_width_dict))



    # mean ews
    ti2mean_ew = []
    for r in range(len(ti2nlte_equiv_width_dict[list(ti2nlte_equiv_width_dict)[0]])):
        ews = []
        for abund in ti2nlte_equiv_width_dict:
            ew = 10**ti2nlte_equiv_width_dict[abund][r]
            ews.append(np.log10(ew/ti2lines[r]))
        ti2mean_ew.append(sum(ews) / len(ti2nlte_equiv_width_dict))


    remv = (np.where(np.asarray(mean_ew) > -8.4))
    print(remv[0])
    print("lines to remove for arcturus", np.asarray(lines)[remv])
    list1, list2 = zip(*sorted(zip(mean_ew, lines)))
    print(list1, "\n", list2)

    remv = (np.where(np.asarray(ti2mean_ew) > -3.62))
    print(remv[0])
    print("lines to remove for arcturus ti2", np.asarray(ti2lines)[remv])
    list1, list2 = zip(*sorted(zip(mean_ew, ti2lines)))
    print(list1, "\n", list2)





    if star == "sun":
        weightti1 = ti1lte['weight']
        weightti2 = ti2lte['weight']
    else:
        weightti1 = np.ones(len(ti1lines))
        weightti2 = np.ones(len(ti2lines))

    """for x in range(len(ti2lines)):
        print(ti2lines[x])
        print(ti2lteabund[x])
        print(ti2lte_threeDdiff[x])
        print("\n")
    exit()"""

    fig = plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("$W_{λ,red}$ (Ti I)")
    ax.set_ylabel("$A(Ti)$")
    ax2.set_xlabel("$W_{λ,red}$ (Ti II)")

    aa = ax.plot(mean_ew, ti1nlteabund, marker='s', label="1D NLTE Ti I", color="blue", linestyle="None",
                 markerfacecolor='none')

    bb = ax2.plot(ti2mean_ew, ti2nlteabund, marker="o", label="1D NLTE Ti II", color="red", linestyle="None",
                  markerfacecolor='none')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')

    print("\n")
    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d) / sum(weightti2)

    y = 0.35
    if star_number == 3:
        y += 0.5
    if star_number == 4:
        y += 0.5
    plt.figtext(0.16, 0.1, "Mean Ti I 1d: " + str(np.round(meantiI1d, 8)))
    plt.figtext(0.16, 0.2, "Mean Ti II 1d: " + str(np.round(meantiII1d, 8)))
    plt.figtext(0.16, 0.3, "Ti I/II diff in 1d: " + str(np.round(meantiII1d - meantiI1d, 8)))

    plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")


    plt.title("NLTE abundance vs equivalent width for " + star)
    lns = aa + bb
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(
        r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\" + star + "\Abundances\\" + k + "Equivalent Width NLTE")

    plt.show()













    # mean ews
    mean_ew = []
    # Get the mean equiv widths for ti1 and 2
    for r in range(len(lte_equiv_width_dict[list(ti2nlte_equiv_width_dict)[0]])):
        ews = []
        for abund in lte_equiv_width_dict:
            ew = 10**lte_equiv_width_dict[abund][r]
            ews.append(np.log10(ew/ti1lines[r]))
        mean_ew.append(sum(ews) / len(lte_equiv_width_dict))

    ews = []
    # mean ews
    ti2mean_ew = []
    for r in range(len(ti2lte_equiv_width_dict[list(ti2nlte_equiv_width_dict)[0]])):
        ews = []
        for abund in ti2lte_equiv_width_dict:
            ew = 10**ti2lte_equiv_width_dict[abund][r]
            ews.append(np.log10(ew/ti2lines[r]))
        ti2mean_ew.append(sum(ews) / len(ti2lte_equiv_width_dict))






    remv = (np.where(np.asarray(mean_ew) > -3.65))
    print(remv[0])
    print("lines to remove for arcturus", np.asarray(lines)[remv])
    list1, list2 = zip(*sorted(zip(mean_ew, lines)))
    print(list1, "\n", list2)




    # same but for lte now
    fig = plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("$W_{λ,red}$ (Ti I)")
    ax.set_ylabel("$A(Ti)$")
    ax2.set_xlabel("$W_{λ,red}$ (Ti II)")

    aa = ax.plot(mean_ew, ti1lteabund, marker='s', label="1D LTE Ti I", color="blue", markerfacecolor='none',
                 linestyle="None")

    bb = ax2.plot(ti2mean_ew, ti2lteabund, marker='o', label="1D LTE Ti II", color="red", markerfacecolor='none',
                  linestyle="None")
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')

    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d) / sum(weightti2)

    y = 0.35
    if star_number == 3:
        y += 0.5
    if star_number == 4:
        y += 0.5

    plt.figtext(0.46, y, "Mean Ti I 1d: " + str(np.round(meantiI1d, 8)))
    plt.figtext(0.46, y - 0.03, "Mean Ti II 1d: " + str(np.round(meantiII1d, 8)))
    plt.figtext(0.46, y - 0.06, "Ti I/II diff in 1d: " + str(np.round(meantiII1d - meantiI1d, 8)))
    plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")

    #plt.title("LTE abundance vs equivalent width for " + star)
    lns = aa + bb
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(
        r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\" + star + "\Abundances\\" + k + "Equivalent Width LTE")

    plt.show()


















def compare_ews_abundances():


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
    ewdict = pickle.load(open(
        r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "TI1EQUIVWIDTHS", "rb"))
    ti2ewdict = pickle.load(open(
        r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\" + star + "\\" + k + "TI2EQUIVWIDTHS", "rb"))




    # collecting ews to mean them
    # mean ews
    ewdict = 10**np.asarray(ewdict)
    ti2ewdict = 10**np.asarray(ti2ewdict)

    mean_ew = []
    saturation_limit = -4.8
    # Get the rediced equiv widths for ti1 and 2
    # list for indexes to remove
    wherelist = []
    for r in range(len(ewdict)):
        ew = ewdict[r]
        #remove oversaturated lines
        if np.log10(ew/ti1lines[r]) > saturation_limit:
            where = np.where(ewdict == ew)[0][0]
            print(np.log10(ew / ti1lines[r]))

            print("here", where)
            wherelist.append(where)
        else:
            mean_ew.append(np.log10(ew / ti1lines[r]))

    print("Ti1 Where")
    print(ti1lines)
    print("Where", wherelist)

    for where in sorted(wherelist, reverse=True):
        np.delete(ewdict, where)
        del (ti1lines[where])
        del (ti1lteabund[where])
        del (ti1nlteabund[where])

    ti2mean_ew = []
    wherelist = []

    for r in range(len(ti2ewdict)):
        ew = ti2ewdict[r]

        if np.log10(ew/ti2lines[r]) > saturation_limit:
            where = np.where(ti2ewdict == ew)[0][0]
            wherelist.append(where)
        else:
            ti2mean_ew.append(np.log10(ew / ti2lines[r]))

    for where in sorted(wherelist, reverse=True):
        np.delete(ti2ewdict, where)
        del (ti2lines[where])
        del (ti2lteabund[where])
        del (ti2nlteabund[where])
    print("Ti2 Where")
    print(ti1lines)
    print("Where", wherelist)




    if star == "NOTUSEDsun":
        weightti1 = ti1lte['weight']
        weightti2 = ti2lte['weight']
    else:
        weightti1 = np.ones(len(ti1lines))
        weightti2 = np.ones(len(ti2lines))

    """for x in range(len(ti2lines)):
        print(ti2lines[x])
        print(ti2lteabund[x])
        print(ti2lte_threeDdiff[x])
        print("\n")
    exit()"""

    fig = plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax.set_xlabel("$W(_{λ,red}$)")
    ax.set_ylabel("A(Ti)")

    plt.plot(mean_ew, ti1nlteabund, marker='s', label="1D NLTE Ti I", color="blue", linestyle="None",
                 markerfacecolor='none')

    plt.plot(ti2mean_ew, ti2nlteabund, marker="o", label="1D NLTE Ti II", color="red", linestyle="None",
                  markerfacecolor='none')

    print("\n")
    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d) / sum(weightti2)

    y = 0.35
    if star_number == 3:
        y += 0.005
    if star_number == 4:
        y += 0.005
    """plt.figtext(0.26, y+0.5, "Mean NLTE Ti I 1D: " + str(np.round(meantiI1d, 8)))
    plt.figtext(0.26, y +0.45, "Mean NLTE Ti II 1D: " + str(np.round(meantiII1d, 8)))
    plt.figtext(0.26, y +0.4, "NLTE Ti I/II diff in 1D: " + str(np.round(meantiII1d - meantiI1d, 8)))"""

    #plt.axhline(np.round(meantiI1d, 8), color="b", ls='--')

    #plt.axhline(np.round(meantiII1d, 8), color="r", ls='--')






    import statistics
    ti1nlteabund = [float(x) for x in ti1nlteabund]
    ti2nlteabund = [float(x) for x in ti2nlteabund]











    if star == "arcturusKarin":
        remv = (np.where(np.asarray(mean_ew) > -4.57))
        print(remv[0])
        print("lines to remove for arcturus", np.asarray(lines)[remv])
    else:
        remv = (np.where(np.asarray(ti1lteabund) > 3.2))
        print(remv[0])
        print("lines to remove for arcturus", np.asarray(lines)[remv])

    list1, list2 = zip(*sorted(zip(mean_ew, lines)))




    """aa2 = ax.plot(mean_ew, ti1lteabund, marker='s', label="1D LTE Ti I", color="blue",
                 linestyle="None")

    bb2 = ax2.plot(ti2mean_ew, ti2lteabund, marker='o', label="1D LTE Ti II", color="red",
                  linestyle="None")    """


    plt.plot(mean_ew, ti1lteabund, marker='s', label="1D LTE Ti I", color="blue",
                 linestyle="None")

    plt.plot(ti2mean_ew, ti2lteabund, marker='o', label="1D LTE Ti II", color="red",
                  linestyle="None")




    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_1d = [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d) / sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d) / sum(weightti2)

    y = 0.35
    if star_number == 3:
        y += 0.5
    if star_number == 4:
        y += 0.5

    """plt.figtext(0.46, y, "Mean LTE Ti I 1d: " + str(np.round(meantiI1d, 8)))
    plt.figtext(0.46, y - 0.03, "Mean LTE Ti II 1d: " + str(np.round(meantiII1d, 8)))
    plt.figtext(0.46, y - 0.06, "LTE Ti I/II diff in 1d: " + str(np.round(meantiII1d - meantiI1d, 8)))"""
    """plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")"""


    minx = plt.xlim()[0]
    maxx = plt.xlim()[1]
    z = np.polyfit(mean_ew, ti1lteabund, 1)
    p = np.poly1d(z)
    yval = p(np.arange(minx,maxx, 0.01))
    plt.plot(np.arange(minx,maxx, 0.01), yval, color="b")


    z = np.polyfit(ti2mean_ew, ti2lteabund, 1)
    p = np.poly1d(z)
    yval = p(np.arange(minx,maxx, 0.01))
    plt.plot(np.arange(minx,maxx, 0.01), yval, color="r")

    z = np.polyfit(mean_ew, ti1nlteabund, 1)
    p = np.poly1d(z)
    yval = p(np.arange(minx,maxx, 0.01))
    plt.plot(np.arange(minx,maxx, 0.01), yval, color="b", ls="--")


    z = np.polyfit(ti2mean_ew, ti2nlteabund, 1)
    p = np.poly1d(z)
    yval = p(np.arange(minx,maxx, 0.01))
    plt.plot(np.arange(minx,maxx, 0.01), yval, color="r", ls="--")


    #plt.title("LTE abundance vs equivalent width for " + star)
    """    lns = aa + bb + aa2+bb2
        labs = [l.get_label() for l in lns]
    """

    plt.legend()
    if star == "sun":
        ax.legend(loc="upper left")
    if star == "sun":
        startitle = "Sun"
    elif star_number == 18:
        startitle = "122563 Modified"
    elif star == "arcturusKarin":
        startitle = "Arcturus"
    else:
        startitle ="HD"+star
    plt.title(startitle)
    lims = plt.xlim()
    plt.ylim([meantiII1d-0.8, meantiII1d+0.6])

    plt.savefig(
        r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\all\\"+star)


    import statistics
    ti1lteabund = [float(x) for x in ti1lteabund]
    ti2lteabund = [float(x) for x in ti2lteabund]

    print("Standard Dev. of Ti I LTE for", star, " is", statistics.stdev(ti1lteabund))
    print("Standard Dev. of Ti II LTE for", star, " is", statistics.stdev(ti2lteabund))

    print("Standard Dev. of Ti I LTE for", star, " is", statistics.stdev(ti1nlteabund))
    print("Standard Dev. of Ti II LTE for", star, " is", statistics.stdev(ti2nlteabund))

    print("Error of Ti I NLTE for", star, " is", statistics.stdev(ti1lteabund)/np.sqrt(len(ti1lteabund)))
    print("Error of Ti II NLTE for", star, " is", statistics.stdev(ti2lteabund)/np.sqrt(len(ti2lteabund)))

    print("Error of Ti I LTE for", star, " is", statistics.stdev(ti1nlteabund)/np.sqrt(len(ti1nlteabund)))
    print("Error of Ti II LTE for", star, " is", statistics.stdev(ti2nlteabund)/np.sqrt(len(ti2nlteabund)))

    print("Error of NLTE ionisation imbalance for", star, " is", np.sqrt(
        statistics.stdev(ti1nlteabund)/np.sqrt(len(ti1nlteabund))**2 + (statistics.stdev(ti2nlteabund)/np.sqrt(len(ti2nlteabund)))**2))
    print("Error of LTE ionisation imbalance for", star, " is", np.sqrt(
        statistics.stdev(ti1lteabund)/np.sqrt(len(ti1lteabund))**2 + (statistics.stdev(ti2lteabund)/np.sqrt(len(ti2lteabund)))**2))
    print("Number of Ti I:", len(ti1nlteabund))
    print("Number of Ti II:", len(ti2nlteabund))

    plt.show()




# compare the abundances vs scott per line for nlte and lte, for ti II or ti I
def compare_scott():
    nltediff = [nlte_val - scott_val for scott_val, nlte_val in zip(literature_abund, nlteabund)]
    ltediff = [lte_val - scott_val for scott_val, lte_val in zip(literature_abund, lteabund)]
    print(nlteabund)
    print(np.mean(nlteabund))
    print(nlte_literature_abund)
    print(np.mean(nlte_literature_abund))
    print("\n")
    print(lteabund)
    print(np.mean(lteabund))
    print(literature_abund)
    print(np.mean(literature_abund))
    plt.scatter(lines, nlteabund, label = "NLTE difference from Literature", edgecolors="red", facecolors='none')
    plt.scatter(lines, lteabund, marker= "s",  label = "LTE difference from Literature", edgecolors = "blue", facecolors='none')
    plt.scatter(lines, literature_abund, marker= "s", label="Scott", edgecolors="black", facecolors = "none")
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("A(Ti I)$_{Mallinson}$ - A(Ti I)$_{Scott}$")
    if ti2:
        plt.title("Abundance calculations in N/LTE for Ti II for star "+str(star))
    if not ti2:
        plt.title("Abundance calculations in N/LTE for Ti I for star "+str(star))
    plt.legend()

    x = 0.6
    if ti2:
        y = 0.34
    else:
        y = 0.64

    if star_number==3:
        if ti2:
            y+=0.7
        y-=0.22
        x-=0.12
    if star_number == 4:
        if not ti2:
            y+=0.2
            x -=0.4
    if star_number == 5:
        if ti2:
            y+=0.65
        else:
            y = 0.9
            x-= 0.1
    print(str(np.round(np.mean(nltediff), 8)))
    #plt.figtext(x, y-0.2, "NLTE mean: " + str(np.round(np.mean(nltediff), 8)))
    #plt.figtext(x, y-0.16, "LTE mean: " + str(np.round(np.mean(ltediff), 8)))

    if ti2:
        plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"LiteratureAbundTi2_differences_nlte_lte")

    elif not ti2:

        plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"LiteratureAbund_differences_nlte_lte")


    plt.show()



# Plots the differences for Ti I or II abundances when calculated in NLTE vs LTE.
def lte_nlte_differences():
    if ti2:
        ti = "Ti II"
    else:
        ti = "Ti I"

    dexdiff = [lteval - nlteval for nlteval, lteval in zip(nlteabund, lteabund)]

    plt.scatter(lines, dexdiff, marker="s", label = "N/LTE abundance difference for "+ti, color = "blue", facecolor='none')


    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Difference in abundance (dex)")

    x = 0.55
    if ti2:
        y=0.25
    else:
        y=0.75
    if star_number==3:
        if ti2: y+=0.5
        y-=0.55
    if star_number == 4:
        x-=0.4
        if ti2:
            y+=0.57
    if star_number == 5:
        if not ti2:
            x-= 0.42
    plt.figtext(x, y, "Mean difference: " + str(np.round(np.mean(dexdiff), 8)))


    plt.title("Difference in N/LTE calculations "+ti+" for star "+star)

    plt.legend()

    if ti2:
        plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"nlte_lte_differences_Ti2")

    elif not ti2:

        plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"nlte_lte_differences_Ti1")


    plt.show()

# testing by adding the 3d difference of LTE from Literature to our results to see the results. comparing nlte vs lte
def test_3d():
    if star != "sun":
        print("Not the sun, no 3D differences to be added.")
        return

    return


    ti2lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTEti2_abundance_comparisons", "rb"))
    ti2nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTEti2_abundance_comparisons", "rb"))
    ti1lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTE_abundance_comparisons", "rb"))
    ti1nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTE_abundance_comparisons", "rb"))

    ti1lines = ti1lte['lines']
    ti1lteabund = ti1lte['our_abund']

    if star == "sun":
        ti1lte_threeDdiff = ti1lte['3d_diff']
        ti1nlte_threeDdiff = ti1nlte['3d_diff']
        ti2lte_threeDdiff = ti2lte['3d_diff']
        ti2nlte_threeDdiff = ti2nlte['3d_diff']

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
    """for x in range(len(ti2lines)):
        print(ti2lines[x])
        print(ti2lteabund[x])
        print(ti2lte_threeDdiff[x])
        print("\n")
    exit()"""



    # added with the 3d differences
    ti1nltediff = [nlte_val+diff for diff, nlte_val in zip(ti1nlte_threeDdiff, ti1nlteabund)]
    ti1ltediff = [lte_val+diff for diff, lte_val in zip(ti1lte_threeDdiff, ti1lteabund)]

    ti2nltediff = [nlte_val+diff for diff, nlte_val in zip(ti2nlte_threeDdiff, ti2nlteabund)]
    ti2ltediff = [lte_val+diff for diff, lte_val in zip(ti2lte_threeDdiff, ti2lteabund)]
    print(ti1nlte_threeDdiff)
    print(ti2nlte_threeDdiff)


    fig=plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("Wavelength for Ti I ($\AA$)")
    ax.set_ylabel("$Log_{10}$(Abundance)")
    ax2.set_xlabel("Wavelength for Ti II ($\AA$)")

    a = ax.plot(ti1lines, ti1nltediff, label = "3D NLTE Ti I", linestyle="dashed", color="red", edgecolors="red", facecolors='none')
    aa= ax.plot(ti1lines, ti1nlteabund, label = "1D NLTE Ti I", color="red", edgecolors="red", facecolors='none')

    b = ax2.plot(ti2lines, ti2nltediff, label = "3D NLTE Ti II", linestyle="dashed" , color="blue", edgecolors="red", facecolors='none')
    bb = ax2.plot(ti2lines, ti2nlteabund, label = "1D NLTE Ti II", color="blue", edgecolors="red", facecolors='none')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')


    meantiII3d = np.mean(ti2nltediff)
    meantiII1d = np.mean(ti2nlteabund)

    meantiI3d = np.mean(ti1nltediff)
    meantiI1d = np.mean(ti1nlteabund)

    print("\n")
    weighted_ti2_3d = [threeDval * weight for threeDval, weight in zip(ti2nltediff, weightti2)]
    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_3d = [threeDval * weight for threeDval, weight in zip(ti1nltediff, weightti1)]
    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]


    # Applying weightings. Remove if not wanted to be weighted.
    meantiI3d = sum(weighted_ti1_3d)/sum(weightti1)
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII3d = sum(weighted_ti2_3d)/sum(weightti2)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)


    plt.figtext(0.5, 0.4, "Mean Ti I 3d: " + str(np.round(meantiI3d, 8)))
    plt.figtext(0.5, 0.35, "Mean Ti I 1d: " + str(np.round(meantiI1d, 8)))

    plt.figtext(0.5, 0.3, "Mean Ti II 3d: " + str(np.round(meantiII3d, 8)))
    plt.figtext(0.5, 0.25, "Mean Ti II 1d: " + str(np.round(meantiII1d, 8)))


    plt.figtext(0.5, 0.17, "Ti I/II diff in 3d: " + str(np.round(meantiII3d-meantiI3d, 8)))
    plt.figtext(0.5, 0.12, "Ti I/II diff in 1d: " + str(np.round(meantiII1d-meantiI1d, 8)))

    plt.title("3d differences of Ti I/II in NLTE for star "+star)
    lns = a + aa + bb + b
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"3d_differences_NLTE")

    plt.show()




    # same but for lte now
    fig=plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("Wavelength for Ti I ($\AA$)")
    ax.set_ylabel("$Log_{10}$(Abundance)")
    ax2.set_xlabel("Wavelength for Ti II ($\AA$)")

    a = ax.plot(ti1lines, ti1nltediff, label = "3D LTE Ti I", linestyle="dashed", color="red", edgecolors="red", facecolors='none')
    aa= ax.plot(ti1lines, ti1lteabund, label = "1D LTE Ti I", color="red", edgecolors="red", facecolors='none')

    b = ax2.plot(ti2lines, ti2ltediff, label = "3D LTE Ti II", linestyle="dashed" , color="blue", edgecolors="red", facecolors='none')
    bb = ax2.plot(ti2lines, ti2lteabund, label = "1D LTE Ti II", color="blue", edgecolors="red", facecolors='none')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')

    meantiII3d = np.mean(ti2ltediff)
    meantiII1d = np.mean(ti2lteabund)

    meantiI3d = np.mean(ti1ltediff)
    meantiI1d = np.mean(ti1lteabund)


    weighted_ti2_3d = [threeDval * weight for threeDval, weight in zip(ti2ltediff, weightti2)]
    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_3d = [threeDval * weight for threeDval, weight in zip(ti1ltediff, weightti1)]
    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI3d = sum(weighted_ti1_3d)/sum(weightti1)
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII3d = sum(weighted_ti2_3d)/sum(weightti2)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)


    plt.figtext(0.5, 0.4, "Mean Ti I 3d: " + str(np.round(meantiI3d, 8)))
    plt.figtext(0.5, 0.35, "Mean Ti I 1d: " + str(np.round(meantiI1d, 8)))

    plt.figtext(0.5, 0.3, "Mean Ti II 3d: " + str(np.round(meantiII3d, 8)))
    plt.figtext(0.5, 0.25, "Mean Ti II 1d: " + str(np.round(meantiII1d, 8)))



    # we apply weightings (so if a weight of 2, we double the value)

    plt.figtext(0.5, 0.17, "Ti I/II diff in 3d: " + str(np.round(meantiII3d-meantiI3d, 8)))
    plt.figtext(0.5, 0.12, "Ti I/II diff in 1d: " + str(np.round(meantiII1d-meantiI1d, 8)))

    plt.title("3d differences of Ti I/II in LTE for star"+star)
    lns = a + aa + bb + b
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"3d_differences_LTE")

    plt.show()


# finding ionization imbalance of 1d N/LTE (As 3D was just a bit of a hack.
def ionization_imbalance():

    ti2lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTEti2_abundance_comparisons", "rb"))
    ti2nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTEti2_abundance_comparisons", "rb"))
    ti1lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTE_abundance_comparisons", "rb"))
    ti1nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTE_abundance_comparisons", "rb"))

    ti1lines = ti1lte['lines']
    ti1lteabund = ti1lte['our_abund']
    ti1nlteabund = ti1nlte['our_abund']

    ti2lines = ti2lte['lines']
    ti2lteabund = ti2lte['our_abund']
    ti2nlteabund = ti2nlte['our_abund']

    # removed Ti II   3348.84   0.12  -1.18   55.4  -4.78   2.79
    # Ti II   3349.03   0.61   0.46   86.1  -4.59   2.64
    # Ti II   3349.40   0.05   0.54  124.0  -4.43   2.81
    # Ti II   3477.18   0.12  -0.95   65.7  -4.72   2.86
    # Ti II   4394.06   1.22  -1.77    9.7  -5.66   2.89
    """ti2lines.extend([3348.84, 3349.03, 3349.40, 3477.18, 4394.06])
    ti2lteabund.extend([2.79, 2.64, 2.81, 2.86, 2.89])
    ti2nlteabund.extend([2.79, 2.64, 2.81, 2.86, 2.89])"""

    if star =="sun":
        weightti1 = ti1lte['weight']
        weightti2 = ti2lte['weight']
    else:
        weightti1 = np.ones(len(ti1lines))
        weightti2 = np.ones(len(ti2lines))

    """for x in range(len(ti2lines)):
        print(ti2lines[x])
        print(ti2lteabund[x])
        print(ti2lte_threeDdiff[x])
        print("\n")
    exit()"""





    fig=plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("Wavelength for Ti I ($\AA$)")
    ax.set_ylabel("$Log_{10}$(Abundance)")
    ax2.set_xlabel("Wavelength for Ti II ($\AA$)")

    import statistics
    ti1nlteabund = [float(x) for x in ti1nlteabund]
    ti2nlteabund = [float(x) for x in ti2nlteabund]

    print("Standard Dev. of Ti I NLTE for", star, " is", statistics.stdev(ti1nlteabund))
    print("Standard Dev. of Ti II NLTE for", star, " is", statistics.stdev(ti2nlteabund))

    aa= ax.plot(ti1lines, ti1nlteabund, marker='s', label = "1D NLTE Ti I", color="blue", linestyle="None", markerfacecolor='none')

    bb = ax2.plot(ti2lines, ti2nlteabund, marker="o", label = "1D NLTE Ti II", color="red", linestyle="None", markerfacecolor='none')
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')



    print("\n")
    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]


    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)
    plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")

    y = 0.35
    if star_number == 3:
        y+=0.5
    if star_number == 4:
        y+=0.5
    plt.figtext(0.46, y, "Mean Ti I: " + str(np.round(meantiI1d, 5)))
    plt.figtext(0.46, y-0.03, "Mean Ti II: " + str(np.round(meantiII1d, 5)))


    plt.figtext(0.46, y-0.06, "Ti I-II: " + str(np.round(meantiI1d-meantiII1d, 5)))

    plt.title("Ionization imbalance of Ti I/II in NLTE for star "+star)
    lns = aa + bb
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"Ionization_imbalance_NLTE")

    plt.show()




    # same but for lte now
    fig=plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("Wavelength for Ti I ($\AA$)")
    ax.set_ylabel("$Log_{10}$(Abundance)")
    ax2.set_xlabel("Wavelength for Ti II ($\AA$)")
    ti1lteabund = [float(x) for x in ti1lteabund]
    ti2lteabund = [float(x) for x in ti2lteabund]

    print("Standard Dev. of Ti I LTE is", statistics.stdev(ti1lteabund))
    print("Standard Dev. of Ti II LTE is", statistics.stdev(ti2lteabund))

    aa= ax.plot(ti1lines, ti1lteabund, marker='s',  label = "1D LTE Ti I", color="blue", markerfacecolor='none', linestyle="None")

    bb = ax2.plot(ti2lines, ti2lteabund,  marker = 'o',label = "1D LTE Ti II", color="red", markerfacecolor='none', linestyle="None")
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')



    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)




    y = 0.35
    if star_number == 3:
        y+=0.5
    if star_number == 4:
        y+=0.5

    plt.figtext(0.46, y, "Mean Ti I: " + str(np.round(meantiI1d, 5)))

    plt.figtext(0.46, y-0.03, "Mean Ti II: " + str(np.round(meantiII1d, 5)))


    plt.figtext(0.46, y-0.06, "Ti I-II: " + str(np.round(meantiI1d-meantiII1d, 5)))
    plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")




    plt.title("Ionization imbalance of Ti I/II in LTE for star "+star)
    lns = aa + bb
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)
    plt.savefig(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\graphs\\"+star+"\Abundances\\"+k+"Ionization_imbalance_LTE")

    plt.show()


# For a thing for the overleaf file
def all_stars():





    ti2lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTEti2_abundance_comparisons", "rb"))
    ti2nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTEti2_abundance_comparisons", "rb"))
    ti1lte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"LTE_abundance_comparisons", "rb"))
    ti1nlte = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\\"+star+"\\" + k+"NLTE_abundance_comparisons", "rb"))

    ti1lines = ti1lte['lines']
    ti1lteabund = ti1lte['our_abund']
    ti1nlteabund = ti1nlte['our_abund']

    ti2lines = ti2lte['lines']
    ti2lteabund = ti2lte['our_abund']
    ti2nlteabund = ti2nlte['our_abund']

    if star =="sun":
        weightti1 = ti1lte['weight']
        weightti2 = ti2lte['weight']
    else:
        weightti1 = np.ones(len(ti1lines))
        weightti2 = np.ones(len(ti2lines))

    """for x in range(len(ti2lines)):
        print(ti2lines[x])
        print(ti2lteabund[x])
        print(ti2lte_threeDdiff[x])
        print("\n")
    exit()"""





    fig=plt.figure(figsize=(8, 6), dpi=100)

    ax = fig.add_subplot(111, label="1")
    ax2 = fig.add_subplot(111, label="2", frame_on=False, sharey=ax)
    ax.set_xlabel("Wavelength for Ti I ($\AA$)")
    ax.set_ylabel("$A(Ti)$")
    ax2.set_xlabel("Wavelength for Ti II ($\AA$)")

    import statistics
    ti1nlteabund = [float(x) for x in ti1nlteabund]
    ti2nlteabund = [float(x) for x in ti2nlteabund]

    print("Standard Dev. of Ti I NLTE is", statistics.stdev(ti1nlteabund))
    print("Standard Dev. of Ti II NLTE is", statistics.stdev(ti2nlteabund))

    aa= ax.plot(ti1lines, ti1nlteabund, marker='s', label = "1D NLTE Ti I", color="blue", linestyle="None", markerfacecolor='none')

    bb = ax2.plot(ti2lines, ti2nlteabund, marker="o", label = "1D NLTE Ti II", color="red", linestyle="None", markerfacecolor='none')
    aa1= ax.plot(ti1lines, ti1lteabund, marker='s',  label = "1D LTE Ti I", color="blue", linestyle="None")

    bb1 = ax2.plot(ti2lines, ti2lteabund,  marker = 'o',label = "1D LTE Ti II", color="red", linestyle="None")

    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')

    y = 0.35
    if star_number == 3:
        y+=0.5
    if star_number == 4:
        y+=0.5

    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2nlteabund, weightti2)]

    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1nlteabund, weightti1)]


    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)
    plt.axhline(np.round(meantiI1d, 8), ls='--', color="b")

    plt.axhline(np.round(meantiII1d, 8),  ls='--', color="r")




    ti1lteabund = [float(x) for x in ti1lteabund]
    ti2lteabund = [float(x) for x in ti2lteabund]

    print("Standard Dev. of Ti I LTE is", statistics.stdev(ti1lteabund))
    print("Standard Dev. of Ti II LTE is", statistics.stdev(ti2lteabund))


    weighted_ti2_1d = [threeDval * weight for threeDval, weight in zip(ti2lteabund, weightti2)]

    weighted_ti1_1d =  [threeDval * weight for threeDval, weight in zip(ti1lteabund, weightti1)]

    # Applying weightings. Remove if not wanted to be weighted.
    meantiI1d = sum(weighted_ti1_1d)/sum(weightti1)
    meantiII1d = sum(weighted_ti2_1d)/sum(weightti2)



    plt.axhline(np.round(meantiI1d, 8), color="b")

    plt.axhline(np.round(meantiII1d, 8), color="r")



    #plt.figtext(0.46, y-0.06, "Ti I/II diff in 1d: " + str(np.round(meantiII1d-meantiI1d, 8)))

    #plt.title("Ionization imbalance of Ti I/II for star "+star)
    """aa2 = plt.plot(ti1lines[0], ti1nlteabund[0], ls = "--", label = "NLTE mean", color = 'blue')
    aa3 = plt.plot(ti1lines[0], ti1nlteabund[0], ls = "--", label = "LTE mean", color = 'red')"""
    lns = aa + bb + bb1 + aa1
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)


    print(star)
    if star == "sun":
        star_name = "Sun"
    elif star =="arcturusKarin":
        star_name = "Arcturus"
    else:
        star_name = str(star)

    plt.title(star_name)
    lns = aa + bb
    labs = [l.get_label() for l in lns]
    ax.legend(lns, labs, loc=0)

    plt.show()







"""compare_ews()
compare_scott()
lte_nlte_differences()"""


compare_ews_abundances()