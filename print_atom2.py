# Putting everything in format that balder requires. Remember balder starts lists at 1, not 0 like python.
# So that's annoying. This proints the atom2, with unmerged lines. o there are more levels, but fewer BB transitions
# as we only trake the experimentally confirmed ones.

import pickle
import numpy as np
import  matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
atomname = r"D:\PhD\TiFeAtoms\Jack\balder\atom2\atom.Ti782.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"

atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data



# prints to file
import sys

# for atom2, instead we find the levels with experimental results
def find_experimental():
    top_indexes = []
    for ref in range(len(atom2['ref'])):
        if "K16" not in atom2['ref'][ref]:
            top_indexes.append(ref)
    return (top_indexes)



# Getting the indexes of the top lines (For grid_points), the top lines themselves and the rates/data they have,
# and the list of equiv widths.
top_indexes = find_experimental()


indexes = top_indexes

# print to notepad
orig_stdout = sys.stdout
f = open('atom2.ti1', 'w')
sys.stdout = f

rtc = 1.0973731569E5  # Conversion factor between Rydberg and cm^(-1)
rte = 13.605698  # Conversion factor between Rydberg and eV

tlist = ["1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "10000", "12000", "20000"]




electrons = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\Final_outputs\full_electron_rate_dict_all_temperatures.pkl", "rb"))
radiative = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\Final_outputs\full_radiative_xsection_nahar+regemorter.pkl", "rb"))
hydrogen = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\Final_outputs\Full_combined_hydrogen_rates.pkl", "rb"))

# now we start printing to the file
print("TI")

        #Abundance  unused                  # AU
print(" 5.04      47.867 ")
# number of levels, number of bb radiuative transitions, number of bf radiative transitions (Not the number of photons)
#print("    ", len(atom1['ec']), "       ", len(atom2['irad']), "       ", len(radiative), "        ", 0)
print("    ", len(atom1['ec']), "       ", len(indexes), "       ", len(radiative), "        ", 0)

print("***\n*** Levels\n***")
# Levels
# For atom2 there will be uncollapsed levels that need identical indexes/notes
for level in range(len(atom1['ec'])):
                                                                                                                                                        # J is unimportant so I just put 0
        print("   ", f"""{atom1['ec'][level]:<22} {atom1['g'][level]:<10}  {"'" + str(atom1['term'][level]) + ":      " + str(atom1['conf'][level]) + "'":<35s} {atom1['ion'][level]:5} {level+1:5} {0:5}""")

# Lines, BB radiative
print("***\n*** Lines\n***")

x = 0
#for bb_rad_transition in range(len(atom2['irad'])):             # to convert it to AA from nm we divide by 10
# In indexes, the indexes of the top x number of transitions.
for bb_rad_transition in indexes:  # to convert it to AA from nm we divide by 10

    if int(atom2['jrad'][bb_rad_transition]) == 12 or int(atom2['jrad'][bb_rad_transition]) == 66:
        if int(atom2['irad'][bb_rad_transition]) == 66 or int(atom2['irad'][bb_rad_transition])== 12:
            pass
        else:
            continue

    else:
        continue
    # Print at 1 level above as balder starts arrays at 1 for some reason
    print(f"{'***':5} {str(bb_rad_transition+1)+':':<5}  {atom2['ALAMB'][bb_rad_transition]/10:11}", "nm")
                                                                                                                                # 0 for atom 1, until we start simplifying stuff. After it is 1 which is fine structure components. Ask karin. IW?#
                                                                                                                                # ga stuff is "broadnening" in anish's, but "dampening" in my email, still fine?
    '''    print("   ", f"{atom2['jrad'][bb_rad_transition]:<10} {atom2['irad'][bb_rad_transition]:<10} {atom2['f'][bb_rad_transition]:<26} "
                 f"{atom2['nq'][bb_rad_transition]:<5} {atom2['qmax'][bb_rad_transition]:<5} {0:<3} {1:<3} "
                 f"{atom2['ga'][bb_rad_transition]:>10} {atom2['gw'][bb_rad_transition]:>5} {atom2['gq'][bb_rad_transition]:>5}")
    '''
    # reducing nq to 1
    print("   ", f"{atom2['jrad'][bb_rad_transition]:<10} {atom2['irad'][bb_rad_transition]:<10} {atom2['f'][bb_rad_transition]:<26} "
                     f"{atom2['nq'][bb_rad_transition]:<5} {atom2['qmax'][bb_rad_transition]} {0:<3} {1:<3} " # nq and qmax and 0 and 1
                     f"{atom2['ga'][bb_rad_transition]:>10} {atom2['gw'][bb_rad_transition]:>5} {atom2['gq'][bb_rad_transition]:>5}")
    print("***")
    x += 1

exit()
# Continua, bf radiative transitions. Remember we must increase the printed levels by 1 to account for Balder starting
# arrays at 1

print("***\n*** Continua\n***")


# Real code for printing all interpolated levels
for bf_rad_transition in range(len(radiative)):
    # Loading interpolated cross sections, in order of original atom1 transitions. Each file is a different level to level transition, consisting of many different wavelengths on a set grid.
    xsections = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\interpolation\interpolated_cross_sections/cross_section_interpolated_transition"+str(bf_rad_transition)+".pkl", "rb"))
    # Loading the wavelengths that gave rise to the interpolated cross sections
    wavelengths = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\balder\interpolation\interpolated_grids\wavelength_grid_interpolated_transition"+str(bf_rad_transition)+".pkl", "rb"))
    # Is still in the same order as the original radiative, so all other information is correct (Target, e.g)
    # Priting transition + 1 as we don't start at 0 in balder.
    print(" ", f"{int(radiative[bf_rad_transition]['target'])+1:< 7} {bf_rad_transition+1:<7} {xsections[0]:<25} {int(np.ceil(len(xsections))):<10}", " -1")
    # So for each photon wavelength/energy, we have an individual cross section.
    for photon in range(0, len(radiative[bf_rad_transition]['wavelength'])):
        print(" ", f"{wavelengths[photon]:<20} {xsections[photon]:<15}")


# Collisions
print("***\nNEWCOL\n***\nTEMP")
print(len(tlist), end="  ")
for temp in tlist:
    print(temp, end="  ")
print("\n***")

# Electron dict doesn't have excitation BB so no need to check for those. and is arranged [icol][jcol][temp] so no
# need to ru nthrough temps like in Hydrogen.
for x in electrons.keys():
    for y in electrons[x]:

        # These are all electron collisions in this dictionary
        if len(electrons[x][y]) == 0:
            continue

        # CI for ionisation collisions
        if x < 459 and y >= 459:
            print("CI")
        else:
            print("CE")
        print("  ", f"{int(x)+1:<4} {y+1:<4}", end = "   ")
        for temperature in electrons[x][y].keys():
            print(f"{electrons[x][y][temperature]:<10}", end = "   ")

        print("")

#' temperature is unimportant here, we just loop through aech key. BB transitions for Ti 1
for x in hydrogen["1000"]:
    for y in hydrogen["1000"][x]:

        # Checking if any rates are above 0. Different method to electrons as the dictionary for hydrogen is structured
        # by temp first.
        no_rates = True
        for temperature in hydrogen.keys():
            if hydrogen[temperature][x][y] > 0:
                no_rates = False
                break

        if no_rates: continue

        # Electrons didn't have any BB excitation saved, only BB de-ex and BF excit so it was unecessary there.
        # But here we can't have excitations during BB transitions.
        if y > x:
            if y < 459 and x < 459:
                continue
            elif y >= 459 and x >= 459:
                continue



        # If it's a charge transfer transition it's ch0. None occur above the level 223 though due to the
        # hydrogenic potential limit of -0.754eV
        if y >= 459 and x < 459:
            print("CH0")
        else:
            print("CH")
        print("  ", f"{int(x) + 1:<4} {y + 1:<4}", end="   ")

        for temperature in hydrogen.keys():
            print(f"{hydrogen[temperature][x][y]:<10}", end="   ")

        print("")

# A file for CT transitions from ti I to ti II
ti2 = pickle.load(open(r"D:\PhD\TiFeAtoms\Jack\ti2grumerrates.pkl", "rb"))
# Now repreated but for CT from Ti I to tI II. Doesn't include Ti 2
for x in ti2['1000']:
    for y in ti2['1000'][x]:
        # Checking if any rates are above 0. Different method to electrons as the dictionary for hydrogen is structured
        # by temp first.
        no_rates = True

        for temperature in ti2.keys():
            if float(ti2[temperature][x][y]) > 0:
                no_rates = False
                break

        if no_rates: continue

        # No BB excitations, deexcitations only.
        if y > x:
            if y < 459 and x < 459:
                continue
            elif y >= 459 and x >= 459:
                continue
        # Excitation CT only.
        if y<x:
            if y < 459 and x >= 459:
                continue

        if y >= 459 and x < 459:
            print("CH0                       comment")
            print("  ", f"{int(x) + 1:<4} {y + 1:<4}", end="   ")
            for temperature in ti2.keys():
                print(f"{ti2[temperature][x][y]:<10}", end="   ")

        print("")


# Now we run BB for Ti II using Andrey Belyaev rates, and one Ti II - > Ti III level transition
andrey_rates = pickle.load(open("D:\PhD\TiFeAtoms\Jack\AndreyCollisionalTiII.pkl", "rb"))
tlist = ["1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"]
print("***\nTEMP")
print(len(tlist), end="  ")
for temp in tlist:
    print(temp, end="  ")
print("\n***")

for x in andrey_rates[1000]:
    for y in andrey_rates[1000][x]:
        # Checking if any rates are above 0. Different method to electrons as the dictionary for hydrogen is structured
        # by temp first.
        no_rates = True
        for temp in andrey_rates.keys():
            temperature = int(temp)
            if float(andrey_rates[temperature][x][y]) > 0:
                no_rates = False
                break

        if no_rates: continue

        # No BB excitations, deexcitations only.
        if y > x:
            # Ti III means this is CT so y>x is good
            if y == 586:
                pass
            elif y >= 459 and x >= 459:
                continue

        # No de-excitations from Ti III
        if x == 586:
            continue
        if y == 586:
            print("CH0                       comment")
        else:
            print("CH                       comment")
        print("  ", f"{int(x) + 1:<4} {y + 1:<4}", end="   ")
        for temperature in andrey_rates.keys():
            print(f"{andrey_rates[temperature][x][y]:<10}", end="   ")
        print("")

print("***\nEND")
sys.stdout = orig_stdout
f.close()