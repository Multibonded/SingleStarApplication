# Putting everything in format that balder requires. Remember balder starts lists at 1, not 0 like python.
# So that's annoying.

import pickle
import numpy as np
import  matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io import fits
atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
#atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Fe463_nq1541144_wH.fits"

atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data


# prints to file
import sys

# We find the indexes of the top transitions, for finding the important BB transitions and assigning nq accordingly.

def index_top_transitions():

    x=0
    # List for the data
    fullline = []
    # List for equiv widths only to determine who gets how many freq points
    ew_sun = []

    # Open the file and seperate lines into arrays for data crunching
    with open(r"D:\PhD\TiFeAtoms\Jack\balder\multi3d\Sun_hd84937_rb_weq.txt") as f:
        for line in f:
            # Skip title
            if x==0:
                x += 1
                continue

            x+=1

            # turns them into floats and adds to a list full of each line in full
            z = []
            # Stops when we reach some weird physical issue. Ones after this point are BF which we don't want now.
            if line.split()[0] == "************":
                break
            for y in line.split():
                z.append(float(y))

            ew_sun.append(abs(z[3]))

            fullline.append(z)

    fullline = np.asarray(fullline)


    cutoff = (4000)

    # Find the top x(cutoff) transitions for each region from multi

    # top 1000 radiative brackets
    s5 = sorted(fullline[:, 5])[-cutoff:]
    # The smallest value there
    s5value = min(s5)
    # Find index where the value of this region was greater than the minimum required for the top 1000 (Will be n>1000)
    s5index = np.where(fullline[:,5] >= s5value)

    s6 = sorted(fullline[:, 6])[-cutoff:]
    s6value = min(s6)
    s6index = np.where(fullline[:,6] >= s6value)

    s7 = sorted(fullline[:, 7])[-cutoff:]
    s7value = min(s7)
    s7index = np.where(fullline[:,7] >= s7value)

    s8 = sorted(fullline[:, 8])[-cutoff:]
    s8value = min(s8)
    s8index = np.where(fullline[:,8] >= s8value)



    # Combine them to get the most powerful transuitions in EACH region, tau 0.4, sun 1.4 etc.
    combined = np.concatenate((s8index,s7index,s6index,s5index), axis = None)
    smalllen = (len(combined))
    # Unique array index of the top 1000 of each region. Obviously different than 4000 as there are many overlapping transit,
    # that are the top in more than one region
    top_indexes = set(combined)
    print("When using the top", cutoff, "transitions, there were", smalllen, "values that made it out of the 4 regions. (Cutoff * 4)"
           "\n But when removing duplicates there were", len(top_indexes),", meaning", len(top_indexes)-cutoff, "values were not the top "
           "in all regions")

    # convert to array of indexes of top lines.
    top_indexes = (np.asarray(list(top_indexes)))
    top_lines = fullline[top_indexes]
    pickle.dump(top_indexes, open("top_indexes", "wb"))
    return top_lines, top_indexes, ew_sun

# Function for finding how many grid points each line will have out of the total allowed
# Output = an array of grid points of the same indexes of top_lines, the output of index_top_transitions
def grid_points(top_indexes, ew_sun):

    # Number of previous wavelength points that were assigned
    total_prev_points = sum(atom2['nq'])
    total_indexed_prev_points = sum(atom2['nq'][top_indexes])
    print(total_prev_points, total_indexed_prev_points)
    print(atom2['nq'][top_indexes])

    # Maximum to allocate
    total_points = 140000

    # Number of points allowed for each transition
    allocated_points = []

    # Index the widths to find the total of all the lines involved. Top 4784 lines (out of 57k) gives 62k/70k
    equiv_widths = np.asarray(ew_sun)[top_indexes]
    og = equiv_widths*1
    print("Using this we have a total Equiv. Width of", sum(equiv_widths), "selected out of", sum(ew_sun))
    #k value to use for smaller line weighting
    k = 0.005
    k_string = str(k).replace(".", "")

    for x in range(len(equiv_widths)):
        # Smaller k upweights smaller values in data
        equiv_widths[x] = 1/(1+np.exp(k*equiv_widths[x]))

    equiv_widths = np.log10(equiv_widths)
    total_width = sum(equiv_widths)

    # Here we modify the values to give a number of frequency points per width/level
    #for width in equiv_widths:
    for w in range(len(equiv_widths)):


        width = equiv_widths[w]
        percentage = (width/total_width)
        level_nq = total_points*percentage
        # round to the nearest odd number
        level_nq = np.ceil(level_nq) // 2 * 2 - 1

        print(og[w], percentage, level_nq)
        allocated_points.append(int(level_nq))
    print(len(np.where(np.asarray(allocated_points) < 28)[0]))
    print(max(allocated_points), min(allocated_points), total_points)


    return allocated_points, k_string, total_points

# Getting the indexes of the top lines (For grid_points), the top lines themselves and the rates/data they have,
# and the list of equiv widths.
top_lines, top_indexes, ew_sun = index_top_transitions()

# Getting the number of grid points each top transition can have. Returns as a list of same length as top lines
grid, k_string, total_points = grid_points(top_indexes, ew_sun)
#rearrange so ascending indexes.
#top_indexes, top_lines, ew_sun = (list(top) for top in zip(*sorted(zip(top_indexes,top_lines, ew_sun))))

atomname = r"D:\PhD\TiFeAtoms\Jack\atom.Ti587_nq2382200.fits"
atom = fits.open(atomname)
atom1 = atom[1].data
atom2 = atom[2].data
indexes = []
# Dictionary to hold nq values [index of atom2] = nqvalue
nq = {}
# Finds indexes in the atom of our top lines. Should be the same indexes I believe?
for x in range(len(top_lines[:, 1])):
    irad = top_lines[x][1]
    jrad = top_lines[x][2]
    index_one = np.where(np.logical_and(irad == atom2['irad'], jrad == (atom2['jrad'])))
    indexes.extend(index_one[0])
    # Making a dictionary with the key of the atom's index to match the nq we want it to have.
    # Must do it this way as we don't loop through the lines themselves but instead loop the new indexes of the atom.
    nq[int(index_one[0])] = grid[x]
indexes = sorted(indexes)
orig_stdout = sys.stdout
# search k = to find the value of k (determines importance of less powerful wavelengths
# The number of freq points used for BB transitions.
print_freq = str(int(total_points/1000))+"k"
abundance = "5.00"
f = open('atom.ti_Nq'+print_freq+'_k'+k_string, 'w')
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
print(" "+abundance+"      47.867 ")
# number of levels, number of bb radiuative transitions, number of bf radiative transitions (Not the number of photons)
#print("    ", len(atom1['ec']), "       ", len(atom2['irad']), "       ", len(radiative), "        ", 0)
print("    ", len(atom1['ec']), "       ", len(indexes), "       ", len(radiative), "        ", 0)

print("***\n*** Levels\n***")
# Levels
for level in range(len(atom1['ec'])):

                                                                        #J is unimportant so I just put 0
    print("   ", f"""{atom1['ec'][level]:<22} {atom1['g'][level]:<10}  {"'" + str(atom1['term'][level]) + ":      " + str(atom1['conf'][level]) + "'":<35s} {atom1['ion'][level]:5} {level+1:5}  {0:5}""")

# Lines, BB radiative
print("***\n*** Lines\n***")

x = 0
#for bb_rad_transition in range(len(atom2['irad'])):             # to convert it to AA from nm we divide by 10
# In indexes, the indexes of the top x number of transitions.
for bb_rad_transition in indexes:  # to convert it to AA from nm we divide by 10

    # Print at 1 level above as balder starts arrays at 1 for some reason. And convert AA to nm
    print(f"{'***':5} {str(bb_rad_transition+1)+':':<5}  {atom2['ALAMB'][bb_rad_transition]/10:11}", "nm")
                                                                                                                                # 0 for atom 1, until we start simplifying stuff. After it is 1 which is fine structure components. Ask karin. IW?#
                                                                                                                                # ga stuff is "broadnening" in anish's, but "dampening" in my email, still fine?
    '''    print("   ", f"{atom2['jrad'][bb_rad_transition]:<10} {atom2['irad'][bb_rad_transition]:<10} {atom2['f'][bb_rad_transition]:<26} "
                 f"{atom2['nq'][bb_rad_transition]:<5} {atom2['qmax'][bb_rad_transition]:<5} {0:<3} {1:<3} "
                 f"{atom2['ga'][bb_rad_transition]:>10} {atom2['gw'][bb_rad_transition]:>5} {atom2['gq'][bb_rad_transition]:>5}")
    '''
    # reducing nq to 1
    print("   ", f"{atom2['jrad'][bb_rad_transition]:<10} {atom2['irad'][bb_rad_transition]:<10} {atom2['f'][bb_rad_transition]:<26} "
                     f"{nq[bb_rad_transition]:<5} {nq[bb_rad_transition]} {0:<3} {1:<3} " # nq and qmax and 0 and 1
                     f"{atom2['ga'][bb_rad_transition]:>10} {atom2['gw'][bb_rad_transition]:>5} {atom2['gq'][bb_rad_transition]:>5}")
    print("***")
    x += 1


# Continua, bf radiative transitions. Remember we must increase the printed levels by 1 to account for Balder starting
# arrays at 1

print("***\n*** Continua\n***")

# One for printing only a single small level transition
"""for bf_rad_transition in range(len(radiative)):
    x = 0
    # size of steps when printing transition. 1 for all photons, 2 for every other, etc. Modifies the length
    # printed of the dictionary specifying the wavelength length. We round up in case it was an odd number that we cut
    # into a decimal with "step"
    steps = 400
    print(" ", f"{int(radiative[bf_rad_transition]['target'])+1:< 7} {bf_rad_transition+1:<7} {radiative[bf_rad_transition]['xsection'][0]:<25} {int(np.ceil(len(radiative[bf_rad_transition]['xsection'])/steps)):<10}", " -1")
    for photon in range(0, len(radiative[bf_rad_transition]['wavelength']), steps):
        x+=1
        wave = radiative[bf_rad_transition]['wavelength'][photon]
        print(" ", f"{radiative[bf_rad_transition]['wavelength'][photon]:<20} {radiative[bf_rad_transition]['xsection'][photon]:<15}")
    break
"""
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