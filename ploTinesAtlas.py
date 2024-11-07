
wl = []
intensity =  []

with open("ftsatlas.txt", "r") as f:
    for line in f:
        wl.append(float(line.split()[0]))
        intensity.append(float(line.split()[1]))
import matplotlib.pyplot as plt

plt.plot(wl, intensity)
plt.show()