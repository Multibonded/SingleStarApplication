

f = open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\122563\wl_122563")


d = open(r"D:\PhD\TiFeAtoms\Jack\balder\Balder_outputs\122563\halved_wl_122563", "w")

y = 0
for x in f:
    if y%2==0:
        d.write(x)

    y+=1