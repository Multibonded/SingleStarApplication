import random
rush = 0
norush = 0
for z in range(1000000):
    a = []
    c = [1, 2, 3, 4]

    for x in range(7):
        b=random.choice(c)
        if b != 1:
            c.remove(b)
        a.append(b)
    if 2 in a: rush+=1
    else: norush+=1
print(rush, norush)
print(norush/1000000)
