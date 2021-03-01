import numpy as np

LS_R = 18000 #mm
print("LS diameter = %d" %(LS_R*2) )

length = 5000 
side = []

xx, zz = [], []
counter = 0

start = 0 - 3*length

for i in range(7):
    for j in range(7):
        for k in range(7):
            x = start + i*length
            y = start + j*length
            z = start + k*length
            if np.sqrt(x**2 + y**2 + z**2) < 17700:
                counter += 1
                print("%d %d %d"%(x, y, z) )
                if y == 0 :
                    xx.append(x)
                    zz.append(z)

gridx = [-3.5*length+length*i for i in range(8)]
gridz = [-3.5*length+length*j for j in range(8)]


print("total cubic volume number = %d" %counter )

import matplotlib.pyplot as plt

fig, ax = plt.subplots()

circle1 = plt.Circle((0, 0), 17700, color='blue', fill=False)
ax.add_patch(circle1)
ax.set_xlim(-20000, 20000)
ax.set_ylim(-20000, 20000)

ax.plot(xx, zz, "o", ms=2, color='r')

for i in gridx:
    plt.vlines(i, gridz[0], gridz[-1], color='r')
for i in gridz:
    plt.hlines(i, gridx[0], gridx[-1], color='r')

plt.show()
