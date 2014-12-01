#!/usr/bin/python

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

if len(sys.argv)!=2:
    print ('Please, provide a file containing the densities.')
    exit()

fdens = sys.argv[1]

data = np.loadtxt(fdens)
x  = data[:,0]
y  = data[:,1]
z  = data[:,2]
d1 = data[:,3]
d2 = data[:,4]
df = data[:,5]

fig1 = plt.figure( figsize=(8, 6) )
gr1 = fig1.add_subplot(1, 1, 1, projection='3d')
gr1.set_xlabel("x [fm]")
gr1.set_ylabel('y [fm]')
gr1.set_zlabel('density')
gr1.plot_wireframe(x, y, d1, rstride=10, cstride=10, color='red', label='SPH density dist.')
#gr1.plot_wireframe(x, y, d2, rstride=1, cstride=1, color='blue', linestyle='dashed', label='Orig. density dist.')
leg = gr1.legend()

fig2 = plt.figure( figsize=(8, 6) )
gr2 = fig2.add_subplot(1, 1, 1, projection='3d')
gr2.set_xlabel("x [fm]")
gr2.set_ylabel('y [fm]')
gr2.set_zlabel('rel. diff.')
gr2.plot_wireframe(x, y, df, rstride=1, cstride=1, color='black')

plt.show()
