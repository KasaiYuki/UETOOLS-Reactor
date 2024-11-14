# Plotting UEDGE computational mesh
#
#-Expects imported modules:
#import matplotlib.pyplot as plt; import numpy as np
#
#
#-Usage:
# execfile("../pylib/plotmesh.py") 
# plotmesh()
#
#-Optional arguments:
#   iso (True/False) - True for equal aspect ratio
#
#
# First coding: MVU, 17-jul-17
#=========================================================#

import numpy as np
import matplotlib.pyplot as plt
import uedge
from uedge import *


def plotmesh(iso=True):

    #fig,ax = plt.subplots(1)

    if (iso):
        plt.axes().set_aspect('equal', 'datalim')
    else:
        plt.axes().set_aspect('auto', 'datalim')


    for iy in np.arange(0,com.ny+2):
        for ix in np.arange(0,com.nx+2):
            plt.plot(com.rm[ix,iy,[1,2,4,3,1]],
                     com.zm[ix,iy,[1,2,4,3,1]], 
                     color="k",
		     lw=0.4)


    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    #fig.suptitle('UEDGE mesh')
    plt.grid(True)

    #read fwall limits
    #fwall = np.loadtxt("First_Wall_Coordinates.txt",skiprows=2)
    #rwall = np.zeros((len(fwall),2))
    #zwall = np.zeros((len(fwall),2))
    #for i in range(0,len(fwall)):
    #    rwall[i]=fwall[i][0]
    #    zwall[i]=fwall[i][1]+com.zshift
    #plt.plot(rwall,zwall,'r')


    plt.show()

#=========================================================#
