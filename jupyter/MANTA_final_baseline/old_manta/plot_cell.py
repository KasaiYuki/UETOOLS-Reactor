import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from uedge import bbb, com, api


def plotCellRotated(ix, iy, edge=1):
    '''
    Plot specified cell to determine if it is healthy and not twisted.
    Cell is rotated to make it easier to see. Increment "edge" variable in case
    it's still difficult to see.
    '''
    corners = (1,3,4)
    rcenter = com.rm[ix, iy, 0]
    zcenter = com.zm[ix, iy, 0]
    rc = com.rm[ix,iy,corners]-rcenter
    zc = com.zm[ix,iy,corners]-zcenter
    #theta = -np.arctan2(zc[1+edge]-zc[0+edge], rc[1+edge]-rc[0+edge])
    #c, s = np.cos(theta), np.sin(theta)
    #R = np.array(((c, -s), (s, c)))
    #cRot = R.dot(np.array([rc, zc]))
    #rc = cRot[0,:]
    #zc = cRot[1,:]
    plt.plot(rc, zc)
    plt.scatter(rc[0], zc[0])
    plt.show()
