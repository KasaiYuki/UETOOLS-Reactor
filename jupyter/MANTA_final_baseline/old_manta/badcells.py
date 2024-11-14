import scipy.optimize
import numpy as np
import shapely.geometry
import uedge
from uedge import *
#from uedge import bbb, com

def badCells():
    '''
    Return all (ix, iy) of cells which are not valid polygons, e.g. twisted, all vertices 0, etc.
    '''
    corners = (1,3,4,2,1)
    bad = []
    for ix in range(com.nx+2):
        for iy in range(com.ny+2):
            rcenter = com.rm[ix, iy, 0]
            zcenter = com.zm[ix, iy, 0]
            p = shapely.geometry.Polygon(zip(com.rm[ix,iy,corners]-rcenter, com.zm[ix,iy,corners]-zcenter))
            if not p.is_valid:
                bad.append((ix, iy))
    return bad
