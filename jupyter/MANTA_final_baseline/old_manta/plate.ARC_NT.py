import uedge
from uedge import *
import numpy as np
from numpy import loadtxt

# for inboard half of mesh:
plate1 = loadtxt("plate1",delimiter=",")

grd.nplate1=len(plate1)
grd.gchange("grd.Mmod",0)

#print len(plate1)
#print plate1[:,1]
print(com.zshift)

grd.rplate1=plate1[:,0]
grd.zplate1=plate1[:,1]+com.zshift

#print grd.zplate1

# for outboard half of mesh:
plate2 = loadtxt("plate2",delimiter=",")

grd.nplate2=len(plate2)
grd.gchange("grd.Mmod",0)

grd.rplate2=plate2[:,0]
grd.zplate2=plate2[:,1]+com.zshift