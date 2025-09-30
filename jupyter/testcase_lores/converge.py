from input import *
from uedge.rundt import *
from uedge.hdf5 import *
import matplotlib.pyplot as plt
from uetools import *
import sys
hdf5_restore('temp_py_puff_last_ii2.hdf5')
bbb.exmain()
Case().populate()
print(bbb.ni[0,0,:])
print(com.nzsp)
print(bbb.isngon)
print(bbb.isnion)
print(bbb.isnicore)