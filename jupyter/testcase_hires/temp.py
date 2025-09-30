import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from uedge import *
from uedge.hdf5 import *
#import plotmesh as pm
from uetools import Case
from uedge.gridue import write_gridue
from Forthon import gchange
c = Case('input.yaml')
print(com.nxleg[0,], com.nxcore[0,], com.nycore[0], com.nysol[0])