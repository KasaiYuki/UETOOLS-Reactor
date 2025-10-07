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

c = Case('centaur_fb_new.yaml')


# bbb.pcoree=3.83e+6
# bbb.pcorei=3.83e+6
c.populate()
c.about.uedge_setup()
# bbb.isbcwdt=1
# bbb.dtphi=1e-11
# c.interpolate.solution(oldgrid="gridue_v15_orthogonal.hdf5", oldsave="fb_centaur_pcore=7.175_ncore=1.56_last_ii2.hdf5", newgrid="gridue_v15_vertical_coarse_1.hdf5", newsavename="interpolate_upscale_v15_vert.hdf5")
print("HELP2")
# c.converge(savefname="fb_centaur_pcore=7.66", dtreal=1e-10)
c.solver.continuation_solve({'pcoree': {'target': 3.5875e+6}, 'pcorei': {'target': 3.5875e+6}}, savedir='pcore_down', dtreal=1e-10, ii1max=500)

