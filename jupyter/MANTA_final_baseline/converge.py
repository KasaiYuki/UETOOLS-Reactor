import subprocess
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
c = Case('reconverged_HELP_v12_p3.93_n2.38_imp4e2_dp0.06_kep0.06_kip0.0226_last_ii2.hdf5_last_ii2.hdf5')
# bbb.pcoree = 1.45e6
# bbb.pcorei = 1.45e6

# bbb.afracs = 1e-2
print(bbb.afracs, bbb.pcorei)

bbb.isbcwdt=0
proc = subprocess.Popen(['caffeinate', "-d", "-u", "-i", "-s"])

try:
    c.converge(savefname='reconverged_HELP_v12_p3.93_n2.38_imp4e2_dp0.06_kep0.06_kip0.0226', dtreal=1e-10)
    # c.continuation_solve({'pcoree': {'target': 1e6}, 'pcorei': {'target': 1e6}}, savedir='pcore_down-imp_7.5e3', dtreal=1e-10, ii1max=500)
    os.system("printf '\7'") #macos plays sound after converging
finally:
    # Terminate caffeinate
    proc.terminate()