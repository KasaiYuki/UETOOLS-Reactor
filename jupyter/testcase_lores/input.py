from uedge import *
from h5py import File
from uedge.hdf5 import hdf5_restore
from uedge.rundt import *

casename = "input"
savename = "temp_py_puff_last_ii2.hdf5"
aph.aphdir = "/Users/holm10/Documents/fusion/uedge/rates/aph"
api.apidir = "/Users/holm10/Documents/fusion/uedge/rates/api"
bbb.restart=1

# ==== SPECIES ====
com.ngsp=2
com.nhsp=2
bbb.isimpon=7
com.nzsp[0]=6

# ==== GRID ====
bbb.mhdgeo=1
com.geometry="snull"
com.isnonog=1
bbb.gengrid=0
bbb.GridFileName="gridue.hdf5"
com.isgriduehdf5=1

# ==== SOLVER ====
bbb.svrpkg="nksol"
bbb.premeth="ilut"

# ==== PHYSICS ====
bbb.oldseec=0

# ==== BOUNDARYCONDITIONS ====
bbb.isnicore[0]=1
bbb.ncore="2.e19" 
bbb.iflcore=0
bbb.tcoree=250
bbb.tcorei=250
bbb.tedge=2
bbb.t_wall=300
bbb.t_plat=500

# ==== ATOMS ====
bbb.ziin[1]=0
bbb.ineudif=2
bbb.isupgon[0]=1
bbb.isngon[0]=0 
bbb.cfnidh=0.2

# ==== RATES ====
com.istabon=10

# ==== CARBON ====
bbb.ngbackg[1]="0"
bbb.afracs = 1e-3
bbb.ingb=2
bbb.istgcon[1]=1
bbb.tgas[1]=1
bbb.rcxighg=0
bbb.n0g[1]="1.e16"
bbb.isch_sput[1]=7
bbb.isph_sput[1]=3
bbb.csfaclb[2:8]=2.191
bbb.csfacrb[2:8]=2.191
bbb.minu[2:8]=[12, 12, 12, 12, 12, 12]
bbb.ziin[2:8]=[1, 2, 3, 4, 5, 6]
bbb.znuclin[2:8]=[6, 6, 6, 6, 6, 6]
bbb.n0[2:8]=['1.e17', '1.e17', '1.e17', '1.e17', '1.e17', '1.e17']
bbb.fphysylb[1,0]=0.5
bbb.fphysyrb[1,0]=0.5
bbb.fchemylb[1,0]=0.1
bbb.fchemyrb[1,0]=0.1
bbb.nzbackg="1.e9"
bbb.inzb=2
bbb.ismctab=2
api.apidir='/home/sseethalla25/UETOOLS-Reactor/uedge-rates/api'
com.mcfilename[0]="C_rates.strahl"
bbb.isnicore[7]=3
bbb.curcore[7]=0
bbb.isnwcono[2:8]=[3, 3, 3, 3, 3, 3]
bbb.isnwconi[2:8]=[3, 3, 3, 3, 3, 3]
bbb.nwomin[2:8]=['1.e7', '1.e7', '1.e7', '1.e7', '1.e7', '1.e7']
bbb.nwimin[2:8]=['1.e7', '1.e7', '1.e7', '1.e7', '1.e7', '1.e7']
bbb.npltsor=1     #number sources on plates; must be <= 10
bbb.sputtr = 1
bbb.sputtrb=1
# Puffs: indices correspond to source # - limited by npltsor

# Outer wall
bbb.igspsorrb[:,0] =   [2,   0, 0, 0, 0, 0, 0, 0, 0, 0]  # Puffed species index (PYTHON OR BASIS INDICES?)
bbb.igasrb[:,0]=      [10,   0, 0, 0, 0, 0, 0, 0, 0, 0]  # Gas currents from outer wall [Amp]
bbb.ygasrb[:,0]=      [0.02,   0, 0, 0, 0, 0, 0, 0, 0, 0]  # Location of wall source: origin set by issorlb
bbb.wgasrb[:,0]=      [0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0]  # Total Cosine widths of source

# ==== DIFFUSIVITIES ====
bbb.difni[:1]=[1]
bbb.kye=1
bbb.kyi=1
bbb.travis[:1]=[1]

# ==== FLUXLIM ====
bbb.flalfe=0.21
bbb.flalfi=0.21
bbb.flalfv=0.5
bbb.flalfgx=1
bbb.flalfgy=1
bbb.flalfgxy=1
bbb.flalftgx=1
bbb.flalftgy=1
bbb.lgmax=0.1
bbb.lgtmax=0.1
bbb.lgvmax=0.1

# ==== DIFFERENCING ====
bbb.methn=33
bbb.methu=33
bbb.methe=33
bbb.methi=33
bbb.methg=66

# ==== RECYCLING ====
bbb.isoldalbarea=0
bbb.recycp[:2]=[0.95, 0.01]
bbb.recycw[1]=0.0001
bbb.recycm=0.1

bbb.allocate()

# # ==== RESTORE SOLUTION ====
# with File("nc20.hdf5") as f:
#     try:
#         savegroup = f["restore/bbb"]
#     except:
#         savegroup = f["bbb"]
#     for var in ["ngs", "nis", "phis", "tes", "tgs", "tis", "ups"]:
#         setattr(bbb, var, savegroup[var][()])

""" 11: RESTORE SAVE FILE """
bbb.restart=1       # Flag to restart from previous case (=1)
hdf5_restore(savename)     # Restore file
print(bbb.ni.shape)