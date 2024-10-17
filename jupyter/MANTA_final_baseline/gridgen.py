from uedge import *
from uedge.gridue import write_gridue
from uetools import Case
from Forthon import gchange
from numpy import array


def setGrid(steps=4, plot=True):
    """ Generates a HDF5/ASCII grid

    Keyword arguments
    -----------------
    steps: int (default = 4)
        How far to progress in the grid generation progress
        =0 : Plot efit
        =1 : Plot flux surfaces
        =2 : Plot orthogonal grid
        =3 : Add plates to orthogonal grid
        =4 : Generate and plot grid with user plate shape
    plot : bool (default = True)
        Plot progress for each step if True
    """

    c=Case()
    # Set the geometry
    bbb.mhdgeo = 1 # =1 use MHD equilibrium files - include in input too
    com.geometry = 'dnbot' #  
    com.isudsym = 1 # up-down symmetric - include in input too

    # Set up EFIT equilibrium
    com.geqdskfname[0] = 'FDC24_TokaMaker_2024-10-16_v8_FS_oval.geqdsk'
    com.aeqdskfname[0] = 'aeqdsk_2'

    '''
    com.geqdskfname[0] = 'MANTA_optimized__11122023'
    com.aeqdskfname[0] ='aeqdsk_orig'
    flx.ycutoff1 = 2.20 #1.50 #-2.3 #0.0 #1.50 # set lower cutoff on the eqdsk data
    '''
    if plot:
        c.grid.plot.efit(com.geqdskfname[0].decode('UTF-8'))

    """ FIRST, SET UP FLUX SURFACES """
    if steps>0:
        # Define flux contour boundaries
        flx.psi0min1 = 0.968 # normalized flux on core bndry
        flx.psi0min2 = 0.93 # normalized flux on pf bndry
        flx.psi0max_inner = 1.078 # Outboard side
        flx.psi0max_outer = 1.06 # Inboard side
        flx.psi0max = flx.psi0max_outer # Don't change this
        # !!! psi0sep HAD TO BE MODIFIED ACCORDING TO THE GEQDSK
        flx.psi0sep = 1.01015 # normalized flux at separatrix
        # Normally, this variable is just set to 1+1e-4 to capture
        # a flux surface just outside the separatrix. In this case,
        # it appears the geqdsk file throws off normalization. 
        # Not sure how serious it is, and it's best to wait until
        # the final equilibrium is agreed upon to investigate what
        # value is needed and how a large value impacts things.    

        # Define number of flux tubes
        com.nysol[0] = 27
        com.nyout[0] = 42 # Asymmetric dnull option, needed for plotting
                          # If error is thrown, modify accordingly
        com.nycore[0] = 17

        flx.alfcy = 2.1 # Radial flux-tube distribtuion: higher value, more 
                        # compressed around separatrix
        flx.altsearch = 1 # PFR search option, no touch

     
        flx.flxrun() # Generate flux surface based on the above
        if plot:
            f = c.grid.plot.flx()
            f.get_axes()[0].plot(com.rseps, com.zseps, 'ro')
            f.get_axes()[0].plot(com.rseps2, com.zseps2, 'go')
            f.get_axes()[0].plot(com.rvsin, com.zvsin, 'bo')
            f.get_axes()[0].plot(com.rvsout, com.zvsout, 'bo')


    """ SECOND, SET UP THE POLOIDAL ROWS """
    if steps>1:
        com.nxleg[0] = [16,24] # No cells in inner/outer leg
        com.nxcore[0] = [14,10] # No cells in inner/outer core region
        grd.kxmesh = 1 # Poloidal cell distribution model: =1 uses slpxt
        grd.slpxt = 1.2 # Cell compression factor at plates: higher values, more compressed
        com.ismmon=3 # Controls mesh generation
                      # =1: compressed distribution of mesh on each surface
                      # =2: standard distribution of mesh on each surface
                      # =3: linear combination of 1 & 2 determined by wtmesh1
        grd.wtmesh1=0.75 # =0 -> ismmon=2; =1 -> ismmon=1
        grd.dmix0=1.0 # Mixing length for upstream and downstream grid distros 
        grd.nsmooth=5 # Number of automated smoothing layers applied to mesh 
       

        grd.grdrun() # generate gridue file 
        write_gridue("gridue_orthogonal.hdf5") # Write HDF5 gridue
        c.reload() # Update Case arrays to reflect zshift
        if plot:
            f = c.plot.gridue('gridue_orthogonal.hdf5', plates=False)


    """ DISTORT MESH TO CONFORM TO PLATES """
    if steps>2:
        # !!! Replace below with your coding !!!
        import plates as pl #imports plates.py
        # Store lengths of arrays
        grd.nplate1 = len(pl.rplate1)
        grd.nplate2 = len(pl.rplate2)
        gchange("Mmod") # Allocate space for arrays
        # Set plate specifying-points
        grd.rplate1 = pl.rplate1 
        grd.zplate1 = array(pl.zplate1) + com.zshift
        grd.rplate2 = pl.rplate2
        grd.zplate2 = array(pl.zplate2) + com.zshift
        # !!! Replace above with your coding !!!

        if plot:
            f.get_axes()[0].plot(grd.rplate1, grd.zplate1, 'r-')
            f.get_axes()[0].plot(grd.rplate2, grd.zplate2, 'b-')

    if steps>3:
        grd.iplate=1 # Use user-defined plate geometry in (r/z)plate(1/2)
        # Set up a strike-point search
        grd.isspnew = 1
        grd.isztest = 0
        grd.rstrike = [2.2152, 2.585] #""" !!! Change !!! """
        # Ensure the rstrike is as close to the plate-intersect with the grid in the previous plot

        grd.grdrun() # generate gridue file 
        write_gridue("gridue_shaped.hdf5") # Write HDF5 gridue

        if plot:
            f = c.plot.gridue('gridue_shaped.hdf5', plates=False)
            f.get_axes()[0].plot(grd.rplate1, grd.zplate1, 'r-')
            f.get_axes()[0].plot(grd.rplate2, grd.zplate2, 'b-')

        c.populate() # Populate all UEDGE arrays to allow subsequent evaluations
                     # Note that plots may not work for subsequent executions
