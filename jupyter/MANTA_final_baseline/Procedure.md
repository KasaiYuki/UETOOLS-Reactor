# Procedure
## Use gridgen.py to create hdf5 gridue
* Aeqdsk needs the x point (rb, zb, rt, zt) coordinates at line 20
* Aeqdsk needs strike point (0, rin, zin, rout, zout) coordinates line 58,59
* Gridgen: iterate flx.psi0sep to get the flux surface through x point and strike point
* Make sure the correct plate geometry is imported
## Converge UEDGE
* Add the hdf5 grid file to the YAML input file
* Name the hdf5 save file correctly (name differently than preconverged file)
* Converge with dtreal=1e-10, if doesn't work, set bbb.gamsec = 1 (for striations issue), bbb.isbcwdt = 1 (for converging errors)
* If the initial fnrm starts increasing rapidly when time-step is ~1e-8
* If converged, change input YAML file to point to new save file, if working in the same notebook shouldn't need to