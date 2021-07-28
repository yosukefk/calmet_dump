# for intel fortran
FFLAGS=-g -traceback -O0 -check bounds
F95=ifort

# for gnu fortran
### FFLAGS=-g 
### F95=gfortran

default: dmp
#default: scn

dmp:
	$(F95) $(FFLAGS) -o calmet_dump dmp_v3.f90
scn:
	$(F95) $(FFLAGS) -o calmet_scan scn.f90
