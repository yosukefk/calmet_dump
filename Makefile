# for intel fortran
NCDDIR=${TACC_NETCDF_DIR}
FFLAGS=-g -traceback -O0 -check bounds -I$(NCDDIR)/include
F95=ifort
LIBS=-L$(NCDDIR)/lib -lnetcdff -lnetcdf


# for gnu fortran
### FFLAGS=-g 
### F95=gfortran

default: dmp
#default: scn

dmp:
	$(F95) $(FFLAGS) -o calmet_dump dmp.f90 coordlib.for $(LIBS)
scn:
	$(F95) $(FFLAGS) -o calmet_scan scn.f90

test_coords:  test_coords.f90
	$(F95) $(FFLAGS) -o test_coords test_coords.f90 coordlib.for
