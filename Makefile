
FC = gfortran
FFLAGS = -O3 -fopenmp
LIBSOFLAGS = -shared -fPIC

ifeq ($(strip $(USE_NC)),True)
NETCDF = -DUSE_NC
ifneq ($(strip $(NC_INCLUDE)),)
NETCDFINC = -I $(NC_INCLUDE)
else
NETCDFINC = -I /usr/include
endif
ifneq ($(strip $(NC_LIB)),)
NETCDFLIB = -L$(NC_LIB) -lnetcdff
else
NETCDFLIB = -L/usr/lib/x86_64-linux-gnu -lnetcdff
endif
else
NETCDF = 
NETCDFINC =
NETCDFLIB = 
endif


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

EXECUTABLE = Run_UFiT
LIBSO = UFiT_Python_Callable.so
SLIP_EXECUTABLE = Run_USlip
SLIP_LIBSO = USlip_Python_Callable.so

EXECUTABLE_SRC = UFiT_Definitions_Fortran.F90 UFiT_User_Functions.F90 UFiT_Functions_Fortran.F90 Run_UFiT.F90
LIBSO_SRC = UFiT_Definitions_Fortran.F90 UFiT_User_Functions.F90 UFiT_Functions_Fortran.F90 UFiT_Python_Callable.F90
SLIP_EXECUTABLE_SRC = UFiT_Definitions_Fortran.F90 UFiT_User_Functions.F90 UFiT_Functions_Fortran.F90 USlip_Functions_Fortran.F90 Run_USlip.F90
SLIP_LIBSO_SRC = UFiT_Definitions_Fortran.F90 UFiT_User_Functions.F90 UFiT_Functions_Fortran.F90 USlip_Functions_Fortran.F90 USlip_Python_Callable.F90

# Rule to build the fortran files

all: 
	$(FC) $(EXECUTABLE_SRC) $(FFLAGS) $(NETCDF) $(NETCDFINC) $(NETCDFLIB) -o $(EXECUTABLE)
	$(FC) $(LIBSOFLAGS) $(LIBSO_SRC) $(FFLAGS) $(NETCDF) $(NETCDFINC) $(NETCDFLIB) -o $(LIBSO)
	$(FC) $(SLIP_EXECUTABLE_SRC) $(FFLAGS) $(NETCDF) $(NETCDFINC) $(NETCDFLIB) -o $(SLIP_EXECUTABLE)
	$(FC) $(LIBSOFLAGS) $(SLIP_LIBSO_SRC) $(FFLAGS) $(NETCDF) $(NETCDFINC) $(NETCDFLIB) -o $(SLIP_LIBSO)

# All the dependencies


