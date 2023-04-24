.POSIX:
.SUFFIXES: .o .cpp .f90
.PHONY: clean all

PREFIX = $(SMARTIES_ROOT)
CXXFLAGS = -O0 -g
FCFLAGS = -O0 -g
 
FC = mpifort
CXX = mpicxx
LINK = $(FC)

M_CXXFLAGS = -I$(PREFIX)/include
M_FCFLAGS = -I$(PREFIX)/include -I.
MPI_EXTRA_LIB = -lmpi_cxx -lstdc++
M_LDFLAGS = -fopenmp -L$(PREFIX)/lib -lsmarties $(MPI_EXTRA_LIB)
#M_CHFLAGS = -lfftw3_mpi -lfftw3 -L$(OPENBLAS_ROOT)/lib -lopenblas -I$(OPENBLAS_ROOT) -I$(FFTW_INCLUDE) -L$(FFTW_LIB)
M_CHFLAGS = -lfftw3_mpi -lfftw3 -lopenblas -I$(FFTW_INCLUDE) -L$(FFTW_LIB)

C = \
mpi.o\
global.o\
mass_flow.o\
interpolation.o\
equations.o\
boundary_conditions.o\
pressure.o\
input_output.o\
statistics.o\
initialization.o\
subgrid.o\
wallmodel.o\
finalization.o\
projection.o\
time_integration.o\
monitor.o\
smarties_stat.o\

O = \
app_main.o\
main.o\
smarties_extern.o\

M = \
app_main\

all: $C $M
$M: $O
.cpp.o:; $(CXX) -c $(CXXFLAGS) $(M_CXXFLAGS) $<
.f90.o:; $(FC) -c $(FCFLAGS) $(M_FCFLAGS) $(M_CHFLAGS) $<
.o:; $(LINK) $C $O $(LDFLAGS) $(M_LDFLAGS) $(M_CHFLAGS) -o $@
clean:; rm -f -- $M $O $C *.mod
