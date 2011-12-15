# If you run "make flags" in your GA build directory, you will get output
# similar to the following. copy-and-paste it here to get the appropriate flag
# values (NOTE: you may need to remove quotation marks from the pasted
# variables.
#
# ** SAMPLE "make flags" output **
# =========================================================================== 
# Suggested compiler/linker options are as follows.
# GA libraries are installed in /Users/d3n000/ga/ga-5-0/bld_openmpi_static/lib
# GA headers are installed in /Users/d3n000/ga/ga-5-0/bld_openmpi_static/include
#
CPPFLAGS=-I/Users/d3n000/ga/ga-5-0/bld_openmpi_static/include
#
LDFLAGS=-L/Users/d3n000/ga/ga-5-0/bld_openmpi_static/lib
#
# For Fortran Programs: 
FFLAGS=-fdefault-integer-8
LIBS=-lga -framework Accelerate
#
# For C Programs: 
CFLAGS=
#LIBS=-lga -framework Accelerate -L/usr/local/lib/gcc/x86_64-apple-darwin10/4.6.0 -L/usr/local/lib/gcc/x86_64-apple-darwin10/4.6.0/../../.. -lgfortran
#
# For C++ Programs: 
CXXFLAGS=
#LIBS=-lga++ -lga -framework Accelerate -L/usr/local/lib/gcc/x86_64-apple-darwin10/4.6.0 -L/usr/local/lib/gcc/x86_64-apple-darwin10/4.6.0/../../.. -lgfortran
# =========================================================================== 

# Compilers
FC = mpif90
CC = mpicc
FFLAGS += -O -g
CFLAGS += -O -g
CPPFLAGS += -DUSE_MPI -DMPI

LINK = $(FC)
LOADER_OPTS = -g

.SUFFIXES: .c .o .h .x

.PHONY: all
all: transp1D.x transp1D-c.x matrix-c.x matrix.x transp1D-c.nb.x

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

.F.o:
	$(FC) $(FFLAGS) $(CPPFLAGS) -c -o $@ $<

.o.x:
	$(LINK) $(LOADER_OPTS) $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) *.o *.x
