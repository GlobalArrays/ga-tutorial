# As of ga-5-1, the variables needed to build GA are determined from the
# ga-config script which is installed with GA.

ifndef GA_CONFIG
error:
	echo "you must set GA_CONFIG e.g. GA_CONFIG=/path/to/ga-config"
	exit 1
endif

CC       = $(shell $(GA_CONFIG) --cc)
F77      = $(shell $(GA_CONFIG) --f77)
CFLAGS   = $(shell $(GA_CONFIG) --cflags)
FFLAGS   = $(shell $(GA_CONFIG) --fflags)
CPPFLAGS = $(shell $(GA_CONFIG) --cppflags)
LDFLAGS  = $(shell $(GA_CONFIG) --ldflags)
LIBS     = $(shell $(GA_CONFIG) --libs)
FLIBS    = $(shell $(GA_CONFIG) --flibs)

# =========================================================================== 

FFLAGS += -O -g
CFLAGS += -O -g
CPPFLAGS += -DUSE_MPI -DMPI

LINK = $(F77)
LOADER_OPTS = -g

PROGRAMS =
PROGRAMS += hello_world-c.x
PROGRAMS += hello_world.x
PROGRAMS += matrix-c.x
PROGRAMS += matrix.x
PROGRAMS += put_get-c.x
PROGRAMS += put_get.x
PROGRAMS += transp1D-c.x
PROGRAMS += transp1D.x

.PHONY: all
all: $(PROGRAMS)

.SUFFIXES: .c .o .h .x

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

.F.o:
	$(F77) $(FFLAGS) $(CPPFLAGS) -c -o $@ $<

.o.x:
	$(LINK) $(LOADER_OPTS) $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) *.o *.x
