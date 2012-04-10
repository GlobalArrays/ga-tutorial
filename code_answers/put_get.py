import mpi4py.MPI # initialize Message Passing Interface
from ga4py import ga # initialize Global Arrays
import numpy as np

NSIZE=100

# 
# Create a Global Array of size NSIZE x NSIZE and fill it with values so that
# the (i,j) element has the value j*NSIZE + i. This corresponds to
# numbering each of the elements consecutively in column-major order.
# 

# Find local processor ID and number of processors
me     = ga.nodeid()
nprocs = ga.nnodes()
if me == 0:
    print "Initialized GA library on %d processes" % nprocs

# Create a GA
dims = (NSIZE,NSIZE)
chunk = (-1,-1)
ld = NSIZE

g_a = ga.create(ga.C_INT, dims, "test_a", chunk)
if me == 0 and g_a:
    print "\nSuccessfully created Global Array"

# Initialize data in GA. Find data owned by neighboring processor

nghbr = (me+1)%nprocs
lo,hi = ga.distribution(g_a, nghbr)

# Create data in local buffer, assign unique value for each data element
patch_shape = hi-lo
a_buf = np.fromfunction(lambda i,j: j*NSIZE + i,
        patch_shape, dtype=ga.dtype(ga.C_INT))
a_buf += lo[1,np.newaxis]
a_buf += lo[np.newaxis,0]*dims[0]

# Copy local data to GA
ga.put(g_a, a_buf, lo, hi)
ga.sync()
if me == 0:
    print "\nCopied values into Global Array from local buffer\n"

# Check data in GA to see if it is correct. Find data owned by this
# processor and then copy it to local buffer
lo,hi = ga.distribution(g_a, me)
b_buf = ga.get(g_a, lo, hi)
if me == 0:
    print "\nCopied values from Global Array to local buffer\n"

# Verify that data is correct
patch_shape = hi-lo
c_buf = np.fromfunction(lambda i,j: j*NSIZE + i,
        patch_shape, dtype=ga.dtype(ga.C_INT))
c_buf += lo[1,np.newaxis]
c_buf += lo[np.newaxis,0]*dims[0]

chk = 1
if not np.all(b_buf == c_buf):
    print "Incorrect value found on process %d" % me
    chk = 0

chk = ga.gop_multiply(chk)
if chk and me == 0:
    print "\nTest of NGA_Put/NGA_Get passed on all processors\n"
elif me == 0:
    print "\nTest of NGA_Put/NGA_Get failed\n"

ga.destroy(g_a)
