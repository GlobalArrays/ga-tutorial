import mpi4py.MPI # initialize Message Passing Interface
from ga4py import ga # initialize Global Arrays

me = ga.nodeid()
nprocs = ga.nnodes()
print "Hello world: My rank is %d on %d processes" % (me, nprocs)

