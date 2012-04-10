### initialize Message Passing Interface
### initialize Global Arrays

me = ga.nodeid()
nprocs = ga.nnodes()
print "Hello world: My rank is %d on %d processes" % (me, nprocs)

