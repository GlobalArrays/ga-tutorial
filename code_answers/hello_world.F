      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"

      integer me, nprocs
      integer ierr
c
c   Initialize communication library
c
#ifdef USE_MPI
      call mpi_init(ierr)
#else
      call pbeginf
#endif
c
c   Initialize GA library
c
      call ga_initialize()
c
c   Find local processor ID and number of processors
c
      me = ga_nodeid()
      nprocs = ga_nnodes()
      write(6,'(a,i4,a,i4,a)') 'Hello world: My rank is ',me,' on ',
     +                       nprocs, ' processors'
      call ga_terminate
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
