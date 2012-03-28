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
c     ### Initialize GA library
c
c   Find local processor ID and number of processors
c
c     ### Assign local processor ID to variable me
c     ### Assign total number of processors to variable nprocs
c
      write(6,'(a,i4,a,i4,a)') 'Hello world: My rank is ',me,' on ',
     +                       nprocs, ' processors'
c     ### Terminate the GA library
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
