#define NSIZE 100
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"

      integer me, nprocs, nghbr, ierr
      integer ndim, dims(2), chunk(2), g_a
      integer i, j, ii, jj, lo(2), hi(2), ld, chk
      integer a_buf(NSIZE,NSIZE), b_buf(NSIZE,NSIZE)
      logical status
c
c   Create a Global Array of size NSIZE x NSIZE and fill it with values so that
c   the (i,j) element has the value (j-1)*NSIZE + i. This corresponds to
c   numbering each of the elements consecutively in column-major order.
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
c     ### initialize the GA library
c
c   Find local processor ID and number of processors
c
c     ### assign the processor ID to the local variable "me"
c     ### assign the total number of processors to the local variable "nprocs"
      if (me.eq.0) then
        write(6,'(a,i4,a)') 'Intialized GA library on ',nprocs,
     +                      ' processes'
      endif
c
c   Create a GA
c
      ndim = 2
      dims(1) = NSIZE
      dims(2) = NSIZE
      chunk(1) = -1
      chunk(2) = -1
      ld = NSIZE
c     ### create a global array of type MT_INT and dimension "ndim", with size
c     ### "dims" and minimum block size "chunk" and assign it to the array handle
c     ### "g_a"
      if (me.eq.0.and.status) then
        write(6,*)
        write(6,'(a)') 'Successfully created Global Array'
        write(6,*)
      endif
c
c   Initialize data in GA. Find data owned by neighboring processor
c
      nghbr = mod(me+1,nprocs)
c     ### find the lowest and highest indices of the local data on processor "nghbr"
c     ### and store them in the arrays "lo" and "hi"
      call nga_distribution(g_a, nghbr, lo, hi)
c
c   Create data in local buffer
c
      do j = lo(2), hi(2)
        jj = j - lo(2) + 1
        do i = lo(1), hi(1)
          ii = i - lo(1) + 1
c
c    Assign unique value for each data element
c
          a_buf(ii,jj) = (j-1)*dims(1) + i
        end do
      end do
c
c   Copy local data to GA
c
c     ### put the contents of local buffer "a_buf(ld,*)" into the global array
c     ### "g_a" in the block denoted by the arrays "lo" and "hi"
c
c     ### synchronize the entire system
      if (me.eq.0) then
        write(6,*)
        write(6,'(a)')
     +        'Copied values into Global Array from local buffer'
        write(6,*)
      endif
c
c   Check data in GA to see if it is correct. Find data owned by this
c   processor and then copy it to local buffer
c
c     ### find the lowest and highest indices on processor "me" and store
c     ### them in the arrays "lo" and "hi"
c
c     ### retrieve the data from "g_a" corresponding to the block represented
c     ### by the arrays "lo" and "hi" and store it in the local buffer "b_buf(ld,*)
      if (me.eq.0) then
        write(6,*)
        write(6,'(a)') 'Copied values from Global Array to local buffer'
        write(6,*)
      endif
c
c   Verify that data is correct
c
      chk = 1
      do j = lo(2), hi(2)
        jj = j - lo(2) + 1
        do i = lo(1), hi(1)
          ii = i - lo(1) + 1
          if (b_buf(ii,jj).ne.(j-1)*dims(1)+i) then
            write(6,'(a,i4,a,i8,a,i8)')
     +          'Incorrect value found on process ',me,
     +          ' actual value: ',b_buf(ii,jj),
     +          ' expected value: ',(j-1)*dims(1)+i
            chk = 0
          endif
        end do
      end do
c
      call ga_igop(1,chk,1,'*')
      if (chk.eq.1.and.me.eq.0) then
        write(6,*)
        write(6,'(a)') 'Test of nga_put/nga_get passed on all processes'
        write(6,*)
      elseif (me.eq.0) then
        write(6,*)
        write(6,'(a)') 'Test of nga_put/nga_get failed'
        write(6,*)
      endif
c
c   Destroy array
c
c     ### destroy the GA corresponding to the handle "g_a"
c
c     ### terminate the GA library
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
