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
      call ga_initialize()
c
c   Find local processor ID and number of processors
c
      me = ga_nodeid()
      nprocs = ga_nnodes()
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
      status = nga_create(MT_INT, ndim, dims, 'test_a', chunk, g_a)
      if (me.eq.0.and.status) then
        write(6,*)
        write(6,'(a)') 'Successfully created Global Array'
        write(6,*)
      endif
c
c   Initialize data in GA. Find data owned by neighboring processor
c
      nghbr = mod(me+1,nprocs)
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
      call nga_put(g_a, lo, hi, a_buf, ld)
      call ga_sync
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
      call nga_distribution(g_a, me, lo, hi)
      call nga_get(g_a, lo, hi, b_buf, ld)
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
      status =  ga_destroy(g_a)
      call ga_terminate
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
