#define   TOTALELEMS 1007
#define   MAXPROC    128
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"

      integer ndim, dims(2), chunk(2), nprocs, me, i, j, k, l
      integer lo(2), hi(2), lo1(2), hi1(2), lo2(2), hi2(2)
      integer lo3(2), hi3(2), ld(2), nelem
      integer g_a, g_b, g_c, g_chk
      double precision a(TOTALELEMS,TOTALELEMS)
      double precision b(TOTALELEMS,TOTALELEMS)
      double precision c(TOTALELEMS,TOTALELEMS)
      double precision atrns(TOTALELEMS,TOTALELEMS)
      double precision rchk
      integer heap, stack, ichk, ierr
      logical status

      heap = 3000000
      stack = 3000000
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
c     ### call initialization routine for GA library
c
c   Find local processor ID and number of processors
c
c     ### assign the local node ID to the integer variable "me"
c
c     ### assign the total number of processors to the
c     ### integer variable "nprocs"
c
c   Allocate memory for GA library
c
      status = ma_init(MT_F_DBL, stack/nprocs, heap/nprocs)
c
c   Configure array dimensions. Force an unequal data distribution.
c
      ndim = 2
      dims(1) = TOTALELEMS
      dims(2) = TOTALELEMS
      ld(1) = TOTALELEMS
      ld(2) = TOTALELEMS
      chunk(1) = TOTALELEMS/nprocs - 1     ! Minimum block dimension on
      chunk(2) = TOTALELEMS/nprocs - 1     ! each processor
c
c   Create global array g_a and then duplicate it to get g_b, g_c, and g_chk
c
c     ### create GA of doubles with dimension "NDIM" and size "dims" with
c     ### minimum block size "chunk" and assign the handle to the
c     ### integer variable "g_a". Then create remaining global arrays,
c     ### assigned to the integer handle "g_b", "g_c", and "g_chk",
c     ### by duplicating g_a. Assign the names "Array A", "Array B", "Array C",
c     ### and "Array Check" to "g_a", "g_b", "g_c", and "g_chk".
c
c   Initialize data in matrices a and b
c
      k = 0
      l = 7
      do j = 1, dims(2)
        do i =  1, dims(1)
          k = k + 1
          k = mod(k,29)
          a(i,j) = dble(k)
          l = l + 1
          l = mod(l,37)
          b(i,j) = dble(l)
        end do
      end do
      lo1(1) = 1
      lo1(2) = 1
      hi1(1) = dims(1)
      hi1(2) = dims(2)
c
c   Copy data to global arrays g_a and g_b
c
      if (me.eq.0) then
c     ### copy the contents of array "a" into the portion of global array
c     ### "g_a" described by "lo1" and "hi1". Similarly, copy the contents
c     ### of the array "b" into corresponding portion of global array "g_b".
c     ### Use the array of strides "ld" to describe the physical layout of
c     ### arrays "a" and "b".
      endif
c
c   Synchronize all processors to make sure everyone has data
c
c     ### synchronize all processors 
c
c   Determine which block of data is locally owned. Note that
c   the same block is locally owned for all GAs.
c
c     ### find out which block of data my node ("me") owns for the global
c     ### array "g_c" and store the contents in the integer arrays "lo" and
c     ### "hi".
c
c   Get the blocks from g_a and g_b needed to compute this block in
c   g_c and copy them into the local buffers a and b.
c
      lo2(1) = lo(1)
      lo2(2) = 1
      hi2(1) = hi(1)
      hi2(2) = dims(2)
c     ### copy the block of data described by the arrays "lo2" and "hi2" from
c     ### the global array "g_a" into the local array "a". Use the array of
c     ### strides "ld" to describe the physical layout of "a".
      lo3(1) = 1
      lo3(2) = lo(2)
      hi3(1) = dims(1)
      hi3(2) = hi(2)
c     ### copy the block of data described by the arrays "lo3" and "hi3" from
c     ### the global array "g_b" into the local array "b". Use the array of
c     ### strides "ld" to describe the physical layout of "b".
c
c   Do local matrix multiplication and store the result in local
c   buffer c. Start by evaluating the transpose of a.
c
      do j = 1, hi2(2) - lo2(2) + 1
        do i = 1, hi2(1) - lo2(1) + 1
          atrns(j,i) = a(i,j)
        end do
      end do
c
c   Multiply a and b to get c
c
      do i = 1, hi(1) - lo(1) + 1
        do j = 1, hi(2) - lo(2) + 1
          c(i,j) = 0.0d00
          do k = 1, dims(1)
            c(i,j) = c(i,j) + atrns(k,i)*b(k,j)
          end do
        end do
      end do
c
c   Copy c back to g_c
c
c     ### copy data from the local array "c" into the block of the global
c     ### array "g_c" described by the integer arrays "lo" and "hi". Use
c     ### the array of strides "ld" to describe the physical layout of "c".
c
c   Verify that g_c equals the correct answer. Use GA operations
c   that support matrix multiplies, matrix additions, and dot
c   products to perform check.
c
      call ga_dgemm('n','n',TOTALELEMS,TOTALELEMS,TOTALELEMS,1.0d00,
     +              g_a,g_b,0.0d00,g_chk)
c
c   Take difference between matrices
c
      call ga_add(1.0d00,g_c,-1.0d00,g_chk,g_chk)
      rchk =  ga_ddot(g_chk,g_chk)
      if (me.eq.0) write(6,100) rchk
  100 format('Normed difference in matrices: ',e12.4)
      if (rchk.lt.0.1.and.me.eq.0) then
        write(6,*) 'Matrix multiply is successful'
      elseif (me.eq.0) then
        write(6,*) 'Matrix multiply is unsuccessful'
      endif
c
c   Deallocate memory for arrays and clean up GA library
c
c     ### destroy the global arrays "g_a", "g_b", "g_c", "g_chk"
c     ### terminate the GA library
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
