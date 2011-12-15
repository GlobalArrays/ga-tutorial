#define   NDIM         1
#define   TOTALELEMS 197
#define   MAXPROC    128
      program main
      implicit none
#include "mafdecls.fh"
#include "global.fh"

      integer dims(3), chunk(3), nprocs, me, i, lo(3), hi(3), lo1(3)
      integer hi1(3), lo2(3), hi2(3), ld(3), nelem
      integer g_a, g_b, a(MAXPROC*TOTALELEMS), b(MAXPROC*TOTALELEMS)
      integer heap, stack, ichk, ierr
      logical status

      heap = 300000
      stack = 300000
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
      if (me.eq.0) write(6,101) nprocs
  101 format('Using ',i4,' processors')
c
c   Allocate memory for GA library
c
      status = ma_init(MT_F_DBL, stack/nprocs, heap/nprocs)
c
c   Configure array dimensions. Force an unequal data distribution.
c
      dims(1) = nprocs*TOTALELEMS + nprocs/2
      ld(1) = MAXPROC*TOTALELEMS
      chunk(1) = TOTALELEMS     ! Minimum data on each processor
c
c   Create global array g_a and then duplicate it to get g_b
c
c     ### create GA of integers with dimension "NDIM" and size "dims" with
c     ### minimum block size "chunk" and assign the handle to the
c     ### integer variable "g_a". Then create a second global array
c     ### assigned to the integer handle "g_b" by duplicating g_a. Assign
c     ### the names "Array A" and "Array B" to "g_a" and "g_b".
c
c   Initialize data in g_a
c
      do i = 1, dims(1)
        a(i) = i
      end do
      lo1(1) = 1
      hi1(1) = dims(1)
c
c   Copy data from local buffer a to global array g_a. Only do this for
c   processor 0.
c
      if (me.eq.0) then
c     ### copy the contents of array "a" into the portion of global array
c     ### "g_a" described by "lo1" and "hi1". Use the array of strides
c     ### "ld" to describe the physical layout of array "a".
      endif
c
c   Synchronize all processors to guarantee that everyone has data
c   before proceeding to the next step.
c
c     ### synchronize all processors
c
c   Start initial phase of inversion by inverting the data held locally on
c   each processor. Start by finding out which data each processor owns.
c
c     ### find out which block of data my node ("me") owns for the global
c     ### array "g_a" and store the contents in the integer arrays "lo" and
c     ### "hi".
c
c   Get locally held data and copy it into local buffer a
c
c     ### use the arrays "lo" and "hi" to copy the locally held block of data
c     ### from the global array "g_a" into the local array "a". Use the array
c     ### of strides "ld" to describe the physical layout of "a".
c
c   Invert local data
c
      nelem = hi(1) - lo(1) + 1
      do i = 1, nelem
        b(i) = a(nelem - i + 1)
      end do
c
c   Do global inversion by copying locally inverted data blocks into
c   their inverted positions in the GA
c
      lo2(1) = dims(1) - hi(1) + 1
      hi2(1) = dims(1) - lo(1) + 1
c     ### copy data from the local array "b" into the block of the global
c     ### array "g_a" described by the integer arrays "lo2" and "hi2". Use
c     ### the array of strides "ld" to describe the physical layout of "b".
c
c   Synchronize all processors to make sure inversion is complete
c
c     ### synchronize all processors
c
c   Check to see if inversion is correct. Start by copying g_a into local
c   buffer a, and g_b into local buffer b.
c
c     ### copy the block of data described by the arrays "lo1" and "hi1" from
c     ### the global array "g_a" into the local array "a". Copy the same block
c     ### of data from "g_b" into the local array "b". Use the array of strides
c     ### "ld" to describe the physical layout of "a" and "b".
      ichk = 0
      do i = 1, dims(1)
        if (a(i).ne.b(dims(1)-i+1) .and. me.eq.0) then
          write(6,111) i,a(i),b(dims(1)-i+1)
  111     format('Mismatch at ',3i8)
          ichk = ichk + 1
        endif
      end do
      if (ichk.eq.0.and.me.eq.0) write(6,*) 'Transpose OK'
c
c   Deallocate memory for arrays and clean up GA library
c
      if (me.eq.0) write(6,*) 'Terminating...'
c     ### destroy global arrays "g_a" and "g_b"
c     ### terminate the GA library
      call ga_terminate
#ifdef USE_MPI
      call mpi_finalize
#else
      call pend
#endif
      stop
      end
