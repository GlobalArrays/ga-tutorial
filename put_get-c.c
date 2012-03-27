#include <stdio.h>
#include <math.h>
#include "ga.h"
#include "macdecls.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

#define NSIZE 100

/**
 * Create a Global Array of size NSIZE x NSIZE and fill it with values so that
 * the (i,j) element has the value j*NSIZE + i. This corresponds to
 * numbering each of the elements consecutively in column-major order.
 */

int main(int argc, char **argv) {
    int me, nprocs, nghbr, ierr;
    int ndim, dims[2], chunk[2], g_a;
    int i, j, ii, jj, lo[2], hi[2], ld, chk;
    int a_buf[NSIZE*NSIZE], b_buf[NSIZE*NSIZE];
    
    /* Initialize Message Passing library */
#ifdef MPI
    MPI_Init(&argc, &argv);   /* initialize MPI */
#else
    PBEGIN_(argc, argv);      /* initialize TCGMSG */
#endif

    /* Initialize GA */
    GA_Initialize();
    
    /* Find local processor ID and number of processors */
    me     = GA_Nodeid();
    nprocs = GA_Nnodes();
    if (me == 0) {
      printf("Initialized GA library on %d processes\n", nprocs); fflush(stdout);
    }

    /* Create a GA */
    ndim = 2;
    dims[0] = NSIZE;
    dims[1] = NSIZE;
    chunk[0] = -1;
    chunk[1] = -1;
    ld = NSIZE;

    g_a = NGA_Create(C_INT, ndim, dims, "test_a", chunk);
    if (me == 0 && g_a) {
      printf("\nSuccessfully created Global Array\n");
    }

    /* Initialize data in GA. Find data owned by neighboring processor */

    nghbr = (me+1)%nprocs;
    NGA_Distribution(g_a, nghbr, lo, hi);

    /* Create data in local buffer */

    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      for (i=lo[0]; i<=hi[0]; i++) {
        ii = i - lo[0];
        /* Assign unique value for each data element */
        a_buf[ii*ld+jj] = i*dims[0] + j;
      }
    }

    /* Copy local data to GA */
    NGA_Put(g_a, lo, hi, a_buf, &ld);
    NGA_Sync();
    if (me == 0) {
      printf("\nCopied values into Global Array from local buffer\n");
    }

    /* Check data in GA to see if it is correct. Find data owned by this
     * processor and then copy it to local buffer */
    NGA_Distribution(g_a, me, lo, hi);
    NGA_Get(g_a, lo, hi, b_buf, &ld);
    if (me == 0) {
      printf("\nCopied values from Global Array to local buffer\n");
    }

    /* Verify that data is correct */
    chk = 1;
    for (j=lo[1]; j<=hi[1]; j++) {
      jj = j - lo[1];
      for (i=lo[0]; i<=hi[0]; i++) {
        ii = i - lo[0];
        /* Assign unique value for each data element */
        if (b_buf[ii*ld+jj] != i*dims[0] + j) {
          printf("Incorrect value found on process %d actual value: %d expected value: %d\n",
              me, b_buf[ii*ld+jj], i*dims[0]+j);
          chk = 0;
        }
      }
    }
    GA_Igop(&chk,1,"*");
    if (chk && me == 0) {
      printf("\nTest of NGA_Put/NGA_Get passed on all processors\n");
    } else if (me == 0) {
      printf("\nTest of NGA_Put/NGA_Get failed\n");
    }

    GA_Destroy(g_a);

    GA_Terminate();
    
#ifdef MPI
    MPI_Finalize();    
#else
    PEND_();
#endif
}
