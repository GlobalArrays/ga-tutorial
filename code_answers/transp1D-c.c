/**
 * transpose of 1-d array.
 * E,g: (1 2 3 4 5 6 7 8 9 10) => (10 9 8 7 6 5 4 3 2 1)
 */
#define   TOTALELEMS   197
#define   MAXPROC      128

#include <stdio.h>
#include <math.h>
#include "ga.h"
#include "macdecls.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

void verify(int g_a, int g_b);

void TRANSPOSE1D() {
    
    int ndim, dims[1], chunk[1], ld[1], lo[1], hi[1];
    int lo1[1], hi1[1], lo2[1], hi2[1];
    int g_a, g_b, a[MAXPROC*TOTALELEMS],b[MAXPROC*TOTALELEMS];
    int nelem, i;    
    int me, nprocs;

    /* Find local processor ID and number of processors */
    me = GA_Nodeid();
    nprocs = GA_Nnodes();
    
    /* Configure array dimensions. Force an unequal data distribution */
    ndim     = 1; /* 1-d transpose */
    dims[0]  = nprocs*TOTALELEMS + nprocs/2;
    ld[0]    = dims[0];
    chunk[0] = TOTALELEMS; /* minimum data on each process */
 
    /* create a global array g_a and duplicate it to get g_b */
    g_a = NGA_Create(C_INT, 1, dims, "array A", chunk);
    if (!g_a) GA_Error("create failed: A", 0);
    if (me==0) printf("  Created Array A\n");
    
    g_b = GA_Duplicate(g_a, "array B");
    if (! g_b) GA_Error("duplicate failed", 0);
    if (me==0) printf("  Created Array B\n");
 
    /* initialize data in g_a */
    if (me==0) {
       printf("  Initializing matrix A\n");
       for(i=0; i<dims[0]; i++) a[i] = i;
       lo[0]  = 0;
       hi[0] = dims[0]-1;
       NGA_Put(g_a, lo, hi, a, ld);
    }

    /* Synchronize all processors to guarantee that everyone has data
       before proceeding to the next step. */
    GA_Sync();

    /* Start initial phase of inversion by inverting the data held locally on
       each processor. Start by finding out which data each processor owns. */
    NGA_Distribution(g_a, me, lo1, hi1);

    /* Get locally held data and copy it into local buffer a  */
    NGA_Get(g_a, lo1, hi1, a, ld);
    
    /* Invert data locally */
    nelem = hi1[0] - lo1[0] + 1;
    for (i=0; i<nelem; i++) b[i] = a[nelem-1-i];
    
    /* Invert data globally by copying locally inverted blocks into
     * their inverted positions in the GA */
    lo2[0] = dims[0] - hi1[0] -1;
    hi2[0] = dims[0] - lo1[0] -1;
    NGA_Put(g_b,lo2,hi2,b,ld);

    /* Synchronize all processors to make sure inversion is complete */
    GA_Sync();

    /* Check to see if inversion is correct */
    if(me == 0) verify(g_a, g_b);
    
    /* Deallocate arrays */
    GA_Destroy(g_a);
    GA_Destroy(g_b);
}

/*
 * Check to see if inversion is correct. Start by copying g_a into local
 * buffer a, and g_b into local buffer b.
 */
void verify(int g_a, int g_b) {

    int i, type, ndim, dims[1], lo[1], hi[1], ld[1];
    int a[MAXPROC*TOTALELEMS],b[MAXPROC*TOTALELEMS];
    
    /* Get dimensions of GA */
    NGA_Inquire(g_a, &type, &ndim, dims);

    lo[0] = 0;
    hi[0] = dims[0]-1;
    NGA_Get(g_a, lo, hi, a, ld);
    NGA_Get(g_b, lo, hi, b, ld);
    
    for(i=0; i<dims[0]; i++)
       if (a[i] != b[dims[0]-i-1]) 
       {
          printf("Mismatch: a[%d]=%d is not equal to b[%d]=%d\n",
                 i, a[i], dims[0]-i-1, b[dims[0]-i-1]);
          GA_Error("verify failed",0);
       }
    
    printf("  Transpose OK\n");
}


int main(int argc, char **argv) {
    int heap=300000, stack=300000;
    int me, nprocs;
    
    /* Step1: Initialize Message Passing library */
#ifdef MPI
    MPI_Init(&argc, &argv);   /* initialize MPI */
#else
    PBEGIN_(argc, argv);      /* initialize TCGMSG */
#endif

    /* Step2: Initialize GA */
    GA_Initialize();
    
    /* Step3: Initialize Memory Allocator (MA) */
    if(! MA_init(C_DBL, stack, heap) ) GA_Error("MA_init failed",stack+heap);

    me     = GA_Nodeid();
    nprocs = GA_Nnodes();
    if(me==0) {
       printf("\nUsing %d processes\n\n", nprocs); fflush(stdout);
    }
    
       
    TRANSPOSE1D();
    
    if(me==0)printf("\nTerminating ..\n");
    GA_Terminate();
    
#ifdef MPI
    MPI_Finalize();    
#else
    PEND_();
#endif
}
