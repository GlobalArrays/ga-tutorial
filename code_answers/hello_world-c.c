#include <stdio.h>
#include <math.h>
#include "ga.h"
#include "macdecls.h"

#ifdef MPI
#include <mpi.h>
#else
#include "sndrcv.h"
#endif

int main(int argc, char **argv) {
    int me, nprocs;
    
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
    printf("Hello world: My rank is %d on %d processes\n", me, nprocs); fflush(stdout);
       
    GA_Terminate();
    
#ifdef MPI
    MPI_Finalize();    
#else
    PEND_();
#endif
}
