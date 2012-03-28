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
    /* ### Initialize the GA library */
    
    /* Find local processor ID and number of processors */
    /* ### assign the processor ID to the variable "me" */
    /* ### assign the total number of processors to the variable "nprocs" */
    printf("Hello world: My rank is %d on %d processes\n", me, nprocs); fflush(stdout);
       
    /* ### terminate the GA library */
    
#ifdef MPI
    MPI_Finalize();    
#else
    PEND_();
#endif
}
