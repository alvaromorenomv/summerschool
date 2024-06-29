#include <cstdio>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int my_id, omp_rank;
    int provided, required=MPI_THREAD_FUNNELED;

    /* TODO: Initialize MPI with thread support. */
    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);

    /* TODO: Find out the MPI rank and thread ID of each thread and print
     *       out the results. */
    #pragma omp parallel
    {
        omp_rank = omp_get_thread_num();
        printf("Hello World! my MPI Rank:%d my Thread ID:%d\n",my_id,omp_rank);
    }

    /* TODO: Investigate the provided thread support level. */
    printf("Provided thread support level is %d\n",provided);

    MPI_Finalize();
    return 0;
}
