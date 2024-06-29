#include <cstdio>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int my_id, omp_rank, ntasks, msg, i;
    int provided, required=MPI_THREAD_FUNNELED;

    MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
    MPI_Comm_rank(MPI_COMM_WORLD,&my_id);
    MPI_Comm_size(MPI_COMM_WORLD,&ntasks);


    #pragma omp parallel private(msg,i,omp_rank)
    {   
        omp_rank = omp_get_thread_num();
        if (my_id == 0){
            for (i = 1;i<ntasks;i++)
            MPI_Send(&omp_rank,1,MPI_INT,i,omp_rank,MPI_COMM_WORLD);
        } else{
            MPI_Recv(&msg,1,MPI_INT,0,omp_rank,MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("my MPI Rank:%d my Thread ID:%d Task 0 colleague ID:%d\n",my_id,omp_rank,msg);
        }

    }


    MPI_Finalize();
    return 0;
}
