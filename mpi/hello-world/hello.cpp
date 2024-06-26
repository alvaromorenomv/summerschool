#include <iostream>
#include <mpi.h>

int main(int argc, char *argv[])
{
    // TODO: say hello! in parallel
    int size,rank,nodelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(MPI_COMM_WORLD, &nodelen);
 
    std::cout << "Hello! " << "My rank is:" << rank << " My node is:" << name << std::endl;

    if (rank == 0){
        std::cout << "Total number of processes:" << size << std::endl;
    }

    MPI_Finalize();
}
