#include <cstdio>
#include <cmath>
#include <mpi.h>

constexpr int n = 840;

int main(int argc, char** argv)
{
  int myid, ntasks;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  printf("Computing approximation to pi with N=%d\n", n);

  int istart = 1;
  int istop = n;
  int midend = n/2;

  double pi = 0.0;
  double  partial_pi = 0.0;
  if (myid == 0){
    for (int i=istart; i <= midend; i++) {
      double x = (i - 0.5) / n;
      pi += 1.0 / (1.0 + x*x);
    }

    MPI_Recv(&partial_pi,1,MPI_DOUBLE,1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);

    pi = partial_pi + pi;
    pi *= 4.0 / n; 
    printf("Approximate pi=%18.16f (exact pi=%10.8f)\n", pi, M_PI);
  } else {
    for (int i=midend + 1; i <= istop; i++) {
      double x = (i - 0.5) / n;
      pi += 1.0 / (1.0 + x*x);
    }
    MPI_Send(&pi,1,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  } 

  
  MPI_Finalize();

}
