#include <cstdio>
#include <omp.h>
int main()
{
    printf("Hello world!\n");
#pragma omp parallel
    {
        int omp_nthreads = omp_get_num_threads();
        int omp_thread_num = omp_get_thread_num();
        printf("Number of threads:%d Thread number: %d\n",
            omp_nthreads,
            omp_thread_num);
    }
    return 0;
}
