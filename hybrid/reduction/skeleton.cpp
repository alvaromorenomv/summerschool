#include <cstdio>
#include <omp.h>

#define NX 102400

int main(void)
{
    long vecA[NX];
    long sum, sumex;
    int i;

    /* Initialization of the vectors */
    for (i = 0; i < NX; i++) {
        vecA[i] = (long) i + 1;
    }

    sumex = (long) NX * (NX + 1) / ((long) 2);
    printf("Arithmetic sum formula (exact):                  %ld\n",
           sumex);

    sum = 0.0;
    /* TODO: Parallelize computation */
    #pragma omp parallel 
    {
        long psum = 0.0;
        #pragma omp for
        for (i = 0; i < NX; i++) {
            psum += vecA[i];
        }
        #pragma omp critical 
            sum += psum;
        
    }
    printf("Sum: %ld\n", sum);

    return 0;
}
