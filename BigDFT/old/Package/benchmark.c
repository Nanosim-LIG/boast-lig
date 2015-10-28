#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <sys/time.h>
#include <boost/preprocessor/repetition/repeat.hpp>
#define __USE_MISC
#include <sys/mman.h>

#include "convolut_kinetic_per_T_k.h"


#define MULTIPLE_OF_CACHE_LINE(x) ((((x) + 7)/8)*8)

#define OUT_STRIDE1(n1) (n1)
#define OUT_STRIDE2(n2) (n2)
#define OUT_STRIDE3(n3) MULTIPLE_OF_CACHE_LINE(n3)
// Magic offset 64 to avoid conflict misses primarily on KNC
#define OUT_STRIDE23(n2, n3) ((n2) * (OUT_STRIDE3(n3)))

// #define OFFSET_OF_OUTPUT(z, y, x, n2, n3) ((z) * (OUT_STRIDE23((n2), (n3))) + (y) * (OUT_STRIDE3((n3))) + x)
#define OUT_STRIDE21(n2, n1) ((n2) * OUT_STRIDE1(n1))

#define OFFSET_OF_OUTPUT(z, y, x, n2, n1) ((x) * OUT_STRIDE21((n2), (n1)) + (y) * OUT_STRIDE1((n1)) + z)

#define _malloc_2M(X) \
  mmap(NULL, (X + 4096), PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED|MAP_POPULATE, -1, 0)

//   mmap(NULL, (X + 4096), PROT_READ|PROT_WRITE, MAP_ANONYMOUS|MAP_SHARED|MAP_HUGETLB|MAP_POPULATE, -1, 0)
main()
{
    // When we integrate these should be defined by the app already.
    const int lowfil = -14; 
    const int lupfil = 14;

    // Input paramts - this is just for this specific instance of. These parameters
    // will change when it is called from inside BigDFT.
    const int n1 = 64;
    const int n2 = 64;
    const int n3 = 64;
    double *hgrid = (double *)malloc(3 * sizeof(double));
    hgrid[0] = 0.321875000000000;
    hgrid[1] = 0.321875000000000;
    hgrid[2] = 0.321875000000000;
    double k1 = 0.0;
    double k2 = 0.0;
    double k3 = 0.0;
    double *kstrten = (double *)malloc(6 * sizeof(double));


    // This is the output array. This is also read in since it might have some input data
    double *y1 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
    double *y2 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
    // This is the output array. This is also read in since it might have some input data
    double *t1 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
    double *t2 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));

    // This is the input array without skirt - this is short lived
    double *t3 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
    double *t4 = (double *)_malloc_2M(OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));

    // Create input
    FILE *fpx1 = fopen("input_x1", "r");
    FILE *fpx2 = fopen("input_x2", "r");
    FILE *fpy1 = fopen("input_y1", "r");
    FILE *fpy2 = fopen("input_y2", "r");


#if 1
    for (int j = 0; j < n3; ++j)
        for (int i = 0; i < n2; ++i)
            for (int i1 = 0; i1 < n1; ++i1)
            {
                // The temporary inputs don't need a skirt hene use 
                fscanf(fpx1, "%lf\n",
                       &t3[OFFSET_OF_OUTPUT(i1, i, j, n2, n3)]);
                fscanf(fpx2, "%lf\n",
                       &t4[OFFSET_OF_OUTPUT(i1, i, j, n2, n3)]);
                fscanf(fpy1, "%lf\n",
                       &t1[OFFSET_OF_OUTPUT(i1, i, j, n2, n3)]);
                fscanf(fpy2, "%lf\n",
                       &t2[OFFSET_OF_OUTPUT(i1, i, j, n2, n3)]);
            }
#endif

    fclose(fpx1);
    fclose(fpx2);
    fclose(fpy1);
    fclose(fpy2);

    // t3, t4 is x[0] and x[1] in FORTRAN
    // t1, t2 is y[0] and y[1] in FORTRAN
    // This is kind of superfluous but just there if we want to wrap the above with
    // for loop for timing
   
    for (int i = 0; i < 1000; ++i) 
    {
        memcpy (y1, t1, OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
        memcpy (y2, t2, OUT_STRIDE1(n1) * OUT_STRIDE23(n2, n3) * sizeof(double));
    
        convolut_kinetic_per_T_k_wrapper(n1, n2, n3, hgrid, t3, t4, y1, y2, 
                                         kstrten, k1, k2, k3);
    }


    // Verification step
    double error_y1 = 0.0;
    double error_y2 = 0.0;
    {

        FILE *fpy1t = fopen("output_y1", "r");
        FILE *fpy2t = fopen("output_y2", "r");

        double y1t, y2t;
        for (int j = 0; j < (n3); ++j)
            for (int i = 0; i < (n2); ++i)
                for (int i1 = 0; i1 < (n1); ++i1)
                {
                    fscanf(fpy1t, "%lf\n", &y1t);
                    fscanf(fpy2t, "%lf\n", &y2t);
                    error_y2 += fabs(y2[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]-y2t);
                    error_y1 += fabs(y1[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]-y1t);
                    if (fabs(y1[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]-y1t) > fabs(y1t*0.1))
                    {
                        printf("Error at (%d, %d, %d): %24.20f %24.20f\n", i1, i, j, y1t, y1[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]);
                        return -1;
                    }
                    if (fabs(y2[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]-y2t) > fabs(y2t*0.1))
                    {
                        printf("Error at (%d, %d, %d): %24.20f %24.20f\n", i1, i, j, y2t, y2[OFFSET_OF_OUTPUT(i1, i, j, (n2), (n3))]);
                        return -1;
                    }
                }
    }
    printf("Error Checking Pass!!!\nChecksum: (%24.20f, %24.20f)\n", error_y1, error_y2);

    exit(1);
}
