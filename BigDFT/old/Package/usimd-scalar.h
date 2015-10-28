#ifndef __USIMD_SCALAR__
#define __USIMD_SCALAR__

#define _MM_LOAD(a)	*a
#define _MM_STORE(dest, src)	*dest = src

#if defined(PRECISION) && (PRECISION == 1)

#define SIMDFPTYPE float
#define ftype float
#define fptype float
#define MPI_FPTYPE MPI_FLOAT
#define VLEN 1
#define fpconst(n) (n ## f)

#elif defined(PRECISION) && (PRECISION == 2)

#define SIMDFPTYPE double
#define ftype double
#define fptype double
#define MPI_FPTYPE MPI_DOUBLE
#define VLEN 1
#define fpconst(n) (n)

#endif

#endif

