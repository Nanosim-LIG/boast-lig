const char * compute_iso_kernel_program = "\
inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#ifndef INDEX2\n\
INDEX2(isize,i,j) i + isize*j\n\
#endif\n\
#ifndef INDEX3\n\
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)\n\
#endif\n\
#ifndef INDEX4\n\
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))\n\
#endif\n\
#ifndef INDEX5\n\
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))\n\
#endif\n\
#ifndef NDIM\n\
#define NDIM 3\n\
#endif\n\
#ifndef NGLLX\n\
#define NGLLX 5\n\
#endif\n\
#ifndef NGLL2\n\
#define NGLL2 25\n\
#endif\n\
#ifndef NGLL3\n\
#define NGLL3 125\n\
#endif\n\
#ifndef NGLL3_PADDED\n\
#define NGLL3_PADDED 128\n\
#endif\n\
#ifndef N_SLS\n\
#define N_SLS 3\n\
#endif\n\
#ifndef IREGION_CRUST_MANTLE\n\
#define IREGION_CRUST_MANTLE 1\n\
#endif\n\
#ifndef IREGION_INNER_CORE\n\
#define IREGION_INNER_CORE 3\n\
#endif\n\
#ifndef IFLAG_IN_FICTITIOUS_CUBE\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#endif\n\
#ifndef R_EARTH_KM\n\
#define R_EARTH_KM 6371.0f\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_INNER_CORE\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#endif\n\
#ifndef COLORING_MIN_NSPEC_OUTER_CORE\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#endif\n\
#ifndef BLOCKSIZE_TRANSFER\n\
#define BLOCKSIZE_TRANSFER 256\n\
#endif\n\
__kernel void compute_iso_kernel(const __global float * epsilondev_xx, const __global float * epsilondev_yy, const __global float * epsilondev_xy, const __global float * epsilondev_xz, const __global float * epsilondev_yz, const __global float * epsilon_trace_over_3, const __global float * b_epsilondev_xx, const __global float * b_epsilondev_yy, const __global float * b_epsilondev_xy, const __global float * b_epsilondev_xz, const __global float * b_epsilondev_yz, const __global float * b_epsilon_trace_over_3, __global float * mu_kl, __global float * kappa_kl, const int NSPEC, const float deltat){\n\
  int ispec;\n\
  int ijk_ispec;\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(ispec < NSPEC){\n\
    ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
    mu_kl[ijk_ispec - 0] = mu_kl[ijk_ispec - 0] + (deltat) * ((epsilondev_xx[ijk_ispec - 0]) * (b_epsilondev_xx[ijk_ispec - 0]) + (epsilondev_yy[ijk_ispec - 0]) * (b_epsilondev_yy[ijk_ispec - 0]) + (epsilondev_xx[ijk_ispec - 0] + epsilondev_yy[ijk_ispec - 0]) * (b_epsilondev_xx[ijk_ispec - 0] + b_epsilondev_yy[ijk_ispec - 0]) + ((epsilondev_xy[ijk_ispec - 0]) * (b_epsilondev_xy[ijk_ispec - 0]) + (epsilondev_xz[ijk_ispec - 0]) * (b_epsilondev_xz[ijk_ispec - 0]) + (epsilondev_yz[ijk_ispec - 0]) * (b_epsilondev_yz[ijk_ispec - 0])) * (2));\n\
    kappa_kl[ijk_ispec - 0] = kappa_kl[ijk_ispec - 0] + (deltat) * (((epsilon_trace_over_3[ijk_ispec - 0]) * (b_epsilon_trace_over_3[ijk_ispec - 0])) * (9));\n\
  }\n\
}\n\
";
