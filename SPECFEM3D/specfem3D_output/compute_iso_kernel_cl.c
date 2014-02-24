const char * compute_iso_kernel_program = "\
static inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#define INDEX2(xsize,x,y) x + (y)*xsize\n\
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)\n\
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))\n\
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))\n\
#define NDIM 3\n\
#define NGLLX 5\n\
#define NGLL2 25\n\
#define NGLL3 125\n\
#define NGLL3_PADDED 128\n\
#define N_SLS 3\n\
#define IREGION_CRUST_MANTLE 1\n\
#define IREGION_INNER_CORE 3\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#define R_EARTH_KM 6371.0f\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#define BLOCKSIZE_TRANSFER 256\n\
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
