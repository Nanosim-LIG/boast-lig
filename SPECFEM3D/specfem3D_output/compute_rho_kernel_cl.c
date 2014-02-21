char * compute_rho_kernel_program = "\
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
__kernel void compute_rho_kernel(const __global int * ibool, const __global float * accel, const __global float * b_displ, __global float * rho_kl, const int NSPEC, const float deltat){\n\
  int ispec;\n\
  int ijk_ispec;\n\
  int iglob;\n\
  ispec = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(ispec < NSPEC){\n\
    ijk_ispec = get_local_id(0) + (NGLL3) * (ispec);\n\
    iglob = ibool[ijk_ispec - 0] - (1);\n\
    rho_kl[ijk_ispec - 0] = rho_kl[ijk_ispec - 0] + (deltat) * ((accel[0 - 0 + (iglob - (0)) * (3)]) * (b_displ[0 - 0 + (iglob - (0)) * (3)]) + (accel[1 - 0 + (iglob - (0)) * (3)]) * (b_displ[1 - 0 + (iglob - (0)) * (3)]) + (accel[2 - 0 + (iglob - (0)) * (3)]) * (b_displ[2 - 0 + (iglob - (0)) * (3)]));\n\
  }\n\
}\n\
";