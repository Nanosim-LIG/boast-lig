char * noise_transfer_surface_to_host_kernel_program = "\
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
__kernel void noise_transfer_surface_to_host_kernel(const __global int * ibelm_top, const int nspec_top, const __global int * ibool, const __global float * displ, __global float * noise_surface_movie){\n\
  int igll;\n\
  int iface;\n\
  igll = get_local_id(0);\n\
  iface = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(iface < nspec_top){\n\
    int i;\n\
    int j;\n\
    int k;\n\
    int ispec;\n\
    int iglob;\n\
    ispec = ibelm_top[iface - 0] - (1);\n\
    k = NGLLX - (1);\n\
    j = (igll) / (NGLLX);\n\
    i = igll - ((j) * (NGLLX));\n\
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);\n\
    noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0] = displ[(iglob) * (3) - 0];\n\
    noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0] = displ[(iglob) * (3) - 0];\n\
    noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0] = displ[(iglob) * (3) - 0];\n\
  }\n\
}\n\
";