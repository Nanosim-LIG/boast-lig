const char * update_potential_kernel_program = "\
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
__kernel void update_potential_kernel(__global float * potential_acoustic, __global float * potential_dot_acoustic, __global float * potential_dot_dot_acoustic, const int size, const float deltat, const float deltatsqover2, const float deltatover2){\n\
  int id;\n\
  id = get_global_id(0) + (get_group_id(1)) * (get_global_size(0));\n\
  if(id < size){\n\
    potential_acoustic[id - 0] = potential_acoustic[id - 0] + (deltat) * (potential_dot_acoustic[id - 0]) + (deltatsqover2) * (potential_dot_dot_acoustic[id - 0]);\n\
    potential_dot_acoustic[id - 0] = potential_dot_acoustic[id - 0] + (deltatover2) * (potential_dot_dot_acoustic[id - 0]);\n\
    potential_dot_dot_acoustic[id - 0] = 0.0f;\n\
  }\n\
}\n\
";
