const char * prepare_boundary_potential_on_device_program = "\
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
__kernel void prepare_boundary_potential_on_device(const __global float * d_potential_dot_dot_acoustic, __global float * d_send_potential_dot_dot_buffer, const int num_interfaces, const int max_nibool_interfaces, const __global int * d_nibool_interfaces, const __global int * d_ibool_interfaces){\n\
  int id;\n\
  int iglob;\n\
  int iloc;\n\
  int iinterface;\n\
  id = get_global_id(0) + (get_global_size(0)) * (get_global_id(1));\n\
  for(iinterface=0; iinterface<=num_interfaces - (1); iinterface+=1){\n\
    if(id < d_nibool_interfaces[iinterface - 0]){\n\
      iloc = id + (max_nibool_interfaces) * (iinterface);\n\
      iglob = d_ibool_interfaces[iloc - 0] - (1);\n\
      d_send_potential_dot_dot_buffer[iloc - 0] = d_potential_dot_dot_acoustic[iglob - 0];\n\
    }\n\
  }\n\
}\n\
";
