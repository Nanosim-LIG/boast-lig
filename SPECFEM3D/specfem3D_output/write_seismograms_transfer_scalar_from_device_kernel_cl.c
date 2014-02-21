char * write_seismograms_transfer_scalar_from_device_kernel_program = "\
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
__kernel void write_seismograms_transfer_scalar_from_device_kernel(const __global int * number_receiver_global, const __global int * ispec_selected_rec, const __global int * ibool, __global float * station_seismo_field, const __global float * desired_field, const int nrec_local){\n\
  int blockID;\n\
  blockID = get_group_id(0) + (get_group_id(1)) * (get_num_groups(0));\n\
  if(blockID < nrec_local){\n\
    int irec;\n\
    int ispec;\n\
    int iglob;\n\
    irec = number_receiver_global[blockID - 0] - (1);\n\
    ispec = ispec_selected_rec[irec - 0] - (1);\n\
    iglob = ibool[get_local_id(0) + (NGLL3) * (ispec) - 0] - (1);\n\
    station_seismo_field[(NGLL3) * (blockID) + get_local_id(0) - 0] = desired_field[iglob - 0];\n\
  }\n\
}\n\
";