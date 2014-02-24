const char * noise_add_source_master_rec_kernel_program = "\
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
__kernel void noise_add_source_master_rec_kernel(const __global int * ibool, const __global int * ispec_selected_rec, const int irec_master_noise, __global float * accel, const __global float * noise_sourcearray, const int it){\n\
  int tx;\n\
  int ispec;\n\
  int iglob;\n\
  tx = get_local_id(0);\n\
  ispec = ispec_selected_rec[irec_master_noise - 0] - (1);\n\
  iglob = ibool[tx + (NGLL3) * (ispec) - 0] - (1);\n\
  atomicAdd(accel + (iglob) * (3) + 0, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 0 - 0]);\n\
  atomicAdd(accel + (iglob) * (3) + 1, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 1 - 0]);\n\
  atomicAdd(accel + (iglob) * (3) + 2, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 2 - 0]);\n\
}\n\
";
