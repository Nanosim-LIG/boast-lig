#define INDEX2(xsize,x,y) x + (y)*xsize
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#define NDIM 3
#define NGLLX 5
#define NGLL2 25
#define NGLL3 125
#define NGLL3_PADDED 128
#define N_SLS 3
#define IREGION_CRUST_MANTLE 1
#define IREGION_INNER_CORE 3
#define IFLAG_IN_FICTITIOUS_CUBE 11
#define R_EARTH_KM 6371.0f
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#define BLOCKSIZE_TRANSFER 256
__global__ void noise_add_source_master_rec_kernel(const int * ibool, const int * ispec_selected_rec, const int irec_master_noise, float * accel, const float * noise_sourcearray, const int it){
  int tx;
  int ispec;
  int iglob;
  tx = threadIdx.x;
  ispec = ispec_selected_rec[irec_master_noise - 0] - (1);
  iglob = ibool[tx + (NGLL3) * (ispec) - 0] - (1);
  atomicAdd(accel + (iglob) * (3) + 0, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 0 - 0]);
  atomicAdd(accel + (iglob) * (3) + 1, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 1 - 0]);
  atomicAdd(accel + (iglob) * (3) + 2, noise_sourcearray[(tx) * (3) + ((NGLL3) * (3)) * (it) + 2 - 0]);
}
