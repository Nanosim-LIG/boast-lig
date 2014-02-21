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
__global__ void compute_add_sources_kernel(float * accel, const int * ibool, const float * sourcearrays, const double * stf_pre_compute, const int myrank, const int * islice_selected_source, const int * ispec_selected_source, const int nsources){
  int ispec;
  int iglob;
  float stf;
  int isource;
  int i;
  int j;
  int k;
  i = threadIdx.x;
  j = threadIdx.y;
  k = threadIdx.z;
  isource = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if(isource < nsources){
    if(myrank == islice_selected_source[isource - 0]){
      ispec = ispec_selected_source[isource - 0] - (1);
      stf = stf_pre_compute[isource - 0];
      iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
      atomicAdd(accel + (iglob) * (3) + 0, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 0, i, j, k, isource) - 0]) * (stf));
      atomicAdd(accel + (iglob) * (3) + 1, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 1, i, j, k, isource) - 0]) * (stf));
      atomicAdd(accel + (iglob) * (3) + 2, (sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 2, i, j, k, isource) - 0]) * (stf));
    }
  }
}
