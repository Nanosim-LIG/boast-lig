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
__global__ void update_accel_acoustic_kernel(float * accel, const int size, const float * rmass){
  int id;
  id = threadIdx.x + (blockIdx.x) * (blockDim.x) + (blockIdx.y) * ((gridDim.x) * (blockDim.x));
  if(id < size){
    accel[id - 0] = (accel[id - 0]) * (rmass[id - 0]);
  }
}
