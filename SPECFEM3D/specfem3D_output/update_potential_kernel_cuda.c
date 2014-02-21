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
__global__ void update_potential_kernel(float * potential_acoustic, float * potential_dot_acoustic, float * potential_dot_dot_acoustic, const int size, const float deltat, const float deltatsqover2, const float deltatover2){
  int id;
  id = threadIdx.x + (blockIdx.x) * (blockDim.x) + (blockIdx.y) * ((gridDim.x) * (blockDim.x));
  if(id < size){
    potential_acoustic[id - 0] = potential_acoustic[id - 0] + (deltat) * (potential_dot_acoustic[id - 0]) + (deltatsqover2) * (potential_dot_dot_acoustic[id - 0]);
    potential_dot_acoustic[id - 0] = potential_dot_acoustic[id - 0] + (deltatover2) * (potential_dot_dot_acoustic[id - 0]);
    potential_dot_dot_acoustic[id - 0] = 0.0f;
  }
}
