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
__global__ void get_maximum_vector_kernel(const float * array, const int size, float * d_max){
  __shared__ float sdata[BLOCKSIZE_TRANSFER + 0 - (1) - (0) + 1];
  int tid;
  int bx;
  int i;
  int s;
  tid = threadIdx.x;
  bx = (blockIdx.y) * (gridDim.x) + blockIdx.x;
  i = tid + (bx) * (blockDim.x);
  sdata[tid - 0] = (i < size ? sqrt((array[(i) * (3) + 0 - 0]) * (array[(i) * (3) + 0 - 0]) + (array[(i) * (3) + 1 - 0]) * (array[(i) * (3) + 1 - 0]) + (array[(i) * (3) + 2 - 0]) * (array[(i) * (3) + 2 - 0])) : 0.0f);
  __syncthreads();
  s = (blockDim.x) / (2);
  while(s > 0){
    if(tid < s){
      if(sdata[tid - 0] < sdata[tid + s - 0]){
        sdata[tid - 0] = sdata[tid + s - 0];
      }
    }
    s = s >> 1;
    __syncthreads();
  }
  if(tid == 0){
    d_max[bx - 0] = sdata[0 - 0];
  }
}
