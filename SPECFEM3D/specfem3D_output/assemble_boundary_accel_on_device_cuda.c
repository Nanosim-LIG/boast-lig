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
__global__ void assemble_boundary_accel_on_device(float * d_accel, const float * d_send_accel_buffer, const int num_interfaces, const int max_nibool_interfaces, const int * d_nibool_interfaces, const int * d_ibool_interfaces){
  int id;
  int iglob;
  int iloc;
  int iinterface;
  id = threadIdx.x + (blockIdx.x) * (blockDim.x) + ((gridDim.x) * (blockDim.x)) * (threadIdx.y + (blockIdx.y) * (blockDim.y));
  for(iinterface=0; iinterface<=num_interfaces - (1); iinterface+=1){
    if(id < d_nibool_interfaces[iinterface - 0]){
      iloc = id + (max_nibool_interfaces) * (iinterface);
      iglob = d_ibool_interfaces[iloc - 0] - (1);
      atomicAdd(d_accel + (iglob) * (3) + 0, d_send_accel_buffer[(iloc) * (3) + 0 - 0]);
      atomicAdd(d_accel + (iglob) * (3) + 1, d_send_accel_buffer[(iloc) * (3) + 1 - 0]);
      atomicAdd(d_accel + (iglob) * (3) + 2, d_send_accel_buffer[(iloc) * (3) + 2 - 0]);
    }
  }
}
