#ifndef INDEX2
#define INDEX2(xsize,x,y) x + (y)*xsize
#endif
#ifndef INDEX3
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)
#endif
#ifndef INDEX4
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#endif
#ifndef INDEX5
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
#endif
#ifndef NDIM
#define NDIM 3
#endif
#ifndef NGLLX
#define NGLLX 5
#endif
#ifndef NGLL2
#define NGLL2 25
#endif
#ifndef NGLL3
#define NGLL3 125
#endif
#ifndef NGLL3_PADDED
#define NGLL3_PADDED 128
#endif
#ifndef N_SLS
#define N_SLS 3
#endif
#ifndef IREGION_CRUST_MANTLE
#define IREGION_CRUST_MANTLE 1
#endif
#ifndef IREGION_INNER_CORE
#define IREGION_INNER_CORE 3
#endif
#ifndef IFLAG_IN_FICTITIOUS_CUBE
#define IFLAG_IN_FICTITIOUS_CUBE 11
#endif
#ifndef R_EARTH_KM
#define R_EARTH_KM 6371.0f
#endif
#ifndef COLORING_MIN_NSPEC_INNER_CORE
#define COLORING_MIN_NSPEC_INNER_CORE 1000
#endif
#ifndef COLORING_MIN_NSPEC_OUTER_CORE
#define COLORING_MIN_NSPEC_OUTER_CORE 1000
#endif
#ifndef BLOCKSIZE_TRANSFER
#define BLOCKSIZE_TRANSFER 256
#endif
__global__ void prepare_boundary_accel_on_device(const float * d_accel, float * d_send_accel_buffer, const int num_interfaces, const int max_nibool_interfaces, const int * d_nibool_interfaces, const int * d_ibool_interfaces){
  int id;
  int iglob;
  int iloc;
  int iinterface;
  id = threadIdx.x + (blockIdx.x) * (blockDim.x) + ((gridDim.x) * (blockDim.x)) * (threadIdx.y + (blockIdx.y) * (blockDim.y));
  for(iinterface=0; iinterface<=num_interfaces - (1); iinterface+=1){
    if(id < d_nibool_interfaces[iinterface - 0]){
      iloc = id + (max_nibool_interfaces) * (iinterface);
      iglob = d_ibool_interfaces[iloc - 0] - (1);
      d_send_accel_buffer[(iloc) * (3) + 0 - 0] = d_accel[(iglob) * (3) + 0 - 0];
      d_send_accel_buffer[(iloc) * (3) + 1 - 0] = d_accel[(iglob) * (3) + 1 - 0];
      d_send_accel_buffer[(iloc) * (3) + 2 - 0] = d_accel[(iglob) * (3) + 2 - 0];
    }
  }
}