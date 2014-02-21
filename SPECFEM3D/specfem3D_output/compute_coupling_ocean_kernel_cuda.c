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
__global__ void compute_coupling_ocean_kernel(float * accel_crust_mantle, const float * rmassx_crust_mantle, const float * rmassy_crust_mantle, const float * rmassz_crust_mantle, const float * rmass_ocean_load, const int npoin_ocean_load, const int * ibool_ocean_load, const float * normal_ocean_load){
  int ipoin;
  int iglob;
  float nx;
  float ny;
  float nz;
  float rmass;
  float force_normal_comp;
  float additional_term_x;
  float additional_term_y;
  float additional_term_z;
  ipoin = threadIdx.x + (blockIdx.x) * (blockDim.x) + ((gridDim.x) * (blockDim.x)) * (threadIdx.y + (blockIdx.y) * (blockDim.y));
  if(ipoin < npoin_ocean_load){
    iglob = ibool_ocean_load[ipoin - 0] - (1);
    nx = normal_ocean_load[INDEX2(NDIM, 0, ipoin) - 0];
    ny = normal_ocean_load[INDEX2(NDIM, 1, ipoin) - 0];
    nz = normal_ocean_load[INDEX2(NDIM, 2, ipoin) - 0];
    force_normal_comp = ((accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)]) * (nx)) / (rmassx_crust_mantle[iglob - 0]) + ((accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)]) * (ny)) / (rmassy_crust_mantle[iglob - 0]) + ((accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)]) * (nz)) / (rmassz_crust_mantle[iglob - 0]);
    rmass = rmass_ocean_load[ipoin - 0];
    additional_term_x = (rmass - (rmassx_crust_mantle[iglob - 0])) * (force_normal_comp);
    additional_term_y = (rmass - (rmassy_crust_mantle[iglob - 0])) * (force_normal_comp);
    additional_term_z = (rmass - (rmassz_crust_mantle[iglob - 0])) * (force_normal_comp);
    accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[0 - 0 + (iglob - (0)) * (3)] + (additional_term_x) * (nx);
    accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[1 - 0 + (iglob - (0)) * (3)] + (additional_term_y) * (ny);
    accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)] = accel_crust_mantle[2 - 0 + (iglob - (0)) * (3)] + (additional_term_z) * (nz);
  }
}
