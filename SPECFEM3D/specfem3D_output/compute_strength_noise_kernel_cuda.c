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
__global__ void compute_strength_noise_kernel(const float * displ, const int * ibelm_top, const int * ibool, const float * noise_surface_movie, const float * normal_x_noise, const float * normal_y_noise, const float * normal_z_noise, float * Sigma_kl, const float deltat, const int nspec_top){
  int iface;
  int ispec;
  int igll;
  int ipoin;
  int i;
  int j;
  int k;
  int iglob;
  float eta;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(iface < nspec_top){
    ispec = ibelm_top[iface - 0] - (1);
    igll = threadIdx.x;
    ipoin = igll + (NGLL2) * (iface);
    k = NGLLX - (1);
    j = (igll) / (NGLLX);
    i = igll - ((j) * (NGLLX));
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    eta = (noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0]) * (normal_x_noise[ipoin - 0]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0]) * (normal_y_noise[ipoin - 0]) + (noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0]) * (normal_z_noise[ipoin - 0]);
    Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] = Sigma_kl[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] + ((deltat) * (eta)) * ((normal_x_noise[ipoin - 0]) * (displ[0 - 0 + (iglob - (0)) * (3)]) + (normal_y_noise[ipoin - 0]) * (displ[1 - 0 + (iglob - (0)) * (3)]) + (normal_z_noise[ipoin - 0]) * (displ[2 - 0 + (iglob - (0)) * (3)]));
  }
}
