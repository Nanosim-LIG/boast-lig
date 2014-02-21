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
__global__ void noise_transfer_surface_to_host_kernel(const int * ibelm_top, const int nspec_top, const int * ibool, const float * displ, float * noise_surface_movie){
  int igll;
  int iface;
  igll = threadIdx.x;
  iface = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(iface < nspec_top){
    int i;
    int j;
    int k;
    int ispec;
    int iglob;
    ispec = ibelm_top[iface - 0] - (1);
    k = NGLLX - (1);
    j = (igll) / (NGLLX);
    i = igll - ((j) * (NGLLX));
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    noise_surface_movie[INDEX3(NDIM, NGLL2, 0, igll, iface) - 0] = displ[(iglob) * (3) - 0];
    noise_surface_movie[INDEX3(NDIM, NGLL2, 1, igll, iface) - 0] = displ[(iglob) * (3) - 0];
    noise_surface_movie[INDEX3(NDIM, NGLL2, 2, igll, iface) - 0] = displ[(iglob) * (3) - 0];
  }
}
