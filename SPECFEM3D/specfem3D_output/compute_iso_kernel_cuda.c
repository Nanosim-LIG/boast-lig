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
__global__ void compute_iso_kernel(const float * epsilondev_xx, const float * epsilondev_yy, const float * epsilondev_xy, const float * epsilondev_xz, const float * epsilondev_yz, const float * epsilon_trace_over_3, const float * b_epsilondev_xx, const float * b_epsilondev_yy, const float * b_epsilondev_xy, const float * b_epsilondev_xz, const float * b_epsilondev_yz, const float * b_epsilon_trace_over_3, float * mu_kl, float * kappa_kl, const int NSPEC, const float deltat){
  int ispec;
  int ijk_ispec;
  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(ispec < NSPEC){
    ijk_ispec = threadIdx.x + (NGLL3) * (ispec);
    mu_kl[ijk_ispec - 0] = mu_kl[ijk_ispec - 0] + (deltat) * ((epsilondev_xx[ijk_ispec - 0]) * (b_epsilondev_xx[ijk_ispec - 0]) + (epsilondev_yy[ijk_ispec - 0]) * (b_epsilondev_yy[ijk_ispec - 0]) + (epsilondev_xx[ijk_ispec - 0] + epsilondev_yy[ijk_ispec - 0]) * (b_epsilondev_xx[ijk_ispec - 0] + b_epsilondev_yy[ijk_ispec - 0]) + ((epsilondev_xy[ijk_ispec - 0]) * (b_epsilondev_xy[ijk_ispec - 0]) + (epsilondev_xz[ijk_ispec - 0]) * (b_epsilondev_xz[ijk_ispec - 0]) + (epsilondev_yz[ijk_ispec - 0]) * (b_epsilondev_yz[ijk_ispec - 0])) * (2));
    kappa_kl[ijk_ispec - 0] = kappa_kl[ijk_ispec - 0] + (deltat) * (((epsilon_trace_over_3[ijk_ispec - 0]) * (b_epsilon_trace_over_3[ijk_ispec - 0])) * (9));
  }
}
