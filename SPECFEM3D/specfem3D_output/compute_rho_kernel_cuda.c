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
__global__ void compute_rho_kernel(const int * ibool, const float * accel, const float * b_displ, float * rho_kl, const int NSPEC, const float deltat){
  int ispec;
  int ijk_ispec;
  int iglob;
  ispec = blockIdx.x + (blockIdx.y) * (gridDim.x);
  if(ispec < NSPEC){
    ijk_ispec = threadIdx.x + (NGLL3) * (ispec);
    iglob = ibool[ijk_ispec - 0] - (1);
    rho_kl[ijk_ispec - 0] = rho_kl[ijk_ispec - 0] + (deltat) * ((accel[0 - 0 + (iglob - (0)) * (3)]) * (b_displ[0 - 0 + (iglob - (0)) * (3)]) + (accel[1 - 0 + (iglob - (0)) * (3)]) * (b_displ[1 - 0 + (iglob - (0)) * (3)]) + (accel[2 - 0 + (iglob - (0)) * (3)]) * (b_displ[2 - 0 + (iglob - (0)) * (3)]));
  }
}
