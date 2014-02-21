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
__global__ void compute_add_sources_adjoint_kernel(const int nrec, float * accel, const float * adj_sourcearrays, const int * ibool, const int * ispec_selected_rec, const int * pre_computed_irec, const int nadj_rec_local){
  int ispec;
  int iglob;
  int irec_local;
  int irec;
  int i;
  int j;
  int k;
  irec_local = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if(irec_local < nadj_rec_local){
    irec = pre_computed_irec[irec_local - 0];
    ispec = ispec_selected_rec[irec - 0] - (1);
    i = threadIdx.x;
    j = threadIdx.y;
    k = threadIdx.z;
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    atomicAdd(accel + (iglob) * (3) + 0, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 0, i, j, k, irec_local) - 0]);
    atomicAdd(accel + (iglob) * (3) + 1, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 1, i, j, k, irec_local) - 0]);
    atomicAdd(accel + (iglob) * (3) + 2, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 2, i, j, k, irec_local) - 0]);
  }
}
