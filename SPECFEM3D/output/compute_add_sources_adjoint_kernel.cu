#ifndef INDEX2
#define INDEX2(isize,i,j) i + isize*j
#endif
#ifndef INDEX3
#define INDEX3(isize,jsize,i,j,k) i + isize*(j + jsize*k)
#endif
#ifndef INDEX4
#define INDEX4(isize,jsize,ksize,i,j,k,x) i + isize*(j + jsize*(k + ksize*x))
#endif
#ifndef INDEX5
#define INDEX5(isize,jsize,ksize,xsize,i,j,k,x,y) i + isize*(j + jsize*(k + ksize*(x + xsize*y)))
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
__global__ void compute_add_sources_adjoint_kernel(float * accel, const int nrec, const float * adj_sourcearrays, const int * ibool, const int * ispec_selected_rec, const int * pre_computed_irec, const int nadj_rec_local){
  int ispec;
  int iglob;
  int irec_local;
  int irec;
  int i;
  int j;
  int k;
  irec_local = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if (irec_local < nadj_rec_local) {
    irec = pre_computed_irec[irec_local];
    ispec = ispec_selected_rec[irec] - (1);
    i = threadIdx.x;
    j = threadIdx.y;
    k = threadIdx.z;
    iglob = ibool[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec)] - (1);
    atomicAdd(accel + (iglob) * (3) + 0, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 0, i, j, k, irec_local)]);
    atomicAdd(accel + (iglob) * (3) + 1, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 1, i, j, k, irec_local)]);
    atomicAdd(accel + (iglob) * (3) + 2, adj_sourcearrays[INDEX5(NDIM, NGLLX, NGLLX, NGLLX, 2, i, j, k, irec_local)]);
  }
}
