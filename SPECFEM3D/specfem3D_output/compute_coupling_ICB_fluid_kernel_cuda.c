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
__global__ void compute_coupling_ICB_fluid_kernel(const float * displ_inner_core, float * accel_inner_core, const float * accel_outer_core, const int * ibool_inner_core, const int * ibelm_top_inner_core, const float * normal_bottom_outer_core, const float * jacobian2D_bottom_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_bottom_outer_core, const float RHO_BOTTOM_OC, const float minus_g_icb, int GRAVITY, const int NSPEC2D_TOP_IC){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  int iglob_oc;
  int iglob_ic;
  float pressure;
  int ispec;
  int ispec_selected;
  float displ_n;
  float nx;
  float ny;
  float nz;
  float weight;
  i = threadIdx.x;
  j = threadIdx.y;
  iface = blockIdx.x + (gridDim.x) * (blockIdx.y);
  if(iface < NSPEC2D_TOP_IC){
    ispec = ibelm_top_inner_core[iface - 0] - (1);
    ispec_selected = ibelm_bottom_outer_core[iface - 0] - (1);
    k = NGLLX - (1);
    k_corresp = 0;
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);
    nx = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];
    ny = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];
    nz = normal_bottom_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];
    weight = (jacobian2D_bottom_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);
    iglob_ic = ibool_inner_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    if(GRAVITY){
      pressure = (RHO_BOTTOM_OC) * ((minus_g_icb) * ((displ_inner_core[(iglob_ic) * (3) - 0]) * (nx) + (displ_inner_core[(iglob_ic) * (3) + 1 - 0]) * (ny) + (displ_inner_core[(iglob_ic) * (3) + 2 - 0]) * (nz)) - (accel_outer_core[iglob_oc - 0]));
    } else {
      pressure = ( - (RHO_BOTTOM_OC)) * (accel_outer_core[iglob_oc - 0]);
    }
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 0, (( - (weight)) * (nx)) * (pressure));
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 1, (( - (weight)) * (ny)) * (pressure));
    atomicAdd(accel_inner_core + (iglob_ic) * (3) + 2, (( - (weight)) * (nz)) * (pressure));
  }
}
