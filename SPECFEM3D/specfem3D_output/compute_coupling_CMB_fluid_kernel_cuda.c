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
__global__ void compute_coupling_CMB_fluid_kernel(const float * displ_crust_mantle, float * accel_crust_mantle, const float * accel_outer_core, const int * ibool_crust_mantle, const int * ibelm_bottom_crust_mantle, const float * normal_top_outer_core, const float * jacobian2D_top_outer_core, const float * wgllwgll_xy, const int * ibool_outer_core, const int * ibelm_top_outer_core, const float RHO_TOP_OC, const float minus_g_cmb, int GRAVITY, const int NSPEC2D_BOTTOM_CM){
  int i;
  int j;
  int k;
  int iface;
  int k_corresp;
  int iglob_oc;
  int iglob_cm;
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
  if(iface < NSPEC2D_BOTTOM_CM){
    ispec = ibelm_bottom_crust_mantle[iface - 0] - (1);
    ispec_selected = ibelm_top_outer_core[iface - 0] - (1);
    k = 0;
    k_corresp = NGLLX - (1);
    iglob_oc = ibool_outer_core[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k_corresp, ispec_selected) - 0] - (1);
    nx = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 0, i, j, iface) - 0];
    ny = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 1, i, j, iface) - 0];
    nz = normal_top_outer_core[INDEX4(NDIM, NGLLX, NGLLX, 2, i, j, iface) - 0];
    weight = (jacobian2D_top_outer_core[INDEX3(NGLLX, NGLLX, i, j, iface) - 0]) * (wgllwgll_xy[INDEX2(NGLLX, i, j) - 0]);
    iglob_cm = ibool_crust_mantle[INDEX4(NGLLX, NGLLX, NGLLX, i, j, k, ispec) - 0] - (1);
    if(GRAVITY){
      pressure = (RHO_TOP_OC) * ((minus_g_cmb) * ((displ_crust_mantle[(iglob_cm) * (3) - 0]) * (nx) + (displ_crust_mantle[(iglob_cm) * (3) + 1 - 0]) * (ny) + (displ_crust_mantle[(iglob_cm) * (3) + 2 - 0]) * (nz)) - (accel_outer_core[iglob_oc - 0]));
    } else {
      pressure = ( - (RHO_TOP_OC)) * (accel_outer_core[iglob_oc - 0]);
    }
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 0, ((weight) * (nx)) * (pressure));
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 1, ((weight) * (ny)) * (pressure));
    atomicAdd(accel_crust_mantle + (iglob_cm) * (3) + 2, ((weight) * (nz)) * (pressure));
  }
}