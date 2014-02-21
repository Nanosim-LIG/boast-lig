char * outer_core_impl_kernel_program = "\
static inline void atomicAdd(volatile __global float *source, const float val) {\n\
  union {\n\
    unsigned int iVal;\n\
    float fVal;\n\
  } res, orig;\n\
  do {\n\
    orig.fVal = *source;\n\
    res.fVal = orig.fVal + val;\n\
  } while (atomic_cmpxchg((volatile __global unsigned int *)source, orig.iVal, res.iVal) != orig.iVal);\n\
}\n\
#define INDEX2(xsize,x,y) x + (y)*xsize\n\
#define INDEX3(xsize,ysize,x,y,z) x + xsize*(y + ysize*z)\n\
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))\n\
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))\n\
#define NDIM 3\n\
#define NGLLX 5\n\
#define NGLL2 25\n\
#define NGLL3 125\n\
#define NGLL3_PADDED 128\n\
#define N_SLS 3\n\
#define IREGION_CRUST_MANTLE 1\n\
#define IREGION_INNER_CORE 3\n\
#define IFLAG_IN_FICTITIOUS_CUBE 11\n\
#define R_EARTH_KM 6371.0f\n\
#define COLORING_MIN_NSPEC_INNER_CORE 1000\n\
#define COLORING_MIN_NSPEC_OUTER_CORE 1000\n\
#define BLOCKSIZE_TRANSFER 256\n\
void compute_element_oc_rotation(const int tx, const int working_element, const float time, const float two_omega_earth, const float deltat, __global float * d_A_array_rotation, __global float * d_B_array_rotation, const float dpotentialdxl, const float dpotentialdyl, float * dpotentialdx_with_rot, float * dpotentialdy_with_rot){\n\
  float two_omega_deltat;\n\
  float cos_two_omega_t;\n\
  float sin_two_omega_t;\n\
  float A_rotation;\n\
  float b_rotation;\n\
  float source_euler_A;\n\
  float source_euler_B;\n\
  cos_two_omega_t = cos((two_omega_earth) * (time));\n\
  sin_two_omega_t = sin((two_omega_earth) * (time));\n\
  two_omega_deltat = (deltat) * (two_omega_earth);\n\
  source_euler_A = (two_omega_deltat) * ((cos_two_omega_t) * (dpotentialdyl) + (sin_two_omega_t) * (dpotentialdxl));\n\
  source_euler_B = (two_omega_deltat) * ((sin_two_omega_t) * (dpotentialdyl) - ((cos_two_omega_t) * (dpotentialdxl)));\n\
  A_rotation = d_A_array_rotation[tx + (working_element) * (NGLL3) - 0];\n\
  b_rotation = d_B_array_rotation[tx + (working_element) * (NGLL3) - 0];\n\
  dpotentialdx_with_rot[0 - 0] = dpotentialdxl + (A_rotation) * (cos_two_omega_t) + (b_rotation) * (sin_two_omega_t);\n\
  dpotentialdy_with_rot[0 - 0] = dpotentialdyl + ( - (A_rotation)) * (sin_two_omega_t) + (b_rotation) * (cos_two_omega_t);\n\
  d_A_array_rotation[tx + (working_element) * (NGLL3) - 0] = d_A_array_rotation[tx + (working_element) * (NGLL3) - 0] + source_euler_A;\n\
  d_B_array_rotation[tx + (working_element) * (NGLL3) - 0] = d_B_array_rotation[tx + (working_element) * (NGLL3) - 0] + source_euler_B;\n\
}\n\
__kernel void outer_core_impl_kernel(const int nb_blocks_to_compute, const int NGLOB, const __global int * d_ibool, const __global int * phase_ispec_inner, const int num_phase_ispec, const int d_iphase, const int use_mesh_coloring_gpu, const __global float * d_potential, __global float * d_potential_dot_dot, const __global float * d_xix, const __global float * d_xiy, const __global float * d_xiz, const __global float * d_etax, const __global float * d_etay, const __global float * d_etaz, const __global float * d_gammax, const __global float * d_gammay, const __global float * d_gammaz, const __global float * d_hprime_xx, const __global float * d_hprimewgll_xx, const __global float * wgllwgll_xy, const __global float * wgllwgll_xz, const __global float * wgllwgll_yz, const int GRAVITY, const __global float * d_xstore, const __global float * d_ystore, const __global float * d_zstore, const __global float * d_d_ln_density_dr_table, const __global float * d_minus_rho_g_over_kappa_fluid, const __global float * wgll_cube, const int ROTATION, const float time, const float two_omega_earth, const float deltat, __global float * d_A_array_rotation, __global float * d_B_array_rotation, const int NSPEC_OUTER_CORE){\n\
  int bx;\n\
  int tx;\n\
  int K;\n\
  int J;\n\
  int I;\n\
  int active;\n\
  int offset;\n\
  int iglob;\n\
  int working_element;\n\
  int l;\n\
  float temp1l;\n\
  float temp2l;\n\
  float temp3l;\n\
  float xixl;\n\
  float xiyl;\n\
  float xizl;\n\
  float etaxl;\n\
  float etayl;\n\
  float etazl;\n\
  float gammaxl;\n\
  float gammayl;\n\
  float gammazl;\n\
  float jacobianl;\n\
  float dpotentialdxl;\n\
  float dpotentialdyl;\n\
  float dpotentialdzl;\n\
  float dpotentialdx_with_rot;\n\
  float dpotentialdy_with_rot;\n\
  float sum_terms;\n\
  float gravity_term;\n\
  float gxl;\n\
  float gyl;\n\
  float gzl;\n\
  float radius;\n\
  float theta;\n\
  float phi;\n\
  float cos_theta;\n\
  float sin_theta;\n\
  float cos_phi;\n\
  float sin_phi;\n\
  float grad_x_ln_rho;\n\
  float grad_y_ln_rho;\n\
  float grad_z_ln_rho;\n\
  int int_radius;\n\
  __local float s_dummy_loc[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp1[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp2[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float s_temp3[NGLL3 + 0 - (1) - (0) + 1];\n\
  __local float sh_hprime_xx[NGLL2 + 0 - (1) - (0) + 1];\n\
  __local float sh_hprimewgll_xx[NGLL2 + 0 - (1) - (0) + 1];\n\
  bx = (get_group_id(1)) * (get_num_groups(0)) + get_group_id(0);\n\
  tx = get_local_id(0);\n\
  K = (tx) / (NGLL2);\n\
  J = (tx - ((K) * (NGLL2))) / (NGLLX);\n\
  I = tx - ((K) * (NGLL2)) - ((J) * (NGLLX));\n\
  active = (tx < NGLL3 && bx < nb_blocks_to_compute ? 1 : 0);\n\
  if(active){\n\
    if(use_mesh_coloring_gpu){\n\
      working_element = bx;\n\
    } else {\n\
      working_element = phase_ispec_inner[bx + (num_phase_ispec) * (d_iphase - (1)) - 0] - (1);\n\
    }\n\
    iglob = d_ibool[(working_element) * (NGLL3) + tx - 0] - (1);\n\
    s_dummy_loc[tx - 0] = d_potential[iglob - 0];\n\
  }\n\
  if(tx < NGLL2){\n\
    sh_hprime_xx[tx - 0] = d_hprime_xx[tx - 0];\n\
    sh_hprimewgll_xx[tx - 0] = d_hprimewgll_xx[tx - 0];\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    temp1l = 0.0f;\n\
    temp2l = 0.0f;\n\
    temp3l = 0.0f;\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      temp1l = temp1l + (s_dummy_loc[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (sh_hprime_xx[(l) * (NGLLX) + I - 0]);\n\
      temp2l = temp2l + (s_dummy_loc[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (sh_hprime_xx[(l) * (NGLLX) + J - 0]);\n\
      temp3l = temp3l + (s_dummy_loc[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprime_xx[(l) * (NGLLX) + K - 0]);\n\
    }\n\
    offset = (working_element) * (NGLL3_PADDED) + tx;\n\
    xixl = d_xix[offset - 0];\n\
    etaxl = d_etax[offset - 0];\n\
    gammaxl = d_gammax[offset - 0];\n\
    xiyl = d_xiy[offset - 0];\n\
    etayl = d_etay[offset - 0];\n\
    gammayl = d_gammay[offset - 0];\n\
    xizl = d_xiz[offset - 0];\n\
    etazl = d_etaz[offset - 0];\n\
    gammazl = d_gammaz[offset - 0];\n\
    jacobianl = (1.0f) / ((xixl) * ((etayl) * (gammazl) - ((etazl) * (gammayl))) - ((xiyl) * ((etaxl) * (gammazl) - ((etazl) * (gammaxl)))) + (xizl) * ((etaxl) * (gammayl) - ((etayl) * (gammaxl))));\n\
    dpotentialdxl = (xixl) * (temp1l) + (etaxl) * (temp2l) + (gammaxl) * (temp3l);\n\
    dpotentialdyl = (xiyl) * (temp1l) + (etayl) * (temp2l) + (gammayl) * (temp3l);\n\
    dpotentialdzl = (xizl) * (temp1l) + (etazl) * (temp2l) + (gammazl) * (temp3l);\n\
    if(ROTATION){\n\
      compute_element_oc_rotation(tx, working_element, time, two_omega_earth, deltat, d_A_array_rotation, d_B_array_rotation, dpotentialdxl, dpotentialdyl,  &dpotentialdx_with_rot,  &dpotentialdy_with_rot);\n\
    } else {\n\
      dpotentialdx_with_rot = dpotentialdxl;\n\
      dpotentialdy_with_rot = dpotentialdyl;\n\
    }\n\
    radius = d_xstore[iglob - 0];\n\
    theta = d_ystore[iglob - 0];\n\
    phi = d_zstore[iglob - 0];\n\
    sin_theta = sincos(theta,  &cos_theta);\n\
    sin_phi = sincos(phi,  &cos_phi);\n\
    int_radius = rint(((radius) * (6371.0f)) * (10.0f)) - (1);\n\
    s_temp1[tx - 0] = (jacobianl) * ((xixl) * (dpotentialdx_with_rot) + (xiyl) * (dpotentialdy_with_rot) + (xizl) * (dpotentialdzl));\n\
    s_temp1[tx - 0] = (jacobianl) * ((etaxl) * (dpotentialdx_with_rot) + (etayl) * (dpotentialdy_with_rot) + (etazl) * (dpotentialdzl));\n\
    s_temp1[tx - 0] = (jacobianl) * ((gammaxl) * (dpotentialdx_with_rot) + (gammayl) * (dpotentialdy_with_rot) + (gammazl) * (dpotentialdzl));\n\
  }\n\
  barrier(CLK_LOCAL_MEM_FENCE);\n\
  if(active){\n\
    temp1l = 0.0f;\n\
    temp2l = 0.0f;\n\
    temp3l = 0.0f;\n\
    for(l=0; l<=NGLLX - (1); l+=1){\n\
      temp1l = temp1l + (s_temp1[(K) * (NGLL2) + (J) * (NGLLX) + l - 0]) * (sh_hprimewgll_xx[(I) * (NGLLX) + l - 0]);\n\
      temp2l = temp2l + (s_temp2[(K) * (NGLL2) + (l) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(J) * (NGLLX) + l - 0]);\n\
      temp3l = temp3l + (s_temp3[(l) * (NGLL2) + (J) * (NGLLX) + I - 0]) * (sh_hprimewgll_xx[(K) * (NGLLX) + l - 0]);\n\
    }\n\
    sum_terms =  - ((wgllwgll_yz[(K) * (NGLLX) + J - 0]) * (temp1l) + (wgllwgll_xz[(K) * (NGLLX) + I - 0]) * (temp2l) + (wgllwgll_xy[(J) * (NGLLX) + I - 0]) * (temp3l));\n\
    if(GRAVITY){\n\
      sum_terms = sum_terms + gravity_term;\n\
    }\n\
    if(use_mesh_coloring_gpu){\n\
      if(NSPEC_OUTER_CORE > 1000){\n\
        d_potential_dot_dot[iglob - 0] = d_potential_dot_dot[iglob - 0] + sum_terms;\n\
      } else {\n\
        atomicAdd(d_potential_dot_dot + iglob, sum_terms);\n\
      }\n\
    } else {\n\
      atomicAdd(d_potential_dot_dot + iglob, sum_terms);\n\
    }\n\
  }\n\
}\n\
";