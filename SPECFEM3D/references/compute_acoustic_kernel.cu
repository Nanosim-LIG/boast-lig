// from compute_kernels_cuda.cu
__global__ void compute_acoustic_kernel(int* ibool,realw* rhostore,realw* kappastore,realw* hprime_xx,realw* d_xix,realw* d_xiy,realw* d_xiz,realw* d_etax,realw* d_etay,realw* d_etaz,realw* d_gammax,realw* d_gammay,realw* d_gammaz,realw* potential_dot_dot_acoustic,realw* b_potential_acoustic,realw* b_potential_dot_dot_acoustic,realw* rho_ac_kl,realw* kappa_ac_kl,realw deltat,int NSPEC) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if( ispec < NSPEC ){

    int ijk = threadIdx.x;

    // local and global indices
    int ijk_ispec = ijk + NGLL3*ispec;
    int ijk_ispec_padded = ijk + NGLL3_PADDED*ispec;
    int iglob = ibool[ijk_ispec] - 1;

    realw accel_elm[3];
    realw b_displ_elm[3];
    realw rhol,kappal;
    realw div_displ,b_div_displ;

    // shared memory between all threads within this block
    __shared__ realw scalar_field_displ[NGLL3];
    __shared__ realw scalar_field_accel[NGLL3];

    // copy field values
    scalar_field_displ[ijk] = b_potential_acoustic[iglob];
    scalar_field_accel[ijk] = potential_dot_dot_acoustic[iglob];
    __syncthreads();

    // displacement vector from backward field
    compute_gradient_kernel(ijk,ispec,scalar_field_displ,b_displ_elm,
                            hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz);

    // acceleration vector
    compute_gradient_kernel(ijk,ispec,scalar_field_accel,accel_elm,
                            hprime_xx,
                            d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz);

    // gets material parameter
    rhol = rhostore[ijk_ispec_padded];

    // density kernel
    rho_ac_kl[ijk_ispec] += deltat * rhol * (accel_elm[0]*b_displ_elm[0] +
                                             accel_elm[1]*b_displ_elm[1] +
                                             accel_elm[2]*b_displ_elm[2]);

    // bulk modulus kernel
    kappal = rhol/ kappastore[ijk_ispec_padded];

    div_displ = kappal * potential_dot_dot_acoustic[iglob];
    b_div_displ = kappal * b_potential_dot_dot_acoustic[iglob];

    kappa_ac_kl[ijk_ispec] += deltat * div_displ * b_div_displ;
  }
}
