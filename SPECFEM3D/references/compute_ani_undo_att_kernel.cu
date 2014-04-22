// from compute_kernels_cuda.cu
__global__ void compute_ani_undo_att_kernel(realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,realw* epsilondev_xz,realw* epsilondev_yz,realw* epsilon_trace_over_3,realw* cijkl_kl,int NSPEC,realw deltat,int* d_ibool,realw* d_b_displ,realw* d_xix,realw* d_xiy,realw* d_xiz,realw* d_etax,realw* d_etay,realw* d_etaz,realw* d_gammax,realw* d_gammay,realw* d_gammaz,realw* d_hprime_xx) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ijk_ispec = threadIdx.x + NGLL3*ispec;

  int tx = threadIdx.x;
  int iglob;

  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw sh_hprime_xx[NGLL2];

  // loads element displacements
  // all threads load their displacement into shared memory
  if( ispec < NSPEC){
    iglob = d_ibool[ijk_ispec]-1;
    // changing iglob indexing to match fortran row changes fast style
    s_dummyx_loc[tx] = d_b_displ[iglob*3];
    s_dummyy_loc[tx] = d_b_displ[iglob*3 + 1];
    s_dummyz_loc[tx] = d_b_displ[iglob*3 + 2];

    // master thread loads hprime
    if( threadIdx.x == 0 ){
      for(int m=0; m < NGLL2; m++){
        // hprime
        sh_hprime_xx[m] = d_hprime_xx[m];
      }
    }
  }

  // synchronizes threads
  __syncthreads();

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    // fully anisotropic kernel contributions
    realw eps_trace_over_3,b_eps_trace_over_3;
    realw prod[21];
    realw epsdev[5];
    realw b_epsdev[5];

    // strain from adjoint wavefield
    epsdev[0] = epsilondev_xx[ijk_ispec];
    epsdev[1] = epsilondev_yy[ijk_ispec];
    epsdev[2] = epsilondev_xy[ijk_ispec];
    epsdev[3] = epsilondev_xz[ijk_ispec];
    epsdev[4] = epsilondev_yz[ijk_ispec];
    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];

    // strain from backward/reconstructed forward wavefield
    compute_element_strain_undo_att(ispec,ijk_ispec,
                                    d_ibool,
                                    s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                                    d_xix,d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gammax,d_gammay,d_gammaz,
                                    sh_hprime_xx,
                                    b_epsdev,&b_eps_trace_over_3);

    // fully anisotropic kernel contributions
    compute_strain_product_cuda(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev);

    // updates full anisotropic kernel
    for(int i=0;i<21;i++){
      cijkl_kl[i + 21*ijk_ispec] += deltat * prod[i];
    }
  }
}
