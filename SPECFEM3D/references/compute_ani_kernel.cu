// from compute_kernels_cuda.cu
__global__ void compute_ani_kernel(realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,realw* epsilondev_xz,realw* epsilondev_yz,realw* epsilon_trace_over_3,realw* b_epsilondev_xx,realw* b_epsilondev_yy,realw* b_epsilondev_xy,realw* b_epsilondev_xz,realw* b_epsilondev_yz,realw* b_epsilon_trace_over_3,realw* cijkl_kl,int NSPEC,realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;

  // handles case when there is 1 extra block (due to rectangular grid)
  if(ispec < NSPEC) {

    int ijk_ispec = threadIdx.x + NGLL3*ispec;

    // fully anisotropic kernel contributions
    realw eps_trace_over_3,b_eps_trace_over_3;
    realw prod[21];
    realw epsdev[5];
    realw b_epsdev[5];

    epsdev[0] = epsilondev_xx[ijk_ispec];
    epsdev[1] = epsilondev_yy[ijk_ispec];
    epsdev[2] = epsilondev_xy[ijk_ispec];
    epsdev[3] = epsilondev_xz[ijk_ispec];
    epsdev[4] = epsilondev_yz[ijk_ispec];

    b_epsdev[0] = b_epsilondev_xx[ijk_ispec];
    b_epsdev[1] = b_epsilondev_yy[ijk_ispec];
    b_epsdev[2] = b_epsilondev_xy[ijk_ispec];
    b_epsdev[3] = b_epsilondev_xz[ijk_ispec];
    b_epsdev[4] = b_epsilondev_yz[ijk_ispec];

    eps_trace_over_3 = epsilon_trace_over_3[ijk_ispec];
    b_eps_trace_over_3 = b_epsilon_trace_over_3[ijk_ispec];

    // fully anisotropic kernel contributions
    compute_strain_product_cuda(prod,eps_trace_over_3,epsdev,b_eps_trace_over_3,b_epsdev);

    // updates full anisotropic kernel
    for(int i=0;i<21;i++){
      cijkl_kl[i + 21*ijk_ispec] += deltat * prod[i];
    }
  }
}
