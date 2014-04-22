// from compute_forces_crust_mantle_cuda.cu
template<int FORWARD_OR_ADJOINT> __global__ void#ifdef USE_LAUNCH_BOUNDS__launch_bounds__(NGLL3_PADDED,LAUNCH_MIN_BLOCKS)#endifcrust_mantle_impl_kernel( int nb_blocks_to_compute,const int* d_ibool,const int* d_ispec_is_tiso,const int* d_phase_ispec_inner,int num_phase_ispec,const int d_iphase,realw deltat,const int use_mesh_coloring_gpu,realw_const_p d_displ,realw_p d_accel,realw_const_p d_xix, realw_const_p d_xiy, realw_const_p d_xiz,realw_const_p d_etax, realw_const_p d_etay, realw_const_p d_etaz,realw_const_p d_gammax, realw_const_p d_gammay, realw_const_p d_gammaz,realw_const_p d_hprime_xx,realw_const_p d_hprimewgll_xx,realw_const_p d_wgllwgll_xy,realw_const_p d_wgllwgll_xz,realw_const_p d_wgllwgll_yz,realw_const_p d_kappavstore,realw_const_p d_muvstore,realw_const_p d_kappahstore,realw_const_p d_muhstore,realw_const_p d_eta_anisostore,const int COMPUTE_AND_STORE_STRAIN,realw_p epsilondev_xx,realw_p epsilondev_yy,realw_p epsilondev_xy,realw_p epsilondev_xz,realw_p epsilondev_yz,realw_p epsilon_trace_over_3,const int ATTENUATION,const int PARTIAL_PHYS_DISPERSION_ONLY,const int USE_3D_ATTENUATION_ARRAYS,realw_const_p one_minus_sum_beta,realw_const_p factor_common,realw_p R_xx, realw_p R_yy, realw_p R_xy, realw_p R_xz, realw_p R_yz,realw_const_p alphaval,realw_const_p betaval,realw_const_p gammaval,const int ANISOTROPY,realw_const_p d_c11store,realw_const_p d_c12store,realw_const_p d_c13store,realw_const_p d_c14store,realw_const_p d_c15store,realw_const_p d_c16store,realw_const_p d_c22store,realw_const_p d_c23store,realw_const_p d_c24store,realw_const_p d_c25store,realw_const_p d_c26store,realw_const_p d_c33store,realw_const_p d_c34store,realw_const_p d_c35store,realw_const_p d_c36store,realw_const_p d_c44store,realw_const_p d_c45store,realw_const_p d_c46store,realw_const_p d_c55store,realw_const_p d_c56store,realw_const_p d_c66store,const int GRAVITY,realw_const_p d_xstore,realw_const_p d_ystore,realw_const_p d_zstore,realw_const_p d_minus_gravity_table,realw_const_p d_minus_deriv_gravity_table,realw_const_p d_density_table,realw_const_p wgll_cube,const int NSPEC_CRUST_MANTLE_STRAIN_ONLY ){

  // block id == spectral-element id
  int bx = blockIdx.y*gridDim.x+blockIdx.x;
  // thread id == GLL point id
  int tx = threadIdx.x;

  int K = (tx/NGLL2);
  int J = ((tx-K*NGLL2)/NGLLX);
  int I = (tx-K*NGLL2-J*NGLLX);

  unsigned short int active;
  int iglob,offset;
  int working_element;

  realw tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l;
  realw xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl;
  realw duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl;
  realw duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl;
  realw duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl;
  realw templ;

  realw fac1,fac2,fac3;
  realw one_minus_sum_beta_use;

  realw sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz;
  realw epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc;
  realw sum_terms1,sum_terms2,sum_terms3;

  // gravity variables
  realw sigma_yx,sigma_zx,sigma_zy;
  realw rho_s_H1,rho_s_H2,rho_s_H3;

#ifndef MANUALLY_UNROLLED_LOOPS
  int l;
#endif

  // shared memory arrays
  __shared__ realw s_dummyx_loc[NGLL3];
  __shared__ realw s_dummyy_loc[NGLL3];
  __shared__ realw s_dummyz_loc[NGLL3];

  __shared__ realw s_tempx1[NGLL3];
  __shared__ realw s_tempx2[NGLL3];
  __shared__ realw s_tempx3[NGLL3];

  __shared__ realw s_tempy1[NGLL3];
  __shared__ realw s_tempy2[NGLL3];
  __shared__ realw s_tempy3[NGLL3];

  __shared__ realw s_tempz1[NGLL3];
  __shared__ realw s_tempz2[NGLL3];
  __shared__ realw s_tempz3[NGLL3];

  // note: using shared memory for hprime's improves performance
  //       (but could tradeoff with occupancy)
  __shared__ realw sh_hprime_xx[NGLL2];
  __shared__ realw sh_hprimewgll_xx[NGLL2];

  // use only NGLL^3 = 125 active threads, plus 3 inactive/ghost threads,
  // because we used memory padding from NGLL^3 = 125 to 128 to get coalescent memory accesses
  active = (tx < NGLL3 && bx < nb_blocks_to_compute) ? 1:0;

  // determines spectral element to work on
  if( active ) {
#ifdef USE_MESH_COLORING_GPU
    working_element = bx;
#else
    //mesh coloring
    if( use_mesh_coloring_gpu ){
      working_element = bx;
    }else{
      // iphase-1 and working_element-1 for Fortran->C array conventions
      working_element = d_phase_ispec_inner[bx + num_phase_ispec*(d_iphase-1)]-1;
    }
#endif
    // local padded index
    offset = working_element*NGLL3_PADDED + tx;

    // global index
    iglob = d_ibool[working_element*NGLL3 + tx]-1;

    // copy displacement from global memory to shared memory
    load_shared_memory_cm<FORWARD_OR_ADJOINT>(&tx,&iglob,d_displ,s_dummyx_loc,s_dummyy_loc,s_dummyz_loc);
  } // active

  // loads hprime's into shared memory
  if( tx < NGLL2 ) {
    // copy hprime from global memory to shared memory
    load_shared_memory_hprime(&tx,d_hprime_xx,sh_hprime_xx);
    // copy hprime from global memory to shared memory
    load_shared_memory_hprimewgll(&tx,d_hprimewgll_xx,sh_hprimewgll_xx);
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if( active ) {

#ifndef MANUALLY_UNROLLED_LOOPS
    tempx1l = 0.f;
    tempx2l = 0.f;
    tempx3l = 0.f;

    tempy1l = 0.f;
    tempy2l = 0.f;
    tempy3l = 0.f;

    tempz1l = 0.f;
    tempz2l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
        fac1 = sh_hprime_xx[l*NGLLX+I];
        tempx1l += s_dummyx_loc[K*NGLL2+J*NGLLX+l]*fac1;
        tempy1l += s_dummyy_loc[K*NGLL2+J*NGLLX+l]*fac1;
        tempz1l += s_dummyz_loc[K*NGLL2+J*NGLLX+l]*fac1;

        fac2 = sh_hprime_xx[l*NGLLX+J];
        tempx2l += s_dummyx_loc[K*NGLL2+l*NGLLX+I]*fac2;
        tempy2l += s_dummyy_loc[K*NGLL2+l*NGLLX+I]*fac2;
        tempz2l += s_dummyz_loc[K*NGLL2+l*NGLLX+I]*fac2;

        fac3 = sh_hprime_xx[l*NGLLX+K];
        tempx3l += s_dummyx_loc[l*NGLL2+J*NGLLX+I]*fac3;
        tempy3l += s_dummyy_loc[l*NGLL2+J*NGLLX+I]*fac3;
        tempz3l += s_dummyz_loc[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else

    tempx1l = s_dummyx_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyx_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempy1l = s_dummyy_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyy_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempz1l = s_dummyz_loc[K*NGLL2+J*NGLLX]*sh_hprime_xx[I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+1]*sh_hprime_xx[NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+2]*sh_hprime_xx[2*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+3]*sh_hprime_xx[3*NGLLX+I]
            + s_dummyz_loc[K*NGLL2+J*NGLLX+4]*sh_hprime_xx[4*NGLLX+I];

    tempx2l = s_dummyx_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyx_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyx_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyx_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempy2l = s_dummyy_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyy_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyy_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyy_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempz2l = s_dummyz_loc[K*NGLL2+I]*sh_hprime_xx[J]
            + s_dummyz_loc[K*NGLL2+NGLLX+I]*sh_hprime_xx[NGLLX+J]
            + s_dummyz_loc[K*NGLL2+2*NGLLX+I]*sh_hprime_xx[2*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+3*NGLLX+I]*sh_hprime_xx[3*NGLLX+J]
            + s_dummyz_loc[K*NGLL2+4*NGLLX+I]*sh_hprime_xx[4*NGLLX+J];

    tempx3l = s_dummyx_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyx_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyx_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyx_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyx_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];

    tempy3l = s_dummyy_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyy_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyy_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyy_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyy_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];

    tempz3l = s_dummyz_loc[J*NGLLX+I]*sh_hprime_xx[K]
            + s_dummyz_loc[NGLL2+J*NGLLX+I]*sh_hprime_xx[NGLLX+K]
            + s_dummyz_loc[2*NGLL2+J*NGLLX+I]*sh_hprime_xx[2*NGLLX+K]
            + s_dummyz_loc[3*NGLL2+J*NGLLX+I]*sh_hprime_xx[3*NGLLX+K]
            + s_dummyz_loc[4*NGLL2+J*NGLLX+I]*sh_hprime_xx[4*NGLLX+K];


#endif

    // compute derivatives of ux, uy and uz with respect to x, y and z
    xixl = d_xix[offset];
    xiyl = d_xiy[offset];
    xizl = d_xiz[offset];
    etaxl = d_etax[offset];
    etayl = d_etay[offset];
    etazl = d_etaz[offset];
    gammaxl = d_gammax[offset];
    gammayl = d_gammay[offset];
    gammazl = d_gammaz[offset];

    duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l;
    duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l;
    duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l;

    duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l;
    duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l;
    duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l;

    duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l;
    duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l;
    duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l;

    // precompute some sums to save CPU time
    duxdxl_plus_duydyl = duxdxl + duydyl;
    duxdxl_plus_duzdzl = duxdxl + duzdzl;
    duydyl_plus_duzdzl = duydyl + duzdzl;
    duxdyl_plus_duydxl = duxdyl + duydxl;
    duzdxl_plus_duxdzl = duzdxl + duxdzl;
    duzdyl_plus_duydzl = duzdyl + duydzl;

    // computes deviatoric strain attenuation and/or for kernel calculations
    if(COMPUTE_AND_STORE_STRAIN) {
      templ = 0.33333333333333333333f * (duxdxl + duydyl + duzdzl); // 1./3. = 0.33333

      // local storage: stresses at this current time step
      epsilondev_xx_loc = duxdxl - templ;
      epsilondev_yy_loc = duydyl - templ;
      epsilondev_xy_loc = 0.5f * duxdyl_plus_duydxl;
      epsilondev_xz_loc = 0.5f * duzdxl_plus_duxdzl;
      epsilondev_yz_loc = 0.5f * duzdyl_plus_duydzl;

      if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) {
        epsilon_trace_over_3[tx] = templ;
      }else{
        epsilon_trace_over_3[tx + working_element*NGLL3] = templ;
      }
    }

    // attenuation
    if(ATTENUATION){
      // use unrelaxed parameters if attenuation
      if( USE_3D_ATTENUATION_ARRAYS ){
        one_minus_sum_beta_use = one_minus_sum_beta[tx+working_element*NGLL3]; // (i,j,k,ispec)
      }else{
        one_minus_sum_beta_use = one_minus_sum_beta[working_element]; // (1,1,1,ispec)
      }
    }

    // computes stresses
    if(ANISOTROPY){
      // full anisotropic case, stress calculations
      compute_element_cm_aniso(offset,
                            d_c11store,d_c12store,d_c13store,d_c14store,d_c15store,d_c16store,d_c22store,
                            d_c23store,d_c24store,d_c25store,d_c26store,d_c33store,d_c34store,d_c35store,
                            d_c36store,d_c44store,d_c45store,d_c46store,d_c55store,d_c56store,d_c66store,
                            ATTENUATION,
                            one_minus_sum_beta_use,
                            duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl,
                            duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                            &sigma_xx,&sigma_yy,&sigma_zz,
                            &sigma_xy,&sigma_xz,&sigma_yz);

    }else{
      if( ! d_ispec_is_tiso[working_element] ){
        // isotropic case
        compute_element_cm_iso(offset,
                            d_kappavstore,d_muvstore,
                            ATTENUATION,
                            one_minus_sum_beta_use,
                            duxdxl,duydyl,duzdzl,
                            duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl,
                            duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                            &sigma_xx,&sigma_yy,&sigma_zz,
                            &sigma_xy,&sigma_xz,&sigma_yz);
      }else{
        // transverse isotropy
        compute_element_cm_tiso(offset,
                              d_kappavstore,d_muvstore,
                              d_kappahstore,d_muhstore,d_eta_anisostore,
                              ATTENUATION,
                              one_minus_sum_beta_use,
                              duxdxl,duxdyl,duxdzl,
                              duydxl,duydyl,duydzl,
                              duzdxl,duzdyl,duzdzl,
                              duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl,
                              iglob,
                              d_ystore,d_zstore,
                              &sigma_xx,&sigma_yy,&sigma_zz,
                              &sigma_xy,&sigma_xz,&sigma_yz);
      }
    } // ! end of test whether isotropic or anisotropic element

    if(ATTENUATION && (! PARTIAL_PHYS_DISPERSION_ONLY ) ){
      // subtracts memory variables if attenuation
      compute_element_cm_att_stress(tx,working_element,
                                    R_xx,R_yy,R_xy,R_xz,R_yz,
                                    &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_xz,&sigma_yz);
    }

    // define symmetric components (needed for non-symmetric dot product and sigma for gravity)
    sigma_yx = sigma_xy;
    sigma_zx = sigma_xz;
    sigma_zy = sigma_yz;

    // jacobian
    jacobianl = 1.0f / (xixl*(etayl*gammazl-etazl*gammayl)
                      - xiyl*(etaxl*gammazl-etazl*gammaxl)
                      + xizl*(etaxl*gammayl-etayl*gammaxl));

    if( GRAVITY ){
      //  computes non-symmetric terms for gravity
      compute_element_cm_gravity(tx,iglob,
                                 d_xstore,d_ystore,d_zstore,
                                 d_minus_gravity_table,d_minus_deriv_gravity_table,d_density_table,
                                 wgll_cube,jacobianl,
                                 s_dummyx_loc,s_dummyy_loc,s_dummyz_loc,
                                 &sigma_xx,&sigma_yy,&sigma_zz,&sigma_xy,&sigma_yx,
                                 &sigma_xz,&sigma_zx,&sigma_yz,&sigma_zy,
                                 &rho_s_H1,&rho_s_H2,&rho_s_H3);
    }

    // form dot product with test vector, non-symmetric form
    s_tempx1[tx] = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl);
    s_tempy1[tx] = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl);
    s_tempz1[tx] = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl);

    s_tempx2[tx] = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl);
    s_tempy2[tx] = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl);
    s_tempz2[tx] = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl);

    s_tempx3[tx] = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl);
    s_tempy3[tx] = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl);
    s_tempz3[tx] = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl);
  }

// synchronize all the threads (one thread for each of the NGLL grid points of the
// current spectral element) because we need the whole element to be ready in order
// to be able to compute the matrix products along cut planes of the 3D element below
  __syncthreads();

  if( active ) {

#ifndef MANUALLY_UNROLLED_LOOPS
    tempx1l = 0.f;
    tempy1l = 0.f;
    tempz1l = 0.f;

    tempx2l = 0.f;
    tempy2l = 0.f;
    tempz2l = 0.f;

    tempx3l = 0.f;
    tempy3l = 0.f;
    tempz3l = 0.f;

    for (l=0;l<NGLLX;l++) {
      fac1 = sh_hprimewgll_xx[I*NGLLX+l];
      tempx1l += s_tempx1[K*NGLL2+J*NGLLX+l]*fac1;
      tempy1l += s_tempy1[K*NGLL2+J*NGLLX+l]*fac1;
      tempz1l += s_tempz1[K*NGLL2+J*NGLLX+l]*fac1;

      // assume hprimewgll_xx == hprimewgll_yy == hprimewgll_zz
      fac2 = sh_hprimewgll_xx[J*NGLLX+l];
      tempx2l += s_tempx2[K*NGLL2+l*NGLLX+I]*fac2;
      tempy2l += s_tempy2[K*NGLL2+l*NGLLX+I]*fac2;
      tempz2l += s_tempz2[K*NGLL2+l*NGLLX+I]*fac2;

      fac3 = sh_hprimewgll_xx[K*NGLLX+l];
      tempx3l += s_tempx3[l*NGLL2+J*NGLLX+I]*fac3;
      tempy3l += s_tempy3[l*NGLL2+J*NGLLX+I]*fac3;
      tempz3l += s_tempz3[l*NGLL2+J*NGLLX+I]*fac3;
    }
#else

    tempx1l = s_tempx1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempx1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempx1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempx1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempx1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempy1l = s_tempy1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempy1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempy1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempy1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempy1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempz1l = s_tempz1[K*NGLL2+J*NGLLX]*sh_hprimewgll_xx[I*NGLLX]
            + s_tempz1[K*NGLL2+J*NGLLX+1]*sh_hprimewgll_xx[I*NGLLX+1]
            + s_tempz1[K*NGLL2+J*NGLLX+2]*sh_hprimewgll_xx[I*NGLLX+2]
            + s_tempz1[K*NGLL2+J*NGLLX+3]*sh_hprimewgll_xx[I*NGLLX+3]
            + s_tempz1[K*NGLL2+J*NGLLX+4]*sh_hprimewgll_xx[I*NGLLX+4];

    tempx2l = s_tempx2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempx2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempx2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempx2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempx2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempy2l = s_tempy2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempy2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempy2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempy2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempy2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempz2l = s_tempz2[K*NGLL2+I]*sh_hprimewgll_xx[J*NGLLX]
            + s_tempz2[K*NGLL2+NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+1]
            + s_tempz2[K*NGLL2+2*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+2]
            + s_tempz2[K*NGLL2+3*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+3]
            + s_tempz2[K*NGLL2+4*NGLLX+I]*sh_hprimewgll_xx[J*NGLLX+4];

    tempx3l = s_tempx3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempx3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempx3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempx3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempx3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

    tempy3l = s_tempy3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempy3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempy3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempy3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempy3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

    tempz3l = s_tempz3[J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX]
            + s_tempz3[NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+1]
            + s_tempz3[2*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+2]
            + s_tempz3[3*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+3]
            + s_tempz3[4*NGLL2+J*NGLLX+I]*sh_hprimewgll_xx[K*NGLLX+4];

#endif

    fac1 = d_wgllwgll_yz[K*NGLLX+J];
    fac2 = d_wgllwgll_xz[K*NGLLX+I];
    fac3 = d_wgllwgll_xy[J*NGLLX+I];

    sum_terms1 = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l);
    sum_terms2 = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l);
    sum_terms3 = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l);

    // adds gravity term
    if( GRAVITY ){
      sum_terms1 += rho_s_H1;
      sum_terms2 += rho_s_H2;
      sum_terms3 += rho_s_H3;
    }

#ifdef USE_MESH_COLORING_GPU
    // no atomic operation needed, colors don't share global points between elements

#ifdef USE_TEXTURES_FIELDS
    d_accel[iglob*3]     = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
    d_accel[iglob*3 + 1] = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
    d_accel[iglob*3 + 2] = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
    d_accel[iglob*3]     += sum_terms1;
    d_accel[iglob*3 + 1] += sum_terms2;
    d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

#else // MESH_COLORING

    //mesh coloring
    if( use_mesh_coloring_gpu ){

      // no atomic operation needed, colors don't share global points between elements
#ifdef USE_TEXTURES_FIELDS
      d_accel[iglob*3]     = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3) + sum_terms1;
      d_accel[iglob*3 + 1] = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3 + 1) + sum_terms2;
      d_accel[iglob*3 + 2] = texfetch_accel_cm<FORWARD_OR_ADJOINT>(iglob*3 + 2) + sum_terms3;
#else
      d_accel[iglob*3]     += sum_terms1;
      d_accel[iglob*3 + 1] += sum_terms2;
      d_accel[iglob*3 + 2] += sum_terms3;
#endif // USE_TEXTURES_FIELDS

    }else{
      // no mesh coloring uses atomic updates

      atomicAdd(&d_accel[iglob*3], sum_terms1);
      atomicAdd(&d_accel[iglob*3 + 1], sum_terms2);
      atomicAdd(&d_accel[iglob*3 + 2], sum_terms3);

      // debug: for testing purposes only: w/out atomic updates
      //d_accel[iglob*3] -= (0.00000001f*tempx1l + 0.00000001f*tempx2l + 0.00000001f*tempx3l);
      //d_accel[iglob*3 + 1] -= (0.00000001f*tempy1l + 0.00000001f*tempy2l + 0.00000001f*tempy3l);
      //d_accel[iglob*3 + 2] -= (0.00000001f*tempz1l + 0.00000001f*tempz2l + 0.00000001f*tempz3l);
    }
#endif // MESH_COLORING

    // update memory variables based upon the Runge-Kutta scheme
    if( ATTENUATION && ( ! PARTIAL_PHYS_DISPERSION_ONLY ) ){
      compute_element_cm_att_memory(tx,working_element,
                                    d_muvstore,
                                    factor_common,alphaval,betaval,gammaval,
                                    R_xx,R_yy,R_xy,R_xz,R_yz,
                                    epsilondev_xx,epsilondev_yy,epsilondev_xy,
                                    epsilondev_xz,epsilondev_yz,
                                    epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,
                                    epsilondev_xz_loc,epsilondev_yz_loc,
                                    d_c44store,ANISOTROPY,USE_3D_ATTENUATION_ARRAYS);
    }

    // save deviatoric strain for Runge-Kutta scheme
    if( COMPUTE_AND_STORE_STRAIN ){
      // fortran: epsilondev_xx(:,:,:,ispec) = epsilondev_xx_loc(:,:,:)
      epsilondev_xx[tx + working_element*NGLL3] = epsilondev_xx_loc;
      epsilondev_yy[tx + working_element*NGLL3] = epsilondev_yy_loc;
      epsilondev_xy[tx + working_element*NGLL3] = epsilondev_xy_loc;
      epsilondev_xz[tx + working_element*NGLL3] = epsilondev_xz_loc;
      epsilondev_yz[tx + working_element*NGLL3] = epsilondev_yz_loc;
    }
  } // active
}
