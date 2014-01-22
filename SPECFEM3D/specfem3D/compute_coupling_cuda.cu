/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  2 . 0
 !               ---------------------------------------
 !
 !          Main authors: Dimitri Komatitsch and Jeroen Tromp
 !    Princeton University, USA and University of Pau / CNRS / INRIA
 ! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
 !                            August 2013
 !
 ! This program is free software; you can redistribute it and/or modify
 ! it under the terms of the GNU General Public License as published by
 ! the Free Software Foundation; either version 2 of the License, or
 ! (at your option) any later version.
 !
 ! This program is distributed in the hope that it will be useful,
 ! but WITHOUT ANY WARRANTY; without even the implied warranty of
 ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ! GNU General Public License for more details.
 !
 ! You should have received a copy of the GNU General Public License along
 ! with this program; if not, write to the Free Software Foundation, Inc.,
 ! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 !
 !=====================================================================
 */

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC - ELASTIC coupling

/* ----------------------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_fluid_cmb_cuda,
              COMPUTE_COUPLING_FLUID_CMB_CUDA)(long* Mesh_pointer_f,
                                               int* FORWARD_OR_ADJOINT) {

  TRACE("compute_coupling_fluid_cmb_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_top_outer_core,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  // launches GPU kernel
  if( *FORWARD_OR_ADJOINT == 1 ){
    compute_coupling_fluid_CMB_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                                        mp->d_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        mp->nspec2D_top_outer_core);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_COUPLING();

    // adjoint simulations
    compute_coupling_fluid_CMB_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        mp->nspec2D_top_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_fluid_CMB_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_fluid_icb_cuda,
              COMPUTE_COUPLING_FLUID_ICB_CUDA)(long* Mesh_pointer_f,
                                               int* FORWARD_OR_ADJOINT) {

  TRACE("compute_coupling_fluid_icb_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_bottom_outer_core,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  // launches GPU kernel
  if( *FORWARD_OR_ADJOINT == 1 ){
    compute_coupling_fluid_ICB_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                                        mp->d_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        mp->nspec2D_bottom_outer_core);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_COUPLING();

    // adjoint simulations
    compute_coupling_fluid_ICB_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        mp->nspec2D_bottom_outer_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_fluid_ICB_kernel");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// ELASTIC - ACOUSTIC coupling

/* ----------------------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(compute_coupling_cmb_fluid_cuda,
              COMPUTE_COUPLING_CMB_FLUID_CUDA)(long* Mesh_pointer_f,
                                               int* FORWARD_OR_ADJOINT) {

  TRACE("compute_coupling_cmb_fluid_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_bottom_crust_mantle,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(5,5,1);

  // launches GPU kernel
  if( *FORWARD_OR_ADJOINT == 1 ){
    compute_coupling_CMB_fluid_kernel<<<grid,threads>>>(mp->d_displ_crust_mantle,
                                                        mp->d_accel_crust_mantle,
                                                        mp->d_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        mp->RHO_TOP_OC,
                                                        mp->minus_g_cmb,
                                                        mp->gravity,
                                                        mp->nspec2D_bottom_crust_mantle);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_COUPLING();

    //  adjoint simulations
    compute_coupling_CMB_fluid_kernel<<<grid,threads>>>(mp->d_b_displ_crust_mantle,
                                                        mp->d_b_accel_crust_mantle,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_crust_mantle,
                                                        mp->d_ibelm_bottom_crust_mantle,
                                                        mp->d_normal_top_outer_core,
                                                        mp->d_jacobian2D_top_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_top_outer_core,
                                                        mp->RHO_TOP_OC,
                                                        mp->minus_g_cmb,
                                                        mp->gravity,
                                                        mp->nspec2D_bottom_crust_mantle);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_CMB_fluid_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */
extern "C"
void FC_FUNC_(compute_coupling_icb_fluid_cuda,
              COMPUTE_COUPLING_ICB_FLUID_CUDA)(long* Mesh_pointer_f,
                                               int* FORWARD_OR_ADJOINT) {

  TRACE("compute_coupling_icb_fluid_cuda");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->nspec2D_top_inner_core,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLLX,NGLLX,1);

  // launches GPU kernel
  if( *FORWARD_OR_ADJOINT == 1 ){
    compute_coupling_ICB_fluid_kernel<<<grid,threads>>>(mp->d_displ_inner_core,
                                                        mp->d_accel_inner_core,
                                                        mp->d_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        mp->RHO_BOTTOM_OC,
                                                        mp->minus_g_icb,
                                                        mp->gravity,
                                                        mp->nspec2D_top_inner_core);
  }else if( *FORWARD_OR_ADJOINT == 3 ){
    // debug
    DEBUG_BACKWARD_COUPLING();

    //  adjoint simulations
    compute_coupling_ICB_fluid_kernel<<<grid,threads>>>(mp->d_b_displ_inner_core,
                                                        mp->d_b_accel_inner_core,
                                                        mp->d_b_accel_outer_core,
                                                        mp->d_ibool_inner_core,
                                                        mp->d_ibelm_top_inner_core,
                                                        mp->d_normal_bottom_outer_core,
                                                        mp->d_jacobian2D_bottom_outer_core,
                                                        mp->d_wgllwgll_xy,
                                                        mp->d_ibool_outer_core,
                                                        mp->d_ibelm_bottom_outer_core,
                                                        mp->RHO_BOTTOM_OC,
                                                        mp->minus_g_icb,
                                                        mp->gravity,
                                                        mp->nspec2D_top_inner_core);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);
  exit_on_cuda_error("compute_coupling_ICB_fluid_cuda");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

/* OCEANS load coupled on free surface */

/* ----------------------------------------------------------------------------------------------- */


/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(compute_coupling_ocean_cuda,
              COMPUTE_COUPLING_OCEAN_CUDA)(long* Mesh_pointer_f,
                                           int* FORWARD_OR_ADJOINT) {

  TRACE("compute_coupling_ocean_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); //get mesh pointer out of fortran integer container

  int blocksize = BLOCKSIZE_TRANSFER;
  int size_padded = ((int)ceil(((double)mp->npoin_oceans)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // uses corrected mass matrices
  if( *FORWARD_OR_ADJOINT == 1 ){
    compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_accel_crust_mantle,
                                                         mp->d_rmassx_crust_mantle,
                                                         mp->d_rmassy_crust_mantle,
                                                         mp->d_rmassz_crust_mantle,
                                                         mp->d_rmass_ocean_load,
                                                         mp->npoin_oceans,
                                                         mp->d_ibool_ocean_load,
                                                         mp->d_normal_ocean_load);
  }else if( *FORWARD_OR_ADJOINT == 3){
    // debug
    DEBUG_BACKWARD_COUPLING();

    // for backward/reconstructed potentials
    compute_coupling_ocean_cuda_kernel<<<grid,threads>>>(mp->d_b_accel_crust_mantle,
                                                         mp->d_b_rmassx_crust_mantle,
                                                         mp->d_b_rmassy_crust_mantle,
                                                         mp->d_b_rmassz_crust_mantle,
                                                         mp->d_rmass_ocean_load,
                                                         mp->npoin_oceans,
                                                         mp->d_ibool_ocean_load,
                                                         mp->d_normal_ocean_load);
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("compute_coupling_ocean_cuda");
#endif
}

