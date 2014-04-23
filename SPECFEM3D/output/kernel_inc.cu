#include "assemble_boundary_accel_on_device.cu"
#include "assemble_boundary_potential_on_device.cu"
#include "prepare_boundary_potential_on_device.cu"
#include "prepare_boundary_accel_on_device.cu"
#include "get_maximum_scalar_kernel.cu"
#include "get_maximum_vector_kernel.cu"
#include "compute_add_sources_adjoint_kernel.cu"
#include "compute_add_sources_kernel.cu"
#include "compute_coupling_fluid_CMB_kernel.cu"
#include "compute_coupling_fluid_ICB_kernel.cu"
#include "compute_coupling_CMB_fluid_kernel.cu"
#include "compute_coupling_ICB_fluid_kernel.cu"
#include "compute_coupling_ocean_kernel.cu"
#include "write_seismograms_transfer_from_device_kernel.cu"
#include "write_seismograms_transfer_scalar_from_device_kernel.cu"
#include "noise_transfer_surface_to_host_kernel.cu"
#include "noise_add_source_master_rec_kernel.cu"
#include "noise_add_surface_movie_kernel.cu"
#include "compute_stacey_acoustic_kernel.cu"
#include "compute_stacey_acoustic_backward_kernel.cu"
#include "compute_stacey_elastic_kernel.cu"
#include "compute_stacey_elastic_backward_kernel.cu"
#include "update_disp_veloc_kernel.cu"
#include "update_potential_kernel.cu"
#include "update_accel_elastic_kernel.cu"
#include "update_veloc_elastic_kernel.cu"
#include "update_accel_acoustic_kernel.cu"
#include "update_veloc_acoustic_kernel.cu"
#include "outer_core_impl_kernel.cu"
#include "inner_core_impl_kernel.cu"
#include "compute_rho_kernel.cu"
#include "compute_iso_kernel.cu"
#include "compute_ani_kernel.cu"
#include "compute_hess_kernel.cu"
#include "compute_acoustic_kernel.cu"
#include "compute_strength_noise_kernel.cu"
#include "crust_mantle_impl_kernel.cu"