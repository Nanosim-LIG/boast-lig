require "./SPECFEM3D_inner_core_impl_kernel.rb"

module BOAST
  def crust_mantle_impl_kernel( ref = true, mesh_coloring = false, textures_fields = false, textures_constants = false, unroll_loops = true, n_gllx = 5, n_gll2 = 25, n_gll3 = 125, n_gll3_padded = 128, r_earth_km = 6371.0, n_sls = 3)
    return BOAST::impl_kernel(:crust_mantle, ref, mesh_coloring, textures_fields, textures_constants, unroll_loops, n_gllx, n_gll2, n_gll3, n_gll3_padded, n_sls, r_earth_km, coloring_min_nspec_inner_core, i_flag_in_fictitious_cube, n_sls)
  end
end
