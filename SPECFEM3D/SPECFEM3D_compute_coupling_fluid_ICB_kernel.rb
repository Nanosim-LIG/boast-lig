require './SPECFEM3D_compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_fluid_ICB_kernel(ref=true)
    return BOAST::compute_coupling_kernel(ref, :fluid_ICB)
  end
end
