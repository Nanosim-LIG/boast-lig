require './SPECFEM3D_compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_ICB_fluid_kernel(ref=true)
    return BOAST::compute_coupling_kernel(ref, :ICB_fluid)
  end
end
