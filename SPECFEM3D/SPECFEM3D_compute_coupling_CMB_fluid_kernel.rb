require './SPECFEM3D_compute_coupling_fluid_CMB_kernel.rb'
module BOAST
  def BOAST::compute_coupling_CMB_fluid_kernel(ref = true)
    BOAST::compute_coupling_kernel(ref, :CMB_fluid)
  end
end