require "./SPECFEM3D_update_disp_veloc_kernel.rb"

module BOAST

  def BOAST::update_veloc_acoustic_kernel(ref = true)
    BOAST::update_kernel(:veloc_acoustic, ref)
  end

end
