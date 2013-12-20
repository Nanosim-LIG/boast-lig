module BOAST
  def BOAST::compute_coupling_fluid_CMB_kernel(ref = true)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_coupling_fluid_CMB_kernel"
    displ_crust_mantle        = Real("displ_crust_mantle",        :dir => :in,  :dim => [ Dim() ])
    accel_outer_core          = Real("accel_outer_core",          :dir => :out, :dim => [ Dim() ])
    ibool_crust_mantle        = Int( "ibool_crust_mantle",        :dir => :in,  :dim => [ Dim() ])
    ibelm_bottom_crust_mantle = Int( "ibelm_bottom_crust_mantle", :dir => :in,  :dim => [ Dim() ])
    normal_top_outer_core     = Real("normal_top_outer_core",     :dir => :in,  :dim => [ Dim() ])
    jacobian2D_top_outer_core = Real("jacobian2D_top_outer_core", :dir => :in,  :dim => [ Dim() ])
    wgllwgll_xy               = Real("wgllwgll_xy",               :dir => :in,  :dim => [ Dim() ])
    ibool_outer_core          = Int( "ibool_outer_core",          :dir => :in,  :dim => [ Dim() ])
    ibelm_top_outer_core      = Int( "ibelm_top_outer_core",      :dir => :in,  :dim => [ Dim() ])
    nspec2D_TOP_OC            = Int( "NSPEC2D_TOP_OC",            :dir => :in)

    ndim =               Int( "NDIM",              :const => 3)
    ngllx =              Int( "NGLLX",             :const => 5)
    p = Procedure(function_name, [displ_crust_mantle,accel_outer_core,ibool_crust_mantle,ibelm_bottom_crust_mantle,normal_top_outer_core,jacobian2D_top_outer_core,wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core,nspec2D_TOP_OC],[ndim, ngllx])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      if (get_lang == CL) then
        if get_default_real_size == 8 then
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
        end
        load "./atomicAdd_f.rb"
      end
      load "./INDEX2.rb"
      load "./INDEX3.rb"
      load "./INDEX4.rb"
      decl p
      decl i = Int("i"), j = Int("j"), k = Int("k")
      decl iface =      Int("iface")
      decl k_corresp =  Int("k_corresp")
      decl iglob_cm =   Int("iglob_cm")
      decl iglob_oc =   Int("iglob_oc")
      decl ispec =      Int("ispec")
      decl ispec_selected = Int("ispec_selected")
      decl *(displ = [Real("displ_x"), Real("displ_y"), Real("displ_z")])
      decl displ_n = Real("displ_n")
      decl *(n = [Real("nx"), Real("ny"), Real("nz")])
      decl weight = Real("weight")

      print i === get_local_id(0)
      print j === get_local_id(1)
      print iface === get_group_id(0) + get_num_groups(0)*get_group_id(1)
      print If( iface < nspec2D_TOP_OC ) {
        print ispec === ibelm_top_outer_core[iface] - 1
        print ispec_selected === ibelm_bottom_crust_mantle[iface] - 1
        print k === ngllx - 1
        print k_corresp === 0
        print iglob_cm === ibool_crust_mantle[INDEX4(ngllx,ngllx,ngllx,i,j,k_corresp,ispec_selected)] - 1
        displ.each_index { |indx| print displ[indx] === displ_crust_mantle[iglob_cm*3+indx] }
        n.each_index { |indx| print n[indx] === normal_top_outer_core[INDEX4(ndim,ngllx,ngllx,indx,i,j,iface)] }
        print displ_n === displ[0]*n[0] + displ[1]*n[1] + displ[2]*n[2]
        print weight === jacobian2D_top_outer_core[INDEX3(ngllx,ngllx,i,j,iface)]*wgllwgll_xy[INDEX2(ngllx,i,j)]
        print iglob_oc === ibool_outer_core[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
        print atomicAdd(accel_outer_core+iglob_oc, weight*displ_n)
      }
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end
end
