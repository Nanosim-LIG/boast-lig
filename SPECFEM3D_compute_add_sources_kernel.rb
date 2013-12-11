module BOAST
  def BOAST::compute_add_sources_kernel
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "compute_add_sources_kernel"
    accel =                  Real("accel",                            :dir => :out,:dim => [ Dim() ] )
    ibool =                  Int("ibool",                             :dir => :in, :dim => [ Dim() ] )
    sourcearrays =           Real("sourcearrays",                     :dir => :in, :dim => [ Dim() ] )
    stf_pre_compute =        Real("stf_pre_compute",      :size => 8, :dir => :in, :dim => [ Dim() ] )
    myrank =                 Int("myrank",                            :dir => :in)
    islice_selected_source = Int("islice_selected_source",            :dir => :in, :dim => [ Dim() ] )
    ispec_selected_source =  Int("ispec_selected_source",             :dir => :in, :dim => [ Dim() ] )
    nsources =               Int("nsources",                          :dir => :in)

    ndim =                   Int("NDIM",                  :const => 3)
    ngllx =                  Int("NGLLX",                 :const => 5)
    p = Procedure(function_name, [accel,ibool,sourcearrays,stf_pre_compute,myrank,islice_selected_source,ispec_selected_source,nsources], [ndim,ngllx])
    if(get_lang == CUDA) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL) then
      @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      if get_default_real_size == 8 then
        @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
      end
      load "./atomicAdd_f.rb"
      load "./INDEX4.rb"
      load "./INDEX5.rb"
      decl p
      ispec =   Int( "ispec")
      iglob =   Int( "iglob")
      stf =     Real("stf")
      isource = Int( "isource")
      i =       Int( "i")
      j =       Int( "j")
      k =       Int( "k")
      decl ispec
      decl iglob
      decl stf
      decl isource
      decl i
      decl j
      decl k
      print i === get_local_id(0)
      print j === get_local_id(1)
      print k === get_local_id(2)
      print isource === get_group_id(0) + get_num_groups(0)*get_group_id(1)

      print If(isource < nsources) {
        print If(myrank == islice_selected_source[isource]) {
          print ispec === ispec_selected_source[isource] - 1
          print stf === stf_pre_compute[isource]
          print iglob === ibool[INDEX4(ngllx,ngllx,ngllx,i,j,k,ispec)] - 1
          (0..2).each { |indx|
            print atomicAdd_f(accel+iglob*3+indx, sourcearrays[INDEX5(ndim,ngllx,ngllx,ngllx,indx,i,j,k,isource)]*stf)
          }
        }
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env(:array_start)
    kernel.procedure = p
    return kernel
  end
end
