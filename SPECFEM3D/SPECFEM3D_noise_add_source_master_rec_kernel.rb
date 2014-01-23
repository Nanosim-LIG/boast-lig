module BOAST

  def BOAST::noise_add_source_master_rec_kernel(ref = true, n_gll3 = 125)
    push_env( :array_start => 0 )
    kernel = CKernel::new

    function_name = "noise_add_source_master_rec_kernel"

    ibool =              Int( "ibool",              :dir => :in, :dim => [Dim()] )
    ispec_selected_rec = Int( "ispec_selected_rec", :dir => :in, :dim => [Dim()] )
    irec_master_noise  = Int("irec_master_noise",   :dir => :in )
    accel =             Real("accel",               :dir => :out,:dim => [Dim()] )
    noise_sourcearray = Real("noise_sourcearray",   :dir => :in, :dim => [Dim()] )
    it =                 Int("it",                  :dir => :in )

    ngll3 = Int("NGLL3", :const => n_gll3)

    p = Procedure(function_name, [ibool, ispec_selected_rec, irec_master_noise, accel, noise_sourcearray, it], [ngll3])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL or get_lang == CUDA) then
      if (get_lang == CL) then
        if get_default_real_size == 8 then
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
        end
        load "atomicAdd_f.rb"
      end
      decl p
        decl tx  = Int("tx")
        decl ispec = Int("ispec")
        decl iglob = Int("iglob")

        print tx === get_local_id(0)
        print ispec === ispec_selected_rec[irec_master_noise] - 1
        print iglob === ibool[tx + ngll3*ispec] - 1

        (0..2).each { |i|
          print atomicAdd(accel + iglob*3 + i, noise_sourcearray[tx*3 + ngll3*3*it + i])
        }
        
      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end

end
