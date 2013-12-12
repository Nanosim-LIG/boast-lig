module BOAST
  def BOAST::assemble_boundary_potential_on_device(ref = true)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "assemble_boundary_potential_on_device"
    num_interfaces                  = Int("num_interfaces",                   :dir => :in)
    max_nibool_interfaces           = Int("max_nibool_interfaces",            :dir => :in)
    d_potential_dot_dot_acoustic    = Real("d_potential_dot_dot_acoustic",    :dir => :out,:dim => [ Dim() ])
    d_send_potential_dot_dot_buffer = Real("d_send_potential_dot_dot_buffer", :dir => :in, :dim => [ Dim(num_interfaces*max_nibool_interfaces) ])
    d_nibool_interfaces             = Int("d_nibool_interfaces",              :dir => :in, :dim => [ Dim(num_interfaces) ])
    d_ibool_interfaces              = Int("d_ibool_interfaces",               :dir => :in, :dim => [ Dim(num_interfaces*max_nibool_interfaces) ])
    p = Procedure(function_name, [d_potential_dot_dot_acoustic,d_send_potential_dot_dot_buffer,num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CUDA or get_lang == CL) then
      if(get_lang == CL) then
        if get_default_real_size == 8 then
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
          @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
        end
        load "./atomicAdd_f.rb"
      end
      decl p
      id         = Int("id")
      iglob      = Int("iglob")
      iloc       = Int("iloc")
      iinterface = Int("iinterface")
      decl id
      decl iglob
      decl iloc
      decl iinterface
      print id === get_global_id(0) + get_global_size(0)*get_global_id(1)
      print For(iinterface, 0, num_interfaces-1) {
        print If(id<d_nibool_interfaces[iinterface]) {
          print iloc === id + max_nibool_interfaces*iinterface
          print iglob === d_ibool_interfaces[iloc] - 1
          print atomicAdd(d_potential_dot_dot_acoustic + iglob, d_send_potential_dot_dot_buffer[iloc])
        }
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

