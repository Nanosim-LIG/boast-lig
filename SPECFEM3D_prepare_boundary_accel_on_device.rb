module BOAST
  def BOAST::prepare_boundary_accel_on_device
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "prepare_boundary_accel_on_device"
    num_interfaces =        Int( "num_interfaces",        :dir => :in)
    max_nibool_interfaces = Int( "max_nibool_interfaces", :dir => :in)
    d_accel =               Real("d_accel",               :dir => :in,  :dim => [ Dim() ] )
    d_send_accel_buffer =   Real("d_send_accel_buffer",   :dir => :out, :dim => [ Dim(num_interfaces*max_nibool_interfaces*3) ] )
    d_nibool_interfaces =   Int( "d_nibool_interfaces",   :dir => :in,  :dim => [ Dim(num_interfaces) ] )
    d_ibool_interfaces =    Int( "d_ibool_interfaces",    :dir => :in,  :dim => [ Dim(num_interfaces*max_nibool_interfaces) ] )
    p = Procedure::new(function_name, [d_accel,d_send_accel_buffer,num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces])
    if(get_lang == CUDA) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CL) then
      if get_default_real_size == 8 then
        @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      end
      decl p
      id =         Int("id")
      iglob =      Int("iglob")
      iloc =       Int("iloc")
      iinterface = Int("iinterface")
      decl id
      decl iglob
      decl iloc
      decl iinterface
      print id === get_global_id(0)+get_global_size(0)*get_global_id(1)
      print For(iinterface, 0, num_interfaces-1) {
        print If(id<d_nibool_interfaces[iinterface]) {
          print iloc === id + max_nibool_interfaces*iinterface
          print iglob === d_ibool_interfaces[iloc] - 1
          (0..2).each { |i|
            print d_send_accel_buffer[iloc*3 + i] === d_accel[iglob*3 + i]
          }
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
