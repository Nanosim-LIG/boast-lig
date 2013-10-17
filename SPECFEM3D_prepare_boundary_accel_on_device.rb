module ConvolutionGenerator
  def ConvolutionGenerator::prepare_boundary_accel_on_device
    old_array_start = $array_start
    $array_start = 0
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "prepare_boundary_accel_on_device"
    num_interfaces = Variable::new("num_interfaces",Int,{:direction => :in})
    max_nibool_interfaces = Variable::new("max_nibool_interfaces",Int,{:direction => :in})
    d_accel = Variable::new("d_accel",Real,{:direction => :in, :dimension => [ Dimension::new() ]})
    d_send_accel_buffer = Variable::new("d_send_accel_buffer",Real,{:direction => :out, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces*3) ]})
    d_nibool_interfaces = Variable::new("d_nibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces) ]})
    d_ibool_interfaces = Variable::new("d_ibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces) ]})
    if kernel.lang == ConvolutionGenerator::CL and ConvolutionGenerator::get_default_real_size == 8 then
      $output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    end
    p = Procedure::new(function_name, [d_accel,d_send_accel_buffer,num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces])
    if(ConvolutionGenerator::get_lang == ConvolutionGenerator::CUDA) then
      $output.print File::read("specfem3D/#{function_name}.cu")
    elsif(ConvolutionGenerator::get_lang == ConvolutionGenerator::CL) then
      p.decl
      id = Variable::new("id", Int)
      iglob = Variable::new("iglob", Int)
      iloc = Variable::new("iloc", Int)
      iinterface = Variable::new("iinterface", Int)
      (id === FuncCall::new("get_global_id",0)+FuncCall::new("get_global_size",0)*FuncCall::new("get_global_id",1)).print
      id.decl
      iglob.decl
      iloc.decl
      iinterface.decl
      f = For::new(iinterface, 0, num_interfaces-1) {
        cond = If::new(id<d_nibool_interfaces[iinterface]) {
          (iloc === id + max_nibool_interfaces*iinterface).print
          (iglob === d_ibool_interfaces[iloc] - 1).print
          (0..2).each { |i|
            (d_send_accel_buffer[iloc*3 + i] === d_accel[iglob*3 + i]).print
          }
        }
        cond.print
      }
      f.print
      p.close
    else
      raise "Unsupported language!"
    end
    kernel.procedure = p
    $array_start = old_array_start
    return kernel
  end
end
