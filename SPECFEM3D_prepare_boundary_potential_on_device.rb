module ConvolutionGenerator
  def ConvolutionGenerator::prepare_boundary_potential_on_device
    old_array_start = $array_start
    $array_start = 0
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "prepare_boundary_potential_on_device"
    num_interfaces = Variable::new("num_interfaces",Int,{:direction => :in})
    max_nibool_interfaces = Variable::new("max_nibool_interfaces",Int,{:direction => :in})
    d_potential_dot_dot_acoustic = Variable::new("d_potential_dot_dot_acoustic",Real,{:direction => :in, :dimension => [ Dimension::new() ]})
    d_send_potential_dot_dot_buffer = Variable::new("d_send_potential_dot_dot_buffer",Real,{:direction => :out, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces) ]})
    d_nibool_interfaces = Variable::new("d_nibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces) ]})
    d_ibool_interfaces = Variable::new("d_ibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces) ]})
    if kernel.lang == ConvolutionGenerator::CL and ConvolutionGenerator::get_default_real_size == 8 then
      $output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    end
    p = Procedure::new(function_name, [d_potential_dot_dot_acoustic,d_send_potential_dot_dot_buffer,num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces])
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
          (d_send_potential_dot_dot_buffer[iloc] === d_potential_dot_dot_acoustic[iglob]).print
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
