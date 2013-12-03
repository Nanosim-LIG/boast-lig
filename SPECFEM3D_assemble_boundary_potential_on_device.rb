module ConvolutionGenerator
  def ConvolutionGenerator::assemble_boundary_potential_on_device
    old_array_start = @@array_start
    @@array_start = 0
    kernel = CKernel::new
    function_name = "assemble_boundary_potential_on_device"
    num_interfaces = Variable::new("num_interfaces",Int,{:direction => :in})
    max_nibool_interfaces = Variable::new("max_nibool_interfaces",Int,{:direction => :in})
    d_potential_dot_dot_acoustic = Variable::new("d_potential_dot_dot_acoustic",Real,{:direction => :out, :dimension => [ Dimension::new() ]})
    d_send_potential_dot_dot_buffer = Variable::new("d_send_potential_dot_dot_buffer",Real,{:direction => :in, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces) ]})
    d_nibool_interfaces = Variable::new("d_nibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces) ]})
    d_ibool_interfaces = Variable::new("d_ibool_interfaces",Int,{:direction => :in, :dimension => [ Dimension::new(num_interfaces*max_nibool_interfaces) ]})
    if kernel.lang == ConvolutionGenerator::CL and ConvolutionGenerator::get_default_real_size == 8 then
      @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      @@output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
    end
    p = Procedure::new(function_name, [d_potential_dot_dot_acoustic,d_send_potential_dot_dot_buffer,num_interfaces,max_nibool_interfaces,d_nibool_interfaces,d_ibool_interfaces])
    if(ConvolutionGenerator::get_lang == ConvolutionGenerator::CUDA) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(ConvolutionGenerator::get_lang == ConvolutionGenerator::CL) then
      type_f = Real::new.decl
      if ConvolutionGenerator::get_default_real_size == 8 then
        type_i = "unsigned long int"
        cmpx_name = "atom_cmpxchg"
      else
        type_i = "unsigned int"
        cmpx_name = "atomic_cmpxchg"
      end
      @@output.print <<EOF
static inline void atomicAdd_f(volatile __global float *source, const float val) {
  union {
    #{type_i} iVal;
    #{type_f} fVal;
  } res, orig;
  do {
    orig.fVal = *source;
    res.fVal = orig.fVal + val;
  } while (#{cmpx_name}((volatile __global #{type_i} *)source, orig.iVal, res.iVal) != orig.iVal);
}
EOF
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
          (FuncCall::new("atomicAdd_f",d_potential_dot_dot_acoustic + iglob, d_send_potential_dot_dot_buffer[iloc])).print
        }
        cond.print
      }
      f.print
      p.close
    else
      raise "Unsupported language!"
    end
    kernel.procedure = p
    @@array_start = old_array_start
    return kernel
  end
end

