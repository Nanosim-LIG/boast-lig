require './BOAST.rb'
require 'rubygems'
require 'narray'
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
      $output.print File::read("specfem3D/assemble_MPI_scalar_cuda.cu")
    elsif(ConvolutionGenerator::get_lang == ConvolutionGenerator::CL) then
      p.decl
      id = Varaiable::new("id", Int)
      iglob = Varaiable::new("iglob", Int)
      iloc = Varaiable::new("iloc", Int)
      iinterface = Varaiable::new("iinterface", Int)
      id.decl
      iglob.decl
      iloc.decl
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

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

#
#n = 1024*1024
#a = NArray.float(n).random
#b = NArray.float(n).random
#c = NArray.float(n)
#c_ref = NArray.float(n)
#
#epsilon = 10e-15
#
#ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
#puts "FORTRAN"
#k = ConvolutionGenerator::vector_add
#puts k.print
#k.run(n,a,b,c_ref)
#ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
#puts "C"
#c.random
#k = ConvolutionGenerator::vector_add
#puts k.print
#k.run(n,a,b,c)
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
#ConvolutionGenerator::set_lang( ConvolutionGenerator::CL )
#puts "CL"
#c.random
#k = ConvolutionGenerator::vector_add
#puts k.print
#k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
#ConvolutionGenerator::set_lang( ConvolutionGenerator::CUDA )
#puts "CUDA"
#c.random
#k = ConvolutionGenerator::vector_add
#puts k.print
#k.build(:LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64")
#k.run(n, a, b, c, :block_number => [rndup(n,32)/32, 1,1], :block_size => [32,1,1] )
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
