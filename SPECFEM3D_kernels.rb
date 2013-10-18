require './BOAST.rb'
require 'rubygems'
require 'narray'

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

kernels = [
:prepare_boundary_potential_on_device,
:assemble_boundary_potential_on_device,
:prepare_boundary_accel_on_device,
:assemble_boundary_accel_on_device,
:get_maximum_scalar_kernel,
:get_maximum_vector_kernel,
:compute_add_sources_adjoint_kernel
]

langs = [ :CUDA, :CL]
ConvolutionGenerator::set_default_real_size(4)

kernels.each { |kern|
  require "./SPECFEM3D_#{kern.to_s}.rb"
  puts kern.to_s
  langs.each { |lang|
    puts lang.to_s
    ConvolutionGenerator::set_lang( ConvolutionGenerator::const_get(lang))
    k = ConvolutionGenerator::method(kern).call
    k.print
  }
}

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
