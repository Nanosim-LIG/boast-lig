Dir.chdir("..")
require './BOAST.rb'
Dir.chdir("SPECFEM3D")
require 'rubygems'
require 'narray'
require './HEADER.rb'

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end

kernels = [
:assemble_boundary_accel_on_device,
:assemble_boundary_potential_on_device,
:prepare_boundary_potential_on_device,
:prepare_boundary_accel_on_device,
:get_maximum_scalar_kernel,
:get_maximum_vector_kernel,
:compute_add_sources_adjoint_kernel,
:compute_add_sources_kernel,
:compute_coupling_fluid_CMB_kernel,
:compute_coupling_fluid_ICB_kernel,
:compute_coupling_CMB_fluid_kernel,
:compute_coupling_ICB_fluid_kernel,
:compute_coupling_ocean_kernel,
:write_seismograms_transfer_from_device_kernel,
:write_seismograms_transfer_scalar_from_device_kernel,
:noise_transfer_surface_to_host_kernel,
:noise_add_source_master_rec_kernel,
:noise_add_surface_movie_kernel,
:compute_stacey_acoustic_kernel,
:compute_stacey_acoustic_backward_kernel,
:compute_stacey_elastic_kernel,
:compute_stacey_elastic_backward_kernel,
:update_disp_veloc_kernel,
:update_potential_kernel,
:update_accel_elastic_kernel,
:update_veloc_elastic_kernel,
:update_accel_acoustic_kernel,
:update_veloc_acoustic_kernel,
:outer_core_impl_kernel,
:inner_core_impl_kernel,
:compute_rho_kernel,
:compute_iso_kernel,
:compute_ani_kernel,
:compute_hess_kernel,
:compute_acoustic_kernel,
:compute_strength_noise_kernel,
:crust_mantle_impl_kernel
]

langs = [ :CUDA, :CL]
BOAST::set_default_real_size(4)
BOAST::set_replace_constants(false)
class Float
  alias_method :old_to_s, :to_s
  def to_s
    s = ""+old_to_s
    s += "f" if BOAST::get_default_real_size == 4
    return s
  end
end

kernels.each { |kern|
  require "./SPECFEM3D_#{kern.to_s}.rb"
  puts kern.to_s
  langs.each { |lang|
    puts lang.to_s
    BOAST::set_lang( BOAST::const_get(lang))
    puts "REF" if lang == :CUDA
    k = BOAST::method(kern).call
    k.print
    if lang == :CUDA then
      puts "Generated"
      k = BOAST::method(kern).call(false)
      k.print
      filename = "#{kern}_cuda.c"
    elsif lang == :CL
      filename = "#{kern}_cl.c"
    end
    f = File::new("./specfem3D_output/"+filename, "w+")
    if lang == :CUDA then
      f.puts k
      k.build( :LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64", :NVCCFLAGS => "-arch sm_20 -O2" )
    elsif lang == :CL then
      s = k.to_s
      res = "const char * #{kern}_program = \"\\\n"
      s.each_line { |line|
        res += line.sub("\n","\\n\\\n")
      }
      res += "\";\n"
      f.print res
      k.build
    end
    f.close
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
#BOAST::set_lang( BOAST::FORTRAN )
#puts "FORTRAN"
#k = BOAST::vector_add
#puts k.print
#k.run(n,a,b,c_ref)
#BOAST::set_lang( BOAST::C )
#puts "C"
#c.random
#k = BOAST::vector_add
#puts k.print
#k.run(n,a,b,c)
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
#BOAST::set_lang( BOAST::CL )
#puts "CL"
#c.random
#k = BOAST::vector_add
#puts k.print
#k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
#BOAST::set_lang( BOAST::CUDA )
#puts "CUDA"
#c.random
#k = BOAST::vector_add
#puts k.print
#k.build(:LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64")
#k.run(n, a, b, c, :block_number => [rndup(n,32)/32, 1,1], :block_size => [32,1,1] )
#diff = (c_ref - c).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}
