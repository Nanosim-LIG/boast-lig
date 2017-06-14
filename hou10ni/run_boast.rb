require 'BOAST'
require './KBoast.rb'
require 'narray_ffi'

include BOAST

class Experiment

 def self.run()

 # set_default_int_size(8)
  set_default_real_size(4)
  k_boast_params = {:kernel => :boast, :preprocessor => false, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none"}  

  set_lang(FORTRAN)
  kernels = {}
  stats_svg = {}
	stats = []
  repeat = 3

  # Creating ref kernel
  kernels[k_boast_params] = KBoast::new(k_boast_params)
  kernels[k_boast_params].generate
	puts "*********** Building kernel ***********\n"
  puts kernels[k_boast_params].kernel
  kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
 
	puts "\n*********** Testing boast kernel with fluid_inner_elt_boast/sub_time_loop_1/ inputs parameters ***********\n"
  inputs = kernels[k_boast_params].kernel.load_ref_inputs()
  outputs = kernels[k_boast_params].kernel.load_ref_outputs()

	repeat.times {
  	inputs.each_key { |key|
			puts key
 			stats.push kernels[k_boast_params].kernel.run(*(inputs[key])) 
 			puts kernels[k_boast_params].kernel.run(*(inputs[key])).inspect 
			puts kernels[k_boast_params].kernel.compare_ref(outputs[key], inputs[key]).inspect
		}
	}
  stats.sort_by! { |a| a[:duration] }
	stats = stats.first



  puts "#{kernels[k_boast_params].kernel.procedure.name}: #{stats[:duration]} s ->  #{stats[:duration]*1.0e3} ms"
  return kernels
 end
end

Experiment.run
