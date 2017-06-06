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
  kernels={}

  # Creating ref kernel
  kernels[k_boast_params] = KBoast::new(k_boast_params)
  kernels[k_boast_params].generate
	puts "*********** Building kernel ***********\n"
  puts kernels[k_boast_params].kernel
  kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
 
	puts "\n*********** Testing boast kernel with fluid_inner_elt_boast/sub_time_loop_1/ inputs parameters ***********\n"
  inputs = kernels[k_boast_params].kernel.load_ref_inputs()
  outputs = kernels[k_boast_params].kernel.load_ref_outputs()


  inputs.each_key { |key|
		puts key
 		puts kernels[k_boast_params].kernel.run(*(inputs[key])).inspect 
		puts kernels[k_boast_params].kernel.compare_ref(outputs[key], inputs[key]).inspect
	}

	puts "Testing done\n"
  return kernels
 end
end

Experiment.run
