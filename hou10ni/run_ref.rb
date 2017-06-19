require 'BOAST'
require './KRef.rb'
require 'narray_ffi'

include BOAST

class Experiment

 def self.run()
  k_ref_params = {:kernel => :ref, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -O2"}  

  set_lang(FORTRAN)
  kernels={}

  # Creating ref kernel
  kernels[k_ref_params] = KRef::new(k_ref_params)
  kernels[k_ref_params].generate
	puts "*********** Building kernel ***********\n"
  puts kernels[k_ref_params].kernel
  kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS], :FCFLAGS => k_ref_params[:FCFLAGS] )
 
  #puts "Parameters = ", kernels[k_ref_params].kernel.procedure.parameters

	# Parameters for loop iteration 1 and MPI Process 0
  puts "** inputs from the application ="
  puts "{ ./fluid_inner_elt_ref/sub_time_loop_1 =>[210025, 0, 1, NArray.int(2577732): "
  puts "[ 5, 1, 4, 538597, 538600, 319233, ...], NArray.int(214812):"
  puts "[ 1, 81, 161, 209, 289, 369, ...],  NArray.float(879732):"
  puts "[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ... ], NArray.float(879732):" 
  puts "[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ... ], NArray.float(879732):" 
  puts "[ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ... ], NArray.float(16952560):"
  puts "[1.877937, 1.1007047E-02, 1.0670856E-02, 1.9009672E-02, 1.1590085E-02, ...], 16952560, 879732, 12, 214811, 214812]} " 

	puts "\n*********** Testing ref kernel with fluid_inner_elt_ref/sub_time_loop_1/ inputs parameters ***********\n"
  inputs = kernels[k_ref_params].kernel.load_ref_inputs()
  outputs = kernels[k_ref_params].kernel.load_ref_outputs()
	puts "** inputs =", inputs
  puts "** outputs=",outputs


  inputs.each_key { |key|
		puts key
 		puts kernels[k_ref_params].kernel.run(*(inputs[key])).inspect 
		puts kernels[k_ref_params].kernel.compare_ref(outputs[key], inputs[key]).inspect
	}

	puts "Testing done\n"
  return kernels
 end
end

Experiment.run
