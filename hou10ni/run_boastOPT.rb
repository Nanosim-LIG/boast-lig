require 'BOAST'
require './KRef.rb'
require './KBoastOPT.rb'
require 'narray_ffi'
require 'yaml'
require 'pp'
require 'csv'
require 'optparse'

include BOAST

class Experiment

	def self.run()

  # set_default_int_size(8)
  set_default_real_size(4)

  set_lang(FORTRAN)
  k={}
  stats={}
	repeat = 50

  nested = 1
	main = 1
	threads = 2 
	flags="-O2"  
	epsilon = 10E-10

  	#k_boast_params = {:kernel => :boast, openmp: true, :optim_nested => nested , optim_main: main, :omp_num_threads => threads,  :LDFLAGS => "-lgfortran -L/usr/lib/ -lblas", :FCFLAGS => "-fimplicit-none #{flags} -fexternal-blas"}  
  	k_boast_params = {:kernel => :boast, openmp: true, :optim_nested => nested , optim_main: main,  :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none #{flags}", :omp_num_threads => threads}  

  	k = KBoastOPT::new(k_boast_params)
  	k.generate
  	k.kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
  	stats[k_boast_params]={:time => []}
		puts k.kernel

	 # input and output parameters to check the kernel
   inputs = k.kernel.load_ref_inputs()
   outputs = k.kernel.load_ref_outputs()

  	inputs.each_key { |key|
			repeat.times{|i|
 				stats[k_boast_params][:time][i] = k.kernel.run(*(inputs[key]))[:duration]
			}
			puts k.kernel.compare_ref(outputs[key], inputs[key], epsilon).inspect
		}

		min = stats[k_boast_params][:time].min
		p min
  	min	

	puts stats
	puts "Testing done\n"

  return stats,k
 end
end


stats,k = Experiment.run

