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

	svg_opt_space = OptimizationSpace::new( optim_nested: 1..10, 
                                    optim_main: 1..10,
                                    OFLAGS: ["-O2", "-O3"], 
																		:omp_num_threads => 1..6
                                    )
	opt_space = OptimizationSpace::new( nested: [1], 
                                    main: [1],
                                    OFLAGS: ["-O2"],
																		:threads => 1..4
                                    )
  set_lang(FORTRAN)
  kernels={}
  stats={}
	repeat = 50

  #.. REF KERNEL ..#

  k_ref_params = {:kernel => :ref, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none -O2"}
  kernels[k_ref_params] = KRef::new(k_ref_params)
  kernels[k_ref_params].generate
  kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS], :FCFLAGS => k_ref_params[:FCFLAGS] ) 
  stats[k_ref_params]={:time => []}


	# input and output parameters to check the kernel
	# found in fluid_inner_elt_ref
  inputs = kernels[k_ref_params].kernel.load_ref_inputs()
  outputs = kernels[k_ref_params].kernel.load_ref_outputs()


	# As many keys as directories in fluid_inner_elt_ref/
  inputs.each_key { |key|
		repeat.times{|i|
 			stats[k_ref_params][:time][i] = kernels[k_ref_params].kernel.run(*(inputs[key]))[:duration]
		}
  }
	p k_ref_params
	p stats[k_ref_params][:time].min



  #.. BOAST KERNEL ..#
  
	optimizer = BruteForceOptimizer::new(opt_space, :randomize => true)
	puts optimizer.optimize { |opt|
		p opt
  	#k_boast_params = {:kernel => :boast, :OPENMP => true, :optim_nested => opt[:nested] , optim_main: opt[:main], :omp_num_threads => opt[:threads],  :LDFLAGS => "-lgfortran -L/usr/lib/ -lblas", :FCFLAGS => "-fimplicit-none #{opt[:OFLAGS]} -fexternal-blas"}  
  	k_boast_params = {:kernel => :boast, :OPENMP => true, :optim_nested => opt[:nested] , optim_main: opt[:main], :omp_num_threads => opt[:threads],  :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none #{opt[:OFLAGS]}"}  

  	kernels[k_boast_params] = KBoastOPT::new(k_boast_params)
  	kernels[k_boast_params].generate
  	kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
  	stats[k_boast_params]={:time => []}


  	inputs.each_key { |key|
			repeat.times{|i|
 				stats[k_boast_params][:time][i] = kernels[k_boast_params].kernel.run(*(inputs[key]))[:duration]
			}
			puts kernels[k_boast_params].kernel.compare_ref(outputs[key], inputs[key]).inspect
		}

		min = stats[k_boast_params][:time].min
		p min
  	min	
	}

	puts "Testing done\n"

  return stats,kernels
 end
end


options = {}

opt_parser = OptionParser.new { |opts|
  opts.banner = "Usage: run.rb --[options]=[value]"

  opts.on("-dVAL", "--data=VAL", "Specify the path where to store the data in yaml format") { |n|
    options[:data_path] = n
  }

  opts.on("-kVAL", "--kernel=VAL", "Specify the path to store the kernel sources") { |n|
    options[:kernel_path] = n
  }
  opts.on("-h", "--help", "Prints this help") {
    puts opts
    exit
  }
}.parse!



stats,kernels = Experiment.run


if options[:data_path] then
  File::open( options[:data_path], "w") { |f|
    f.print YAML::dump(stats)
  }
end


if options[:kernel_path] then
  File::open( options[:kernel_path], "w") { |f|
    kernels.each { |key, value|
      f.puts value.kernel if key[:kernel] == :boast
    }
  }
end
