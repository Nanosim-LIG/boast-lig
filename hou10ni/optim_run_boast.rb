require 'BOAST'
require './KBoast.rb'
require 'narray_ffi'

include BOAST

class Experiment

 def self.run()

 
  #set_default_int_size(8)
  set_default_real_size(4)

	# Parameters - optim = number of unroll of the nested loop per iteration
  opts1 = {:kernel => :boast, :optim_nested => 1, :optim_main => 1, :preprocessor => false}  
  opts2 = {:kernel => :boast, :optim_nested => 5, :optim_main => 25, :preprocessor => false}  
  opts3 = {:kernel => :boast, :optim_nested => 5, :optim_main => 50, :preprocessor => false}  

  set_lang(FORTRAN)
  kernels = {}
	stats = {}
  repeat = 10

  # Creating boast kernels depending of the optimization #

  kernels[opts1] = KBoast::new(opts1)
  kernels[opts1].generate
  #puts kernels[opts1].kernel
  kernels[opts1].kernel.build(:LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none -O3" )

  kernels[opts2] = KBoast::new(opts2)
  kernels[opts2].generate
  #puts kernels[opts2].kernel
  kernels[opts2].kernel.build(:LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none -O3" )

  kernels[opts3] = KBoast::new(opts3)
  kernels[opts3].generate
  #puts kernels[opts3].kernel
  kernels[opts3].kernel.build(:LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none -O3" )


  stats[opts1] = {:time => []}
  stats[opts2] = {:time => []}
  stats[opts3] = {:time => []}


	# Test boast kernels with fluid_inner_elt_boast/sub_time_loop_1/ inputs parameters 
  inputs = kernels[opts1].kernel.load_ref_inputs()
  outputs = kernels[opts1].kernel.load_ref_outputs()


	repeat.times { |i|
  	inputs.each_key { |key|
 			stats[opts1][:time][i] =  kernels[opts1].kernel.run(*(inputs[key]))[:duration] 
 			puts kernels[opts1].kernel.run(*(inputs[key])).inspect 
			puts kernels[opts1].kernel.compare_ref(outputs[key], inputs[key]).inspect
 
			stats[opts2][:time][i] =  kernels[opts2].kernel.run(*(inputs[key]))[:duration] 
 			puts kernels[opts2].kernel.run(*(inputs[key])).inspect 
			puts kernels[opts2].kernel.compare_ref(outputs[key], inputs[key]).inspect

 			stats[opts3][:time][i] =  kernels[opts3].kernel.run(*(inputs[key]))[:duration] 
 			puts kernels[opts3].kernel.run(*(inputs[key])).inspect 
			puts kernels[opts3].kernel.compare_ref(outputs[key], inputs[key]).inspect
		}
	}

  return stats, kernels
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
      f.puts value.kernel if key[:optim_nested] == 1
    }
  }
end
