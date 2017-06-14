require 'BOAST'
require './KRef.rb'
require './KBoast.rb'
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

  k_ref_params = {:kernel => :ref, :preprocessor => false, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check"}
  k_boast_params = {:kernel => :boast, :preprocessor => false, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none -O3"}  

  set_lang(FORTRAN)
  kernels={}
  stats={}

  # Creating ref kernel
  kernels[k_ref_params] = KRef::new(k_ref_params)
  kernels[k_ref_params].generate
	puts "* Building Ref kernel ...\n"
  #puts kernels[k_ref_params].kernel
  kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS], :FCFLAGS => k_ref_params[:FCFLAGS] ) 
  stats[k_ref_params]={:time => []}

  # Creating boast kernel
  kernels[k_boast_params] = KBoast::new(k_boast_params)
  kernels[k_boast_params].generate
	puts "* Building BOAST kernel ...\n"
  #puts kernels[k_boast_params].kernel
  kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
  stats[k_boast_params]={:time => []}

  # Run (inputs are the same for ref and boast kernels) 
	puts "\n* Testing boast and ref versions of the kernel with fluid_inner_elt_boast/sub_time_loop_1/ inputs parameters ...\n"
  inputs = kernels[k_boast_params].kernel.load_ref_inputs()
  outputs = kernels[k_boast_params].kernel.load_ref_outputs()

10.times{|i|
  inputs.each_key { |key|
		#puts "KEY=  #{inputs[key]}"
 		#puts kernels[k_ref_params].kernel.run(*(inputs[key])).inspect 
 		stats[k_ref_params][:time][i] = kernels[k_ref_params].kernel.run(*(inputs[key]))[:duration]
 		stats[k_boast_params][:time][i] = kernels[k_boast_params].kernel.run(*(inputs[key]))[:duration]
		puts kernels[k_boast_params].kernel.compare_ref(outputs[key], inputs[key]).inspect
	}
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
