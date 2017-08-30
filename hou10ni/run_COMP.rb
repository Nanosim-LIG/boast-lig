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


  #generer la sortie avec la ref et comparer avec ce qu'a boast
  # utiliser des numeros pour les inputs

  set_lang(FORTRAN)
  kernels={}
  stats={}
	repeat = 1

  #.. REF KERNEL ..#

  k_ref_params = {:kernel => :ref, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none -O2"}
  kernels[k_ref_params] = KRef::new(k_ref_params)
  kernels[k_ref_params].generate
  kernels[k_ref_params].kernel.build(:LDFLAGS => k_ref_params[:LDFLAGS], :FCFLAGS => k_ref_params[:FCFLAGS] ) 
  stats[k_ref_params]={:time => []}

  inputs = kernels[k_ref_params].kernel.load_ref_inputs()
  outputs = kernels[k_ref_params].kernel.load_ref_outputs()

  inputs.each_key { |key|
 		stats[k_ref_params].push kernels[k_ref_params].kernel.run(*(inputs[key]))
  }

  #.. BOAST KERNEL ..#
  
  	k_boast_params = {:kernel => :boast, :OPENMP => true, :optim_nested => 1 , optim_main: 1, :omp_num_threads => 1,  :LDFLAGS => "-lgfortran -L/usr/lib/ -lblas", :FCFLAGS => "-fimplicit-none -O2 -fexternal-blas"}  

  	kernels[k_boast_params] = KBoast::new(k_boast_params)
  	kernels[k_boast_params].generate
  	kernels[k_boast_params].kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
  	stats[k_boast_params]={:time => []}

  	inputs.each_key { |key|
			stats[k_boast_params].push kernels[k_boast_params].kernel.run(*(inputs[key]))
		}


	puts "Testing done\n"

	## TODO: compare ref and boast result (P_new_ref and P_new_boast )
	#diff = (c_ref - c).abs
	#diff.each { |elem|
  	#raise "Warning: residue too big: #{elem}" if elem > epsilon
	#}

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
