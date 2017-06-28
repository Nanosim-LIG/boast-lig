require 'BOAST'
require './KBoast.rb'
require 'narray_ffi'

include BOAST

optim_nested = 3
optim_main = 1
num_threads = 1

#set_default_int_size(8)
set_default_real_size(4)

# Parameters - optim = number of unroll of the nested loop per iteration
k_boast_params = {:kernel => :boast, omp_num_threads: num_threads, optim_nested: optim_nested, optim_main: optim_main, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fbounds-check -fimplicit-none -O2"}  

set_lang(FORTRAN)
stats = []
repeat = 5

# Creating boast kernel
k = KBoast::new(k_boast_params)
k.generate
puts k.kernel
k.kernel.build(:OPENMP => true, :LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
 
inputs = k.kernel.load_ref_inputs()
outputs = k.kernel.load_ref_outputs()
puts "** inputs =", inputs
puts "** outputs=",outputs


inputs.each_key { |key|
	repeat.times {
 		stats.push k.kernel.run(*(inputs[key])) 
	}
	puts k.kernel.compare_ref(outputs[key], inputs[key]).inspect
}

stats.sort_by! { |a| a[:duration] }
p stats
stats = stats.first

puts "#{k.kernel.procedure.name}: optim_nested = #{k_boast_params[:optim_nested]}  #{stats[:duration]} s ->  #{stats[:duration]*1.0e3} ms"

