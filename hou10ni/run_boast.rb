require 'BOAST'
require './KBoast.rb'
require 'narray_ffi'

include BOAST

optim_nested = 3
optim_main = 1
num_threads = 1

#set_default_int_size(8)
set_default_real_size(4)

#k_boast_params = {:kernel => :boast, :LDFLAGS => "-lgfortran -L/usr/lib/ -lblas", :FCFLAGS => "-fimplicit-none -fexternal-blas -O3"}  
k_boast_params = {:kernel => :boast, :LDFLAGS => "-lgfortran", :FCFLAGS => "-fimplicit-none -O3"}  

set_lang(FORTRAN)
stats = []
repeat = 1

# Creating boast kernel
k = KBoast::new(k_boast_params)
k.generate
puts k.kernel
k.kernel.build(:LDFLAGS => k_boast_params[:LDFLAGS], :FCFLAGS => k_boast_params[:FCFLAGS] )
 
inputs = k.kernel.load_ref_inputs()
#outputs = k.kernel.load_ref_outputs()
#puts "** inputs =", inputs
#puts "** outputs=",outputs



inputs.each_key { |key|
	k.kernel.run(*(inputs[key]))
	repeat.times {
 		stats.push k.kernel.run(*(inputs[key])) 
	}

	#puts k.kernel.compare_ref(outputs[key], inputs[key]).inspect
}

stats.sort_by! { |a| a[:duration] }
p stats
stats = stats.first

#puts "#{k.kernel.procedure.name}: #{stats[:duration]*1.0e3} ms    #{ stats[:reference_return][:cost] / (stats[:duration]*1.0e9)} GFlops"
puts "#{k.kernel.procedure.name}: #{stats[:duration]*1.0e3} ms    #{ 68914896 / (stats[:duration]*1.0e9)} GFlops"

