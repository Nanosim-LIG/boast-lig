require './GenericConvolution.rb'
require './WaveletFilters.rb'
require './MagicFilter.rb'
openmp = true
def MagicFilter1d(filter, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(filter, :ld => true)
  conv_operation.optimize(optims) if optims

  p, subops = conv_operation.procedure

  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op|
    BOAST::pr op
    puts "chosen:"+ op.name
  }
  BOAST::pr p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
  return kernel
end

n1 = 32
n2 = 32
n3 = 32
optims = GenericOptimization::new(:repeat => 9, :unroll_range => [8,8], :mod_arr_test => true, :tt_arr_test => true, :dimensions => [n1,n2,n3], :vector_length => [2,4], :unrolled_dim_index_test => true, :openmp => openmp)
conv_filter = ConvolutionFilter::new('sfrf',SYM8_MF.reverse,8)
ksym8 = MagicFilter1d( conv_filter, optims )
#puts ksym8

input = ANArray.float(64,n1,n2,n3).random!
work1 = ANArray.float(64,n1,n2,n3)
work2 = ANArray.float(64,n1,n2,n3)
output = ANArray.float(64,n1,n2,n3)
output_ref = ANArray.float(64,n1,n2,n3)
epsilon = 10e-15

k_ref = magicfilter_ref
stats = k_ref.run(n1, n2*n3, input, output_ref)
stats = k_ref.run(n2, n1*n3, output_ref, work1)
stats = k_ref.run(n3, n2*n1, work1, output_ref)
k_ref.dump_source

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3

ksym8.build(:openmp => openmp)
stats0 = ksym8.run(3, 0, n, BC::PERIODIC, n, n, input, work1)
stats0 = ksym8.run(3, 0, n, BC::PERIODIC, n, n, input, work1)
stats1 = ksym8.run(3, 1, n, BC::PERIODIC, n, n, work1, work2)
stats2 = ksym8.run(3, 2, n, BC::PERIODIC, n, n, work2, output)
puts "#{ksym8.procedure.name} d0: #{stats0[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats0[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name} d1: #{stats1[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats1[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name} d2: #{stats2[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats2[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name}: #{(stats0[:duration]+stats1[:duration]+stats2[:duration])*1.0e3} ms #{3*32*n1*n2*n3 / ((stats0[:duration]+stats1[:duration]+stats2[:duration])*1.0e9)} GFlops"
ksym8.dump_source

diff = (output_ref - output).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
