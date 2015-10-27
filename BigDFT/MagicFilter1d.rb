require './GenericConvolution.rb'
require './WaveletFilters.rb'
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
optims = GenericOptimization::new(:unroll_range => 6, :mod_arr_test => true, :tt_arr_test => true, :dimensions => [n1,n2,n3])
conv_filter = ConvolutionFilter::new('sfrf',SYM8_MF,8)
ksym8 = MagicFilter1d( conv_filter, optims )
puts ksym8

input = ANArray.float(256,n1,n2,n3).random!
work1 = ANArray.float(256,n1,n2,n3)
work2 = ANArray.float(256,n1,n2,n3)
output = ANArray.float(256,n1,n2,n3)

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3

ksym8.build(:openmp => true)
stats0 = ksym8.run(3, 0, n, BC::PERIODIC, n, n, input, work1)
stats0 = ksym8.run(3, 0, n, BC::PERIODIC, n, n, input, work1)
stats1 = ksym8.run(3, 1, n, BC::PERIODIC, n, n, work1, work2)
stats2 = ksym8.run(3, 2, n, BC::PERIODIC, n, n, work2, output)
puts "#{ksym8.procedure.name} d0: #{stats0[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats0[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name} d1: #{stats1[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats1[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name} d2: #{stats2[:duration]*1.0e3} ms #{32*n1*n2*n3 / (stats2[:duration]*1.0e9)} GFlops"
puts "#{ksym8.procedure.name}: #{(stats0[:duration]+stats1[:duration]+stats2[:duration])*1.0e3} ms #{3*32*n1*n2*n3 / ((stats0[:duration]+stats1[:duration]+stats2[:duration])*1.0e9)} GFlops"
