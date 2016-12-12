require './GenericConvolution2.rb'
require './WaveletFilters.rb'
def Wavelet1d(wavelet_filter, direction, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(wavelet_filter, :wavelet => direction, :a => true, :a_y => true, :ld => true, :narr => true)
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

optims = GenericOptimization::new(:unroll_range => [6,6], :mod_arr_test => true, :tt_arr_test => true)
wave_filter_decompose = WaveletFilterDecompose::new("sym#{SYM8_LP.length/2}", SYM8_LP)
wave_filter_recompose = WaveletFilterRecompose::new("sym#{SYM8_LP.length/2}", SYM8_LP)
ksym8 = Wavelet1d( wave_filter_decompose, :decompose, optims )
puts ksym8
ksym8i = Wavelet1d( wave_filter_recompose, :recompose, optims )
puts ksym8i
