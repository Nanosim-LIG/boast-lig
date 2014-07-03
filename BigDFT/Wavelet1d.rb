require './GenericConvolution.rb'
require './WaveletFilters.rb'
def Wavelet1d(wavelet_filter, direction, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(wavelet_filter, :wavelet => direction, :a => true, :a_y => true, :ld => true, :narr => true)
  conv_operation.optimize(optims) if optims

  p, subops = conv_operation.procedure

  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op|
    BOAST::print op
    puts "chosen:"+ op.name
  }
  BOAST::print p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
  return kernel
end

optims = nil #GenericOptimization::new()#:unroll_range => 6, :mod_arr_test => true, :tt_arr_test => true)
wave_filter = WaveletFilter::new("sym#{SYM8.length/2}", SYM8)
ksym8 = Wavelet1d( wave_filter, :decompose, optims )
ksym8.print
ksym8i = Wavelet1d( wave_filter, :recompose, optims )
ksym8i.print
