require './GenericConvolution.rb'
def Wavelet(wavelet_filter, direction, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator::new(wavelet_filter, :transpose => 1, :work => true, :wavelet => direction)
  conv_operation.optimize(optims)

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
