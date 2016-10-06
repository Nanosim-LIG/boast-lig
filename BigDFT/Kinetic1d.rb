require './GenericConvolution.rb'
require './WaveletFilters.rb'

openmp = true

def Kinetic1d(filter, optims=GenericOptimization::new)
  conv_operation = GenericConvolutionOperator1d::new(filter, :kinetic => true, :ld => true, :narr => true, :a_x => true,:a_y=>true, :a => true, :dot_in => true)
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

conv_filter = ConvolutionFilter::new('sym8_d2',SYM8_D2,14)
kkin8 = Kinetic1d( conv_filter, optims )

puts kkin8
kkin8.build(:openmp => openmp)
kkin8.dump_source
