require './GenericConvolution.rb'
module BOAST
  def BOAST::Wavelet(wavelet_filter, direction, optims=GenericOptimization::new)
    conv_operation = GenericConvolutionOperator::new(wavelet_filter, :transpose => 1, :work => true, :wavelet => direction)
    conv_operation.optimize(optims)

    p, subops= conv_operation.procedure()

    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang
    if BOAST::get_lang == C then
      @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      @@output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      @@output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end

    subops.each_value { |op|
      print op
      puts "chosen:"+ op.name
    }
    print p

    kernel.procedure = p
    kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
    return kernel
  end
end
