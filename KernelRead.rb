require "./BOAST.rb"
require 'narray'
module ConvolutionGenerator
  def ConvolutionGenerator::kernel_read_ref( unrolled = 1, size = 4 )
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::C)
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::C
    function_name = "kernel_read_ref"
    function_name += "_#{unrolled}l_#{size*8}b"
    m_start = Variable::new("m_start",Int,{:direction => :in})
    m_cycles = Variable::new("m_cycles",Int,{:direction => :in})
    m_stride = Variable::new("m_stride",Int,{:direction => :in})
    buffer_size = Variable::new("buffer_size",Int,{:direction => :in})
    sum = Variable::new("sum",Int);
    buffer = Variable::new("buffer", Int, {:direction => :in, :size => size, :dimension => [ Dimension::new(0,buffer_size-1) ]})
    i = Variable::new("i",Int);
    j = Variable::new("j",Int);
    p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], [], {:return => sum } ) {
      i.decl
      j.decl
      sum.decl
      (sum === 0).print
      For::new(i, 1, m_cycles*m_stride) {
        For::new(j, m_start, buffer_size+m_start-1, m_stride*unrolled) {
          unrolled.times { |k|
            (sum === sum + buffer[j+m_stride*k]).print
          }
        }.print
      }.print
    }
    p.print 
    kernel.procedure = p
    return kernel
  end

end

