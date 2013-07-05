require "./BOAST.rb"
require 'narray'
module ConvolutionGenerator
  def ConvolutionGenerator::kernel_read_ref( unrolled = false, size = 4 )
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::C)
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::C
    function_name = "kernel_read_ref"
    function_name += "_unroll" if unrolled
    m_start = Variable::new("m_start",Int,{:direction => :in})
    m_cycles = Variable::new("m_cycles",Int,{:direction => :in})
    m_stride = Variable::new("m_stride",Int,{:direction => :in})
    buffer_size = Variable::new("buffer_size",Int,{:direction => :in})
    buffer = Variable::new("buffer", Int, {:direction => :in, :size => size, :dimension => [ Dimension::new(0,buffer_size-1) ]})
    i = Variable::new("i",Int);
    j = Variable::new("j",Int);
    sum = Variable::new("sum",Int);
    p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer]) {
      i.decl
      j.decl
      sum.decl
      For::new(i, 1, m_cycles*m_stride) {
        For::new(j, m_start, buffer_size+m_start-1, m_stride) {
          (sum === sum + buffer[j]).print
        }.print
      }.print
    }.print 
    
    kernel.procedure = p
    return kernel
  end

end

k = ConvolutionGenerator::kernel_read_ref
puts k.print
