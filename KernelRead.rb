require "./BOAST.rb"
require 'narray'

module ConvolutionGenerator
  def ConvolutionGenerator::kernel_read_vectorized( unrolled = 1, elem_size = 8, length = 2, machine = "sse3")
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::C)
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::C
    function_name = "kernel_read_ref"
    function_name += "_#{unrolled}l_#{elem_size*8}b#{length}v"
    m_start = Variable::new("m_start",Int,{:direction => :in})
    m_cycles = Variable::new("m_cycles",Int,{:direction => :in})
    m_stride = Variable::new("m_stride",Int,{:direction => :in})
    buffer_size = Variable::new("buffer_size",Int,{:direction => :in})
    type_name=nil
    header=nil
    if machine == "sse3" then
      type_name = "__m128i" 
      header="immintrin.h"
    elsif machine == "neon" then
      type_name = "int64x2_t"
      header="arm_neon.h"
    end
    resultV = Variable::new("resultV",Int, {:size => 8, :dimension => [Dimension::new(0,1)], :local => true})
    sum = Variable::new("sum",Int)
    sumV = Variable::new("sumV",CustomType, {:type_name => type_name, :size => elem_size*length})
    buffer = Variable::new("buffer", CustomType, {:type_name => type_name, :direction => :in, :size => elem_size*length, :dimension => [ Dimension::new(0,buffer_size-1) ]})
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    $output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
    p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], [], {:return => sum , :headers => [header]}) {
      i.decl
      j.decl
      resultV.decl
      sumV.decl
      sum.decl
      if machine == "sse3" then
        $output.print "  sumV = _mm_set_epi32(0,0,0,0);\n"
      elsif machine == "neon" then
             $output.print "  sumV = vmovq_n_s64(0);\n"
      end
      For::new(i, 1, m_cycles*m_stride) {
        For::new(j, m_start, buffer_size + m_start - m_stride*unrolled, m_stride*unrolled) {
          unrolled.times { |k|
             if machine == "sse3" then
               (sumV === FuncCall::new("_mm_add_epi64", sumV, buffer[j+m_stride*k])).print
             elsif machine == "neon" then
                    (sumV === FuncCall::new("vaddq_s64", sumV, buffer[j+m_stride*k])).print
             end
                }
        }.print
        For::new(j, buffer_size + m_start - FuncCall::new( "modulo", buffer_size+m_start, m_stride*unrolled),  buffer_size + m_start - 1, m_stride) {
           if machine == "sse3" then
             (sumV === FuncCall::new("_mm_add_epi64", sumV, buffer[j])).print
           elsif machine == "neon" then
                  (sumV === FuncCall::new("vaddq_s64", sumV, buffer[j])).print
           end
        }.print
      }.print
      if machine == "sse3" then
        $output.print "  _mm_store_si128((__m128i *) resultV, sumV);\n"
      elsif machine == "neon" then
             $output.print "  vst1q_s64(resultV, sumV);\n"
      end
      (sum === resultV[0] + resultV[1]).print
    }
    p.print 
    kernel.procedure = p
    return kernel
  end
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
    sum = Variable::new("sum",Int)
    buffer = Variable::new("buffer", Int, {:direction => :in, :size => size, :dimension => [ Dimension::new(0,buffer_size-1) ]})
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    $output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
    p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], [], {:return => sum } ) {
      i.decl
      j.decl
      sum.decl
      (sum === 0).print
      For::new(i, 1, m_cycles*m_stride) {
        For::new(j, m_start, buffer_size + m_start - m_stride*unrolled, m_stride*unrolled) {
          unrolled.times { |k|
            (sum === sum + buffer[j+m_stride*k]).print
          }
        }.print
        For::new(j, buffer_size + m_start - FuncCall::new( "modulo", buffer_size+m_start, m_stride*unrolled),  buffer_size + m_start - 1, m_stride) {
            (sum === sum + buffer[j]).print
        }.print
      }.print
    }
    p.print 
    kernel.procedure = p
    return kernel
  end

end

