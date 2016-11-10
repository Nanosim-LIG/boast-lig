require "BOAST"
require 'narray'

include BOAST

def kernel_read_vectorized( unrolled = 1, elem_size = 8, length = 2, machine = "sse3")
  push_env(:lang => C)
  kernel = CKernel::new
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
  elsif machine == "avx" then
    type_name = "__m256d"
    header="immintrin.h"
  end
  resultV=nil
  if machine == "sse3" or machine == "neon" then
    resultV = Variable::new("resultV",Int, {:size => 8, :dimension => [Dimension::new(0,1)], :local => true})
  elsif machine == "avx" then
    resultV = Variable::new("resultV",Real, {:size => 8, :dimension => [Dimension::new(0,3)], :local => true})
  end
  sum = Variable::new("sum",Int)
  sumV = Variable::new("sumV",CustomType, {:type_name => type_name, :size => elem_size*length})
  buffer = Variable::new("buffer", CustomType, {:type_name => type_name, :direction => :in, :size => elem_size*length, :dimension => [ Dimension::new(0,buffer_size-1) ]})
  i = Variable::new("i",Int)
  j = Variable::new("j",Int)
  get_output.print "#include <#{header}>\n"
  get_output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
  p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], :return => sum , :headers => [header]) {
    i.decl
    j.decl
    resultV.decl
    sumV.decl
    if machine == "sse3" then
      get_output.print "  sumV = _mm_set_epi32(0,0,0,0);\n"
    elsif machine == "neon" then
      get_output.print "  sumV = vmovq_n_s64(0);\n"
    elsif machine == "avx" then
      get_output.print "  _mm256_set1_pd(0);\n"
    end
    pr For::new(i, 1, m_cycles*m_stride) {
      pr For::new(j, m_start, buffer_size + m_start - m_stride*unrolled, step: m_stride*unrolled) {
        unrolled.times { |k|
           if machine == "sse3" then
             pr sumV === FuncCall::new("_mm_add_epi64", sumV, buffer[j+m_stride*k])
           elsif machine == "neon" then
             pr sumV === FuncCall::new("vaddq_s64", sumV, buffer[j+m_stride*k])
           elsif machine == "avx" then
             pr sumV === FuncCall::new("_mm256_add_pd", sumV, buffer[j+m_stride*k])
           end
              }
      }
      pr For::new(j, buffer_size + m_start - FuncCall::new( "modulo", buffer_size+m_start, m_stride*unrolled),  buffer_size + m_start - 1, step: m_stride) {
         if machine == "sse3" then
           pr sumV === FuncCall::new("_mm_add_epi64", sumV, buffer[j])
         elsif machine == "neon" then
           pr sumV === FuncCall::new("vaddq_s64", sumV, buffer[j])
         elsif machine == "avx" then
           pr sumV === FuncCall::new("_mm256_add_pd", sumV, buffer[j])
         end
      }
    }
    if machine == "sse3" then
      get_output.print "  _mm_store_si128((__m128i *) resultV, sumV);\n"
      pr sum === resultV[0] + resultV[1]
    elsif machine == "neon" then
      get_output.print "  vst1q_s64(resultV, sumV);\n"
      pr sum === resultV[0] + resultV[1]
    elsif machine == "avx" then
      get_output.print "  _mm256_store_pd(resultV, sumV);\n"
      pr sum === resultV[0] + resultV[1] + resultV[2] + resultV[3]
    end
  }
  pr p 
  kernel.procedure = p
  pop_env(:lang)
  return kernel
end
def kernel_read_ref( unrolled = 1, size = 4 )
  push_env(:lang => C)
  kernel = CKernel::new
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
  get_output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
  p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], :return => sum ) {
    i.decl
    j.decl
    sum.decl
    (sum === 0).print
    pr For::new(i, 1, m_cycles*m_stride) {
      pr For::new(j, m_start, buffer_size + m_start - m_stride*unrolled, step: m_stride*unrolled) {
        unrolled.times { |k|
          pr sum === sum + buffer[j+m_stride*k]
        }
      }
      pr For::new(j, buffer_size + m_start - FuncCall::new( "modulo", buffer_size+m_start, m_stride*unrolled),  buffer_size + m_start - 1, step: m_stride) {
          pr sum === sum + buffer[j]
      }
    }
  }
  pr p 
  pop_env(:lang)
  kernel.procedure = p
  return kernel
end

