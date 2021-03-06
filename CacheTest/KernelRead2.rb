require "BOAST"
require 'narray'
include BOAST
def kernel_read_vectorized2( unrolled = 1, elem_size = 8, length = 2)
  push_env( :array_start => 0 )
  push_env( :lang => C )
  kernel = CKernel::new

  function_name = "kernel_read_ref"
  function_name += "_#{unrolled}l_#{elem_size*8}b#{length}v"
  m_start = Int("m_start", {:direction => :in})
  m_cycles = Int("m_cycles", {:direction => :in})
  m_stride = Int("m_stride", {:direction => :in})
  buffer_size = Int("buffer_size", {:direction => :in})

  buffer = Int("buffer", { :direction => :in, :size => elem_size, :vector_length => length, :dimension => [ Dim(buffer_size) ]})
  header="immintrin.h"
  get_output.print "#include <immintrin.h>\n"
  get_output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
  register_funccall("modulo")
  sum = Int("sum")
  pr p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], :return => sum, :headers => [header]) {

    resultV = Int("resultV", {:size => elem_size, :dimension => [Dim(length)], :local => true})
    sumV = Int("sumV", {:size => elem_size, :vector_length => length})
    vec = Int("vec", {:size => elem_size, :vector_length => length})
    decl i = Int("i")
    decl j = Int("j")
    decl elem_n = Int("elem_n")
    decl resultV
    decl sumV

    pr elem_n === buffer_size / (elem_size*length)
#      (0...length).each { |indx|
#        print resultV[indx] === 0
#      }


    pr Set(sumV, 0)

    pr For(i, 1, m_cycles*m_stride) {
      pr For(j, m_start, elem_n + m_start - m_stride*unrolled, step: m_stride*unrolled) {
        unrolled.times { |k|
          pr sumV === sumV + buffer[j+m_stride*k]
        }
      }
      pr For(j, elem_n + m_start - modulo(elem_n+m_start, m_stride*unrolled),  elem_n + m_start - 1, step: m_stride) {
        pr sumV === sumV + buffer[j]
      }
    }

    rV = resultV.dereference
    rV.alignment = length * elem_size
    pr rV === sumV

    pr sum === 0
    (0...length).each { |indx|
      pr sum === sum + resultV[indx]
    }
  }
  pop_env( :lang )
  pop_env( :array_start )
  kernel.procedure = p
  return kernel
end
