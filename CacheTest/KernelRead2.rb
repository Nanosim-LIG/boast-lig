require "BOAST"
require 'narray'

module BOAST
  def BOAST::kernel_read_vectorized2( unrolled = 1, elem_size = 8, length = 2)
    push_env( :array_start => 0 )
    push_env( :lang => C )
    kernel = CKernel::new

    function_name = "kernel_read_ref"
    function_name += "_#{unrolled}l_#{elem_size*8}b#{length}v"
    m_start = Int("m_start", {:direction => :in})
    m_cycles = Int("m_cycles", {:direction => :in})
    m_stride = Int("m_stride", {:direction => :in})
    buffer_size = Int("buffer_size", {:direction => :in})
    header="immintrin.h"

    buffer = Int("buffer", { :direction => :in, :size => elem_size, :vector_length => length, :dimension => [ Dim(buffer_size) ]})
    @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
    sum = Int("sum")
    print p = Procedure::new(function_name, [m_start, m_cycles, m_stride, buffer_size, buffer], [], {:return => sum , :headers => [header]}) {

      resultV = Int("resultV", {:size => elem_size, :dimension => [Dim(length)], :local => true})
      sumV = Int("sumV", {:size => elem_size, :vector_length => length})
      vec = Int("vec", {:size => elem_size, :vector_length => length})
      i = Int("i")
      j = Int("j")
      i.decl
      j.decl
      resultV.decl
      sumV.decl
      sum.decl

      (0...length).each { |indx|
        print resultV[indx] === 0
      }

      print sumV === resultV.dereference

      print For(i, 1, m_cycles*m_stride) {
        print For(j, m_start, buffer_size + m_start - m_stride*unrolled, m_stride*unrolled) {
          unrolled.times { |k|
            print sumV === sumV + buffer[j+m_stride*k]
          }
        }
        print For(j, buffer_size + m_start - modulo(buffer_size+m_start, m_stride*unrolled),  buffer_size + m_start - 1, m_stride) {
          print sumV === sumV + buffer[j]
        }
      }

      print resultV.dereference === sumV

      print sum === 0
      (0...length).each { |indx|
        print sum === sum + resultV[indx]
      }
    }
    pop_env( :lang )
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
