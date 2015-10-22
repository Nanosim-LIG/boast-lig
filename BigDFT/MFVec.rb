require 'BOAST'
include BOAST

FILTER = [ "8.4334247333529341094733325815816e-7",
          "-0.1290557201342060969516786758559028e-4",
           "0.8762984476210559564689161894116397e-4",
          "-0.30158038132690463167163703826169879e-3",
           "0.174723713672993903449447812749852942e-2",
          "-0.942047030201080385922711540948195075e-2",
           "0.2373821463724942397566389712597274535e-1",
           "0.612625895831207982195380597e-1",
           "0.9940415697834003993178616713",
          "-0.604895289196983516002834636e-1",
          "-0.2103025160930381434955489412839065067e-1",
           "0.1337263414854794752733423467013220997e-1",
          "-0.344128144493493857280881509686821861e-2",
           "0.49443227688689919192282259476750972e-3",
          "-0.5185986881173432922848639136911487e-4",
           "2.72734492911979659657715313017228e-6" ]

set_lang( C )
#set_replace_constants(false)
set_boast_inspect( true )
set_array_start(0)

def comp_tile_n1( n1 = 1, n2 = 1, vector_length = 1, real_size = 8 )
  raise "Invalid vector size" if n1%vector_length != 0
  arr = BOAST::ConstArray::new(FILTER, BOAST::Real)
  filter = Real( "filt", :const => arr, :dim => Dim(0,FILTER.length-1) )
  n1_stride_input = Int( "n1_stride_in", :dir => :in)
  n1_stride_output = Int( "n1_stride_out", :dir => :in)
  input =  Real( "input", :dir => :in, :align => vector_length * real_size, :dim => Dim() )
  output = Real( "output", :dir => :out, :align => vector_length * real_size, :dim => Dim() )
  return (1..n1).collect { |n|
    Procedure( "comp_n1_#{n}_#{n2}_#{vector_length}_#{real_size}", [n1_stride_input, n1_stride_output, input, output], [filter], :inline => true, :local => true ) {
      f = Real( "f", :vector_length => vector_length )
      decl f
      tt = ((n+vector_length-1)/vector_length).times.collect { |i|
        n2.times.collect { |j|
          Real( "tt#{i}#{j}", :vector_length => vector_length )
        }
      }
      xx = ((n+vector_length-1)/vector_length).times.collect { |i|
        n2.times.collect { |j|
          Real( "xx#{i}#{j}", :vector_length => vector_length )
        }
      }
      decl *tt.flatten
      decl *xx.flatten

      tt.flatten.each { |t| pr t.set( 0.0 ) }
      FILTER.length.times { |k|
        pr f.set(filter[k])
        ((n+vector_length-1)/vector_length).times.collect { |i|
          n2.times.collect { |j|
            if (n-i*vector_length) / vector_length > 0 then
              pr xx[i][j] === input[n1_stride_input*j+vector_length*i+k]
            else
              pr xx[i][j] === MaskLoad(input[n1_stride_input*j+vector_length*i+k], [true]*(n-i*vector_length)+[false]*((i+1)*vector_length-n), xx[i][j])
            end
            pr tt[i][j] === FMA(xx[i][j], f, tt[i][j])
          }
        }
      }
      ((n+vector_length-1)/vector_length).times.collect { |i|
        n2.times.collect { |j|
          out_index = output[n1_stride_output*j+vector_length*i]
          out_index.align = vector_length * real_size
          if (n-i*vector_length) / vector_length > 0 then
            pr out_index === tt[i][j]
          else
            pr MaskStore( out_index, tt[i][j], [-1]*(n-i*vector_length)+[0]*((i+1)*vector_length-n))
          end
        }
      }
    }
  }
end

def mfvec( n1 = 1, n2 = 1, n3 = 1, vector_length = 1, real_size = 8 )

  input_n1_size = 128+16
  output_n1_size = 128
  input_n2_size = 128+16
  output_n2_size = 128
  input_n3_size = 128+16
  output_n3_size = 128
  input_n1_stride = input_n1_size
  output_n1_stride = output_n1_size
  input_n2_stride = input_n2_size
  output_n2_stride = output_n2_size
  input_n3_stride = input_n3_size
  output_n3_stride = output_n3_size

  kern = CKernel::new
  input     = Real( "input",     :dir => :in,  :align => 32, :dim => [ Dim(input_n1_stride),  Dim(input_n2_stride *      input_n3_stride) ] )
  output_d1 = Real( "output_d1", :dir => :out, :align => 32, :dim => [ Dim(output_n1_stride), Dim(input_n2_stride *      input_n3_stride) ] )
  input_d2  = Real( "output_d1", :dir => :out, :align => 32, :dim => [ Dim(output_n1_stride), Dim(input_n2_stride),  Dim(input_n3_stride) ] )
  output_d2 = Real( "output_d2", :dir => :out, :align => 32, :dim => [ Dim(output_n1_stride), Dim(output_n2_stride), Dim(input_n3_stride) ] )
  input_d3  = Real( "output_d2", :dir => :out, :align => 32, :dim => [ Dim(output_n1_stride), Dim(output_n2_stride), Dim(input_n3_stride) ] )
  output    = Real( "output",    :dir => :out, :align => 32, :dim => [ Dim(output_n1_stride), Dim(output_n2_stride), Dim(output_n3_stride)] )

  subs_n1 = comp_tile_n1( n1, 1, vector_length, real_size)
  subs_n2 = comp_tile_n2( n2, vector_length, real_size)
  subs_n3 = comp_tile_n2( n3, vector_length, real_size)
  get_output.print <<EOF
#include <immintrin.h>
static inline __m128i _mm_setr_epi64x(long long __q1, long long __q0)
{
  return _mm_set_epi64x(__q0, __q1);
}
EOF

  p = Procedure( "mfsd1", [input,output] ) {
    i = Int "i"
    j = Int "j"
    tt = Real "tt"
    decl i,j,k
    decl tt
    conditions_lambdas_n1 = (1...n1).to_a.reverse.collect { |i|
      [ i, lambda { pr subs_n1[i-1].call(input_n1_stride, output_n1_stride, input[output_n1_size-1-i,j].address, output_d1[output_n1_size-1-i,j].address ) } ]
    }
    conditions_lambdas_n2 = (1...n2).to_a.reverse.collect { |k|
      [ i, lambda {
        pr For( i, 0, output_n2_size-1 ) {
          pr subs_n2[k-1].call(output_n1_stride, output_n1_stride, input_d2[output_n1_size-1-k,i,j].address, output_d2[output_n1_size-1-k,i,j].address )
        }
      } ]
    }
    conditions_lambdas_n3 = (1...n3).to_a.reverse.collect { |k|
      [ i, lambda {
        pr For( i, 0, output_n2_size-1 ) {
          pr subs_n3[i-1].call(input_n1_stride, output_n1_stride, input_d2[output_n1_size-1-k,j,i].address, output_d2[output_n1_size-1-k,j,i].address )
        }
      } ]
    }
    conditions_lambdas.flatten!
    pr For( j, 0, output_n2_size*output_n3_size-1, :step => 1 ) {
      pr For( i, 0, output_n1_size-n1, :step => n1 ) {
        pr subs_n1[-1].call(input_n1_stride, output_n1_stride, input[i, j].address, output_d1[i,j].address)
      }
      pr Case(output_n1_size % n1, *conditions_lambdas_n1)
    }
    pr For( j, 0, output_n3_size-1 ) {
      pr For( k, 0, output_n1_size-1, :step => n2 ) {
        pr For( i, 0, output_n2_size-1 ) {
          pr subs_n2[-1].call(output_n1_stride, output_n1_stride, input_d2[k, i, j], output_d2[k, i, j])
        }
      }
      pr Case(output_n1_size % n2, *conditions_lambdas_n2)
    }
    pr For( j, 0, output_n2_size-1) {
      pr For( k, 0, output_n1_size-1, :step => n3 ) {
        pr For( i, 0, output_n3_size-1 ) {
          pr subs_n3[-1].call(output_n1_stride*output_n2_stride, output_n1_stride*output_n2_stride, input_d3[k, j, i], output_d3[k, j, i])
        }
      }
      pr Case(output_n1_size % n3, *conditions_lambdas_n3)
    }
  }

  subs.each { |sub| pr sub }
  pr p
  kern.procedure = p

  return kern

end

k = mfvec(8,8,8,4)

puts k
k.build
puts k.maqao_analysis

