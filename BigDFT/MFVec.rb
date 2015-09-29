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

def comp_tile_n1( n1 = 1, n2 = 1, vector_length = 1, real_size = 8 )
  raise "Invalid vector size" if n1%vector_length != 0
  arr = BOAST::ConstArray::new(FILTER, BOAST::Real)
  filter = Real( "filt", :const => arr, :dim => Dim(0,FILTER.length-1) )
  n1_stride_input = Int( "n1_stride_in", :dir => :in)
  n1_stride_output = Int( "n1_stride_out", :dir => :in)
  input =  Real( "input", :dir => :in, :align => vector_length * real_size, :dim => Dim() )
  output = Real( "output", :dir => :out, :align => vector_length * real_size, :dim => Dim() )
  p = Procedure( "comp_n1_#{n1}_#{n2}_#{vector_length}_#{real_size}", [n1_stride_input, n1_stride_output, input, output], [filter], :inline => true, :local => true ) {
    f = Real( "f", :vector_length => vector_length )
    decl f
    tt = (n1/vector_length).times.collect { |i|
      n2.times.collect { |j|
        Real( "tt#{i}#{j}", :vector_length => vector_length )
      }
    }
    xx = (n1/vector_length).times.collect { |i|
      n2.times.collect { |j|
        Real( "xx#{i}#{j}", :vector_length => vector_length )
      }
    }
    decl *tt.flatten
    decl *xx.flatten
    (n1/vector_length).times.collect { |i|
      n2.times.collect { |j|
        pr tt[i][j].set( 0.0 )
      }
    }
    FILTER.length.times { |k|
      pr f.set(filter[k])
      (n1/vector_length).times.collect { |i|
        n2.times.collect { |j|
          pr xx[i][j] === input[n1_stride_input*j+vector_length*i+k]
          pr tt[i][j] === tt[i][j] + xx[i][j]*f
        }
      }
    }
    (n1/vector_length).times.collect { |i|
      n2.times.collect { |j|
        out_index = output[n1_stride_output*j+vector_length*i]
        out_index.align = vector_length * real_size
        pr out_index === tt[i][j]
      }
    }
  }
  return p 
end

def mfvec( n1 = 1, n2 = 1, vector_length = 1, real_size = 8 )

  kern = CKernel::new
  input  = Real( "input",  :dir => :in, :align => 32, :dim => [Dim(128+16), Dim(128*128)] )
  output = Real( "output", :dir => :out, :align => 32, :dim => [Dim(128), Dim(128*128)] )

  sub = comp_tile_n1( n1, n2, vector_length, real_size)
  get_output.puts "#include <immintrin.h>"

  p = Procedure( "mfs", [input,output] ) {
    i = Int "i"
    j = Int "j"
    tt = Real "tt"
    decl i,j
    decl tt
    pr For( j, 1, 128*128, :step => n2 ) {
      pr For( i, 1, 128, :step => n1 ) {
        pr sub.call(128+16, 128, input[i, j].address, output[i,j].address)
      }
    }
  }

  pr sub
  pr p
  kern.procedure = p

  return kern

end

k = mfvec(8,1,4)

puts k
k.build
puts k.maqao_analysis

