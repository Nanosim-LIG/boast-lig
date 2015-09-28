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
set_replace_constants(false)

def mfvec(use_sub = false)
  tile_x = 1
  tile_y = 1

  kern = CKernel::new
  arr = BOAST::ConstArray::new(FILTER, BOAST::Real)
  filter = Real( "filt", :const => arr, :align => 32, :dim => Dim(0,15) )
  input  = Real( "input",  :dir => :in, :align => 32, :dim => [Dim(128+16), Dim(128*128)] )
  output = Real( "output", :dir => :out, :align => 32, :dim => [Dim(128), Dim(128*128)] )

  sub_input = Real( "sub_input", :dir => :in, :align => 32, :dim => [Dim(0,15)] )
  sub_tt = Real( "sub_tt", :dir => :out )
  sub = Procedure( "comp1", [sub_tt, sub_input], [filter], :inline => true, :local => true ) {
    pr sub_tt === 0.0
    k = Int "k"
    For( k, 0, 15 ) {
      pr sub_tt === sub_tt + sub_input[k]*filter[k]
    }.unroll
  }

  p = Procedure( "mfs", [input,output], [filter] ) {
    i = Int "i"
    j = Int "j"
    tt = Real "tt"
    decl i,j
    decl tt
    pr For( j, 1, 128*128, :step => tile_y ) {
      pr For( i, 1, 128 ) {
        if use_sub then
          pr sub.call(tt.address, input[i, j].address)
        else
          k = Int "k"
          decl k
          pr tt === 0.0
          For( k, 0, 15 ) {
            pr tt === tt + input[i+k,j]*filter[k]
          }.unroll
        end
        pr output[i, j] === tt
      }
    }
  }

  pr sub
  pr p
  kern.procedure = p

  return kern

end

k = mfvec

puts k
k.build
puts k.maqao_analysis

