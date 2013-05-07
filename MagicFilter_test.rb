require'./MagicFilter.rb'

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



n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
epsilon = 10e-15



ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
k = ConvolutionGenerator::magicfilter_per_ref
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false).run(32, 32, " "*32*32*8, " "*32*32*8)
1.upto(2) { |i|
  k = ConvolutionGenerator::MagicFilter(FILTER,8,i,false)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

k = ConvolutionGenerator::magicfilter_per_ref(true)
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false).run(32, 32, " "*32*32*8, " "*32*32*8)
1.upto(2) { |i|
  k = ConvolutionGenerator::MagicFilter(FILTER,8,i,true)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n2,n3,n1+15)
output = NArray.float(n2,n3,n1+15)
k = ConvolutionGenerator::magicfilter_per_ref(false, true)
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(2) { |i|
  k = ConvolutionGenerator::MagicFilter(FILTER,8,i,false,true)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

input = NArray.float(n1+15,n2,n3).random
output_ref = NArray.float(n2,n3,n1)
output = NArray.float(n2,n3,n1)
k = ConvolutionGenerator::magicfilter_per_ref(true, true)
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(2) { |i|
  k = ConvolutionGenerator::MagicFilter(FILTER,8,i,true,true)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    puts "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
