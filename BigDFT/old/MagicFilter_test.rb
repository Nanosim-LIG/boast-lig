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
work = NArray.float(n1,n2,n3)
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
epsilon = 10e-15

k = BOAST::magicfilter_per_ref
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n2, n1*n3, output_ref, work)
stats = k.run(n3, n2*n1, work, output_ref)

puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(2) { |i|
  k = BOAST::MagicFilter(FILTER,8,i,false)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n2, n1*n3, output, work)
  stats = k.run(n3, n2*n1, work, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
1.upto(1) { |i|
  k = BOAST::MF3d(FILTER,8,i)
  k.build(:openmp => true)
  begin
    stats = k.run(n1, n2, n3, input, output, work)
    stats = k.run(n1, n2, n3, input, output, work)
  rescue Exception => e
    puts e.inspect
  end
  diff = (output_ref - output).abs
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
#dfgdfag
k = BOAST::magicfilter_per_ref(true)
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n2, n1*n3, output_ref, work)
stats = k.run(n3, n2*n1, work, output_ref)

puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

#BOAST::MagicFilter(FILTER,8,0,false).run(32, 32, " "*32*32*8, " "*32*32*8)
1.upto(1) { |i|
  k = BOAST::MF3d(FILTER.reverse,7,i)
  #k.print
  k.build(:openmp => true)
  begin
    stats = k.run(n1, n2, n3, input, output, work)
    stats = k.run(n1, n2, n3, input, output, work)
  rescue Exception => e
    puts e.inspect
  end
  diff = (output_ref - output).abs
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }

#  k = BOAST::MagicFilter(FILTER,8,i,true)
#  stats = k.run(n1, n2*n3, input, output)
#  stats = k.run(n1, n2*n3, input, output)
#  diff = (output_ref - output).abs
#  diff.each { |elem|
#    raise "Warning: residue too big: #{elem}" if elem > epsilon
#  }
#  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

###here follow the Free BC cases

input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n2,n3,n1+15)
output = NArray.float(n2,n3,n1+15)
k = BOAST::magicfilter_per_ref(false, true)
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(2) { |i|
  k = BOAST::MagicFilter(FILTER,8,i,false,true)
  k.build
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs
  
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
#redo the 3d case for the grow operation
input = NArray.float(n1,n2,n3).random
work = NArray.float(n1+15,n2+15,n3+15)
output_ref = NArray.float(n1+15,n2+15,n3+15)
output = NArray.float(n1+15,n2+15,n3+15)
epsilon = 10e-15
puts 'again, grow'
k = BOAST::magicfilter_per_ref(false,true)
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n1, n2*n3, input, output_ref)
stats = k.run(n2, (n1+15)*n3, output_ref, work)
stats = k.run(n3, (n2+15)*(n1+15), work, output_ref)

puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(1) { |i|
  k = BOAST::MagicFilter(FILTER,8,i,false,true)
  #k.print
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n2, n3*(n1+15), output, work)
  stats = k.run(n3, (n2+15)*(n1+15), work, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
#rtret
1.upto(1) { |i|
  k = BOAST::MF3d(FILTER,8,i,[BOAST::BC::GROW]*3)
  #k.print
  k.build(:openmp => true)
  begin
    stats = k.run(n1, n2, n3, input, output, work)
    stats = k.run(n1, n2, n3, input, output, work)
  rescue Exception => e
    puts e.inspect
  end
  diff = (output_ref - output).abs
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}



input = NArray.float(n1+15,n2,n3).random
output_ref = NArray.float(n2,n3,n1)
output = NArray.float(n2,n3,n1)
k = BOAST::magicfilter_per_ref(true, true)
stats = k.run(n1, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(2) { |i|
  k = BOAST::MagicFilter(FILTER,8,i,true,true)
  stats = k.run(n1, n2*n3, input, output)
  stats = k.run(n1, n2*n3, input, output)
  diff = (output_ref - output).abs

  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}

#redo the 3d case for the shrink operation
input = NArray.float(n1+15,n2+15,n3+15).random
work = NArray.float(n1+15,n2+15,n3+15)
output_ref = NArray.float(n1,n2+15,n3+15)
output = NArray.float(n1,n2+15,n3+15)
epsilon = 10e-15
puts 'again, shrink'
k = BOAST::magicfilter_per_ref(true,true)
stats = k.run(n1, (n2+15)*(n3+15), input, output_ref)
stats = k.run(n1, (n2+15)*(n3+15), input, output_ref)
stats = k.run(n2, n1*(n3+15), output_ref, work)
stats = k.run(n3, n2*n1, work, output_ref)

puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
1.upto(1) { |i|
  k = BOAST::MagicFilter(FILTER,8,i,true,true)
  stats = k.run(n1, (n2+15)*(n3+15), input, output)
  stats = k.run(n1, (n2+15)*(n3+15), input, output)
  stats = k.run(n2, n1*(n3+15), output, work)
  stats = k.run(n3, n2*n1, work, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
#rtret
1.upto(6) { |i|
  k = BOAST::MF3d(FILTER.reverse,7,i,[BOAST::BC::SHRINK]*3)
  #k.print
  k.build(:openmp => true)
  begin
    stats = k.run(n1, n2, n3, input, output, work)
    stats = k.run(n1, n2, n3, input, output, work)
  rescue Exception => e
    puts e.inspect
  end
  diff = (output_ref - output).abs
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
