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

m = 2
n1 = 124
n2 = 132
n3 = 130
pad_in = 0
pad_out = 2
input = NArray.float(n1+pad_in,n2+pad_in,n3+pad_in,m).random
(1..(m-1)).each { |i|
  input[ true, true, true, i] = input[ true, true, true, 0]
}
work1 = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
work2 = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
epsilon = 10e-15

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3

ld_in = NArray.int(3)
ld_in[0] = n1+pad_in
ld_in[1] = n2+pad_in
ld_in[2] = n3+pad_in

ld_out = NArray.int(3)
ld_out[0] = n1+pad_out
ld_out[1] = n2+pad_out
ld_out[2] = n3+pad_out

bc = NArray.int(3)

k_ref = BOAST::magicfilter_per_ref
stats = k_ref.run(n1, n2*n3, input, output_ref)
stats = k_ref.run(n2, n1*n3, output_ref, work1)
stats = k_ref.run(n3, n2*n1, work1, output_ref)

puts "#{k_ref.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

conv_filter = BOAST::ConvolutionFilter::new('sfrf',FILTER,8)
optims = BOAST::GenericOptimization::new(:unroll_range => 6, :mod_arr_test => true,:tt_arr_test => true)
k = BOAST::MFG(conv_filter, optims)
#k.print
k.build(:openmp => true)
bc[0] = BOAST::BC::PERIODIC
bc[1] = BOAST::BC::PERIODIC
bc[2] = BOAST::BC::PERIODIC

repeat = 5
begin
  stats_a = []
  repeat.times { 
    stats_a.push k.run(3, n, bc, ld_in, ld_out, m, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

m.times { |i|
  diff = (output_ref - output[0..(n1-1), 0..(n2-1), 0..(n3-1), i]).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{k.cost(n, bc, m) / (stats[:duration]*1.0e9)} GFlops"

puts 'again, grow'
input = NArray.float(n1,n2,n3,m).random
work1 = NArray.float(n1+15,n2+15,n3+15,m)
work2 = NArray.float(n1+15,n2+15,n3+15,m)
output_ref = NArray.float(n1+15,n2+15,n3+15)
output = NArray.float(n1+15,n2+15,n3+15,m)
(1..(m-1)).each { |i|
  input[ true, true, true, i ] = input[ true, true, true, 0 ]
}

k_ref = BOAST::magicfilter_per_ref(false,true)
stats = k_ref.run(n1, n2*n3, input, work2)
stats = k_ref.run(n2, (n1+15)*n3, work2, work1)
stats = k_ref.run(n3, (n2+15)*(n1+15), work1, output_ref)

puts "#{k_ref.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

bc[0] = BOAST::BC::GROW
bc[1] = BOAST::BC::GROW
bc[2] = BOAST::BC::GROW

ld_out[0] = n1+15
ld_out[1] = n2+15
ld_out[2] = n3+15

begin
  stats_a = []
  repeat.times { 
    stats_a.push k.run(3, n, bc, ld_in, ld_out, m, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

m.times { |i|
  diff = (output_ref - output[0..(n1+15-1), 0..(n2+15-1), 0..(n3+15-1), i]).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{k.cost(n, bc, m) / (stats[:duration]*1.0e9)} GFlops"


puts 'again, shrink'
input = NArray.float(n1+15,n2+15,n3+15,m).random
work1 = NArray.float(n1+15,n2+15,n3+15,m)
work2 = NArray.float(n1+15,n2+15,n3+15,m)
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3,m)
(1..(m-1)).each { |i|
  input[ true, true, true, i ] = input[ true, true, true, 0 ]
}

k_ref = BOAST::magicfilter_per_ref(true,true)
stats = k_ref.run(n1, (n2+15)*(n3+15), input, work2)
stats = k_ref.run(n2, n1*(n3+15), work2, work1)
stats = k_ref.run(n3, n2*n1, work1, output_ref)
puts "#{k_ref.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

conv_filter = BOAST::ConvolutionFilter::new('rfsf',FILTER.reverse,7)
k = BOAST::MFG(conv_filter, optims)
k.build(:openmp => true)

bc[0] = BOAST::BC::SHRINK
bc[1] = BOAST::BC::SHRINK
bc[2] = BOAST::BC::SHRINK

ld_in[0] = n1+15
ld_in[1] = n2+15
ld_in[2] = n3+15

ld_out[0] = n1
ld_out[1] = n2
ld_out[2] = n3

begin
  stats_a = []
  repeat.times { 
    stats_a.push k.run(3, n, bc, ld_in, ld_out, m, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

m.times { |i|
  diff = (output_ref - output[0..(n1-1), 0..(n2-1), 0..(n3-1), i]).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{k.cost(n, bc, m) / (stats[:duration]*1.0e9)} GFlops"

puts 'again, periodic inverse'
input = NArray.float(n1+pad_in,n2+pad_in,n3+pad_in,m).random
(1..(m-1)).each { |i|
  input[ true, true, true, i] = input[ true, true, true, 0]
}
work1 = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
work2 = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1+pad_out,n2+pad_out,n3+pad_out,m)
epsilon = 10e-15

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3

ld_in = NArray.int(3)
ld_in[0] = n1+pad_in
ld_in[1] = n2+pad_in
ld_in[2] = n3+pad_in

ld_out = NArray.int(3)
ld_out[0] = n1+pad_out
ld_out[1] = n2+pad_out
ld_out[2] = n3+pad_out

bc = NArray.int(3)

k_ref = BOAST::magicfilter_per_ref( true )
stats = k_ref.run(n1, n2*n3, input, output_ref)
stats = k_ref.run(n2, n1*n3, output_ref, work1)
stats = k_ref.run(n3, n2*n1, work1, output_ref)

puts "#{k_ref.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"

bc[0] = BOAST::BC::PERIODIC
bc[1] = BOAST::BC::PERIODIC
bc[2] = BOAST::BC::PERIODIC

repeat = 5
begin
  stats_a = []
  repeat.times { 
    stats_a.push k.run(3, n, bc, ld_in, ld_out, m, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

m.times { |i|
  diff = (output_ref - output[0..(n1-1), 0..(n2-1), 0..(n3-1), i]).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{k.cost(n, bc, m) / (stats[:duration]*1.0e9)} GFlops"


