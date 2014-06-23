require './Wavelet.rb'

#L =[-0.12940952255092145,
#    0.22414386804185735,
#    0.836516303737469,
#    0.48296291314469025]


#L= [0.0014009155259146807,
#    0.0006197808889855868,
#    -0.013271967781817119,
#    -0.01152821020767923,
#    0.03022487885827568,
#    0.0005834627461258068,
#    -0.05456895843083407,
#    0.238760914607303,
#    0.717897082764412,
#    0.6173384491409358,
#    0.035272488035271894,
#    -0.19155083129728512,
#    -0.018233770779395985,
#    0.06207778930288603,
#    0.008859267493400484,
#    -0.010264064027633142,
#    -0.0004731544986800831,
#    0.0010694900329086053 ]

L = [-0.0033824159510050025955,
     -0.00054213233180001068935,
      0.031695087811525991431,
      0.0076074873249766081919,
     -0.14329423835127266284,
     -0.061273359067811077843,
      0.48135965125905339159,
      0.77718575169962802862,
      0.36444189483617893676,
     -0.051945838107881800736,
     -0.027219029917103486322,
      0.049137179673730286787,
      0.0038087520138944894631,
     -0.014952258337062199118,
     -0.00030292051472413308126,
      0.0018899503327676891843]


optims = BOAST::GenericOptimization::new(:unroll_range => 6, :mod_arr_test => true, :tt_arr_test => true)
wave_filter = BOAST::WaveletFilter::new("sym#{L.length/2}", L)
n1 = 62
n2 = 66
n3 = 65
input   = NArray.float(n1*2, n2*2, n3*2).random
work1   = NArray.float(n1*2, n2*2, n3*2)
work2   = NArray.float(n1*2, n2*2, n3*2)
output  = NArray.float(n1*2, n2*2, n3*2)
output2 = NArray.float(n1*2, n2*2, n3*2)

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3
bc = NArray.int(3)
bc[0] = BOAST::BC::PERIODIC
bc[1] = BOAST::BC::PERIODIC
bc[2] = BOAST::BC::PERIODIC

epsilon = 10e-12

k1 = BOAST::Wavelet(wave_filter, :decompose, optims)
k1.build(:openmp => true)
k2 = BOAST::Wavelet(wave_filter, :recompose, optims)
k2.build(:openmp => true)

repeat = 5
begin
  stats_a = []
  repeat.times {
    stats_a.push k1.run(3, n, bc, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k1.procedure.name}: #{stats[:duration]*1.0e3} #{k1.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"

begin
  stats_a = []
  repeat.times {
    stats_a.push k2.run(3, n, bc, output, output2, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k2.procedure.name}: #{stats[:duration]*1.0e3} #{k2.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"

diff = (input - output2).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

puts "again grow then shrink"

input   = NArray.float(n1*2, n2*2, n3*2).random
work1   = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
work2   = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output  = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output2 = NArray.float(n1*2, n2*2, n3*2)

bc[0] = BOAST::BC::GROW
bc[1] = BOAST::BC::GROW
bc[2] = BOAST::BC::GROW

begin
  stats_a = []
  repeat.times {
    stats_a.push k1.run(3, n, bc, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k1.procedure.name}: #{stats[:duration]*1.0e3} #{k1.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"


bc[0] = BOAST::BC::SHRINK
bc[1] = BOAST::BC::SHRINK
bc[2] = BOAST::BC::SHRINK

begin
  stats_a = []
  repeat.times {
    stats_a.push k2.run(3, n, bc, output, output2, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k2.procedure.name}: #{stats[:duration]*1.0e3} #{k2.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"

diff = (input - output2).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

puts "again shrink then grow"

input   = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2).random
work1   = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
work2   = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output  = NArray.float(n1*2, n2*2, n3*2)
output2 = NArray.float(n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)


bc[0] = BOAST::BC::SHRINK
bc[1] = BOAST::BC::SHRINK
bc[2] = BOAST::BC::SHRINK

begin
  stats_a = []
  repeat.times {
    stats_a.push k1.run(3, n, bc, input, output, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k1.procedure.name}: #{stats[:duration]*1.0e3} #{k1.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"

bc[0] = BOAST::BC::GROW
bc[1] = BOAST::BC::GROW
bc[2] = BOAST::BC::GROW
n[0] = n1
n[1] = n2
n[2] = n3



repeat = 5
begin
  stats_a = []
  repeat.times {
    stats_a.push k2.run(3, n, bc, output, output2, work1, work2)
  }
rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first
puts "#{k2.procedure.name}: #{stats[:duration]*1.0e3} #{k2.cost(n, bc) / (stats[:duration]*1.0e9)} GFlops"
lowbound = L.length
highbound = -L.length
area = [lowbound..highbound]*3
diff = (input[*area] - output2[*area]).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

