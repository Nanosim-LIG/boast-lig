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
$options = { :unroll => 6, :step => 1, :precision => 8 }
$parser = OptionParser::new do |opts|
  opts.on("-s", "--step VAL", "Unroll step") { |step|
    $options[:step] = step.to_i
  }
  opts.on("-d", "--dump", "Dump binary") { |step|
    $options[:dump] = true
  }
  opts.on("-u", "--unroll VAL", "Unroll bound") { |unroll|
    $options[:unroll] = unroll.to_i
  }
  opts.on("-f", "--force_unroll", "Force outer unroll") {
    $options[:force_unroll] = true
  }
  opts.on("-o", "--output", "output code") {
    $options[:output] = true
  }
  opts.on("-p", "--precision VAL", "Precision of the code") { |precision|
    $options[:precision] = precision.to_i
  }
  opts.parse!
end

case $options[:precision]
when 4
  type = NArray::SFLOAT
  epsilon = 10e-6
when 8
  type = NArray::FLOAT
  epsilon = 10e-12
else
  raise "Unsupported precision!"
end

BOAST::default_real_size = $options[:precision]

if $options[:force_unroll] then
  unroll_range = [$options[:unroll],$options[:unroll]]
else
  unroll_range = [1,$options[:unroll],$options[:step]]
end

optims = GenericOptimization::new(:unroll_range => unroll_range, :mod_arr_test => true, :tt_arr_test => true, :unroll_inner_test => true)
wave_filter = WaveletFilter::new("sym#{L.length/2}", L)
n1 = 62
n2 = 66
n3 = 65
input   = NArray::new(type, n1*2, n2*2, n3*2).random
work1   = NArray::new(type, n1*2, n2*2, n3*2)
work2   = NArray::new(type, n1*2, n2*2, n3*2)
output  = NArray::new(type, n1*2, n2*2, n3*2)
output2 = NArray::new(type, n1*2, n2*2, n3*2)

n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3
bc = NArray.int(3)
bc[0] = BC::PERIODIC
bc[1] = BC::PERIODIC
bc[2] = BC::PERIODIC


k1 = Wavelet(wave_filter, :decompose, optims)
k1.build(:openmp => true)
k2 = Wavelet(wave_filter, :recompose, optims)
k2.build(:openmp => true)

if $options[:output] then
  suffix = ".c" if BOAST::get_lang == BOAST::C
  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
  File::open("wavelets#{suffix}","w") { |f|
    f.puts k1
    f.puts k2
  }
end
if $options[:dump] then
  k1.dump_binary
  k2.dump_binary
end

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

input   = NArray::new(type, n1*2, n2*2, n3*2).random
work1   = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
work2   = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output  = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output2 = NArray::new(type, n1*2, n2*2, n3*2)

bc[0] = BC::GROW
bc[1] = BC::GROW
bc[2] = BC::GROW

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


bc[0] = BC::SHRINK
bc[1] = BC::SHRINK
bc[2] = BC::SHRINK

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

input   = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2).random
work1   = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
work2   = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)
output  = NArray::new(type, n1*2, n2*2, n3*2)
output2 = NArray::new(type, n1*2 + L.length - 2, n2*2 + L.length - 2, n3*2 + L.length - 2)


bc[0] = BC::SHRINK
bc[1] = BC::SHRINK
bc[2] = BC::SHRINK

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

bc[0] = BC::GROW
bc[1] = BC::GROW
bc[2] = BC::GROW
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
lowbound = L.length - 2
highbound = -L.length + 1
area = [lowbound..highbound]*3
diff = (input[*area] - output2[*area]).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

