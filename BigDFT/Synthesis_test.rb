require './Synthesis.rb'

FILTER = ["0.0018899503327676891843",
          "-0.00030292051472413308126",
          "-0.014952258337062199118",
          "0.0038087520138944894631",
          "0.049137179673730286787",
          "-0.027219029917103486322",
          "-0.051945838107881800736",
          "0.36444189483617893676",
          "0.77718575169962802862",
          "0.48135965125905339159",
          "-0.061273359067811077843",
          "-0.14329423835127266284",
          "0.0076074873249766081919",
          "0.031695087811525991431",
          "-0.00054213233180001068935",
          "-0.0033824159510050025955"]

puts "Periodic Boundary Conditions"
n1 = 124
n2 = 132
n3 = 130
repetition = 10
out_fortran = File::new("synthesys_kernels.f90","w+");
out_c = File::new("synthesys_kernels.c","w+");

input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
input3D = NArray.float(n1,n2,n3).random
output3D_ref = NArray.float(n1,n2,n3)
output3D = NArray.float(n1,n2,n3)
temp3D = NArray.float(n1,n2,n3)
epsilon = 10e-15
BOAST::set_lang( BOAST::FORTRAN )
k = BOAST::synthesis_per_ref
out_fortran.puts k
k.build( :verbose => true )
stats = k.run(n1/2, n2*n3, input, output_ref)
puts "Reference"
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
puts "Reference 3D"
k3D = BOAST::synthesis3D_per(k)
out_fortran.puts k3D
#k3D.build(:verbose => true)
stats = k3D.run(n1,n2,n3, input3D, output3D, temp3D)
puts "#{k3D.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3*3 / (stats[:duration]*1.0e9)} GFlops"
puts "FORTRAN"
(0..8).each{ |unroll|
  k = BOAST::synthesis(FILTER,7,unroll,false)
  out_fortran.puts k
  k.build
  stats = k.run(n1/2, n2*n3, input, output)
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
puts "C"
BOAST::set_lang( BOAST::C )
(0..8).each{ |unroll|
  k = BOAST::synthesis(FILTER,7,unroll,false)
  out_c.puts k
  k.build
  stats = k.run(n1/2, n2*n3, input, output)
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
puts "Free Boundary Conditions"
n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n2,n3,n1+14)
output = NArray.float(n2,n3,n1+14)
epsilon = 10e-15
BOAST::set_lang( BOAST::FORTRAN )
k = BOAST::synthesis_free_ref
out_fortran.puts k
k.build(:openmp => true)
durations = []
repetition.times {
  stats = k.run(n1/2, n2*n3, input, output_ref)
  durations.push( stats[:duration] )
}
duration = durations.sort.first
puts "Reference"
puts "#{k.procedure.name}: #{duration*1.0e3} #{32*n1*n2*n3 / (duration*1.0e9)} GFlops"
puts "FORTRAN OpenMP"
(1..8).each{ |unroll|
  durations = []
  k = BOAST::synthesis(FILTER,7,unroll,true)
  out_fortran.puts k
  k.build( :openmp => true )
  repetition.times {
    stats = k.run(n1/2, n2*n3, input, output)
    durations.push( stats[:duration] )
  }
  diff = (output_ref - output).abs
  diff.each { |elem|
   raise  "Warning: residue too big: #{elem}" if elem > epsilon
  }
  duration = durations.sort.first
  puts "#{k.procedure.name}: #{duration*1.0e3} #{32*n1*n2*n3 / (duration*1.0e9)} GFlops"
}
puts "C OpenMP"
BOAST::set_lang( BOAST::C )
(1..8).each{ |unroll|
  durations = []
  k = BOAST::synthesis(FILTER,7,unroll,true)
  out_c.puts k
#  k.print if unroll == 0
  k.build( :openmp => true )
  repetition.times {
    stats = k.run(n1/2, n2*n3, input, output)
    durations.push( stats[:duration] )
  }
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  duration = durations.sort.first
  puts "#{k.procedure.name}: #{duration*1.0e3} #{32*n1*n2*n3 / (duration*1.0e9)} GFlops"
}
puts "FORTRAN"
BOAST::set_lang( BOAST::FORTRAN )
(1..8).each{ |unroll|
  durations = []
  k = BOAST::synthesis(FILTER,7,unroll,true)
  k.build
  repetition.times {
    stats = k.run(n1/2, n2*n3, input, output)
    durations.push( stats[:duration] )
  }
  diff = (output_ref - output).abs
  diff.each { |elem|
   raise  "Warning: residue too big: #{elem}" if elem > epsilon
  }
  duration = durations.sort.first
  puts "#{k.procedure.name}: #{duration*1.0e3} #{32*n1*n2*n3 / (duration*1.0e9)} GFlops"
}
puts "C"
BOAST::set_lang( BOAST::C )
(1..8).each{ |unroll|
  durations = []
  k = BOAST::synthesis(FILTER,7,unroll,true)
  k.build
  repetition.times {
    stats = k.run(n1/2, n2*n3, input, output)
    durations.push( stats[:duration] )
  }
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  duration = durations.sort.first
  puts "#{k.procedure.name}: #{duration*1.0e3} #{32*n1*n2*n3 / (duration*1.0e9)} GFlops"
}
