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

#n1 = 124
#n2 = 132
#n3 = 130
#input = NArray.float(n1,n2,n3).random
#output_ref = NArray.float(n1,n2,n3)
#output = NArray.float(n1,n2,n3)
#epsilon = 10e-15
#ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
#k = ConvolutionGenerator::synthesis_per_ref
#stats = k.run(n1/2, n2*n3, input, output_ref)
#puts "#{k.procedure.name}: #{stats["duration"]*1.0e3} #{32*n1*n2*n3 / (stats["duration"]*1.0e9)} GFlops"
#
#(0..8).each{ |unroll|
#  k = ConvolutionGenerator::Synthesis(FILTER,7,unroll,false)
#  k.print if unroll == 0
#  k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
#  stats = k.run(n1/2, n2*n3, input, output)
#  stats = k.run(n1/2, n2*n3, input, output)
#  diff = (output_ref - output).abs
#  diff.each { |elem|
#    raise "Warning: residue too big: #{elem}" if elem > epsilon
#  }
#  puts "#{k.procedure.name}: #{stats["duration"]*1.0e3} #{32*n1*n2*n3 / (stats["duration"]*1.0e9)} GFlops"
#}
#ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
#(0..8).each{ |unroll|
#  k = ConvolutionGenerator::Synthesis(FILTER,7,unroll,false)
#
#  stats = k.run(n1/2, n2*n3, input, output)
#  stats = k.run(n1/2, n2*n3, input, output)
#  diff = (output_ref - output).abs
#  diff.each { |elem|
#    raise "Warning: residue too big: #{elem}" if elem > epsilon
#  }
#  puts "#{k.procedure.name}: #{stats["duration"]*1.0e3} #{32*n1*n2*n3 / (stats["duration"]*1.0e9)} GFlops"
#}

n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n2,n3,n1+14)
output = NArray.float(n2,n3,n1+14)
epsilon = 10e-15
ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
k = ConvolutionGenerator::synthesis_free_ref
k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
stats = k.run(n1/2, n2*n3, input, output_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
(0..8).each{ |unroll|
  k = ConvolutionGenerator::Synthesis(FILTER,7,unroll,true)
#  k.print if unroll == 0
#  k.build({:FC => 'gfortran-4.6',:CC => 'gcc',:FCFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  stats = k.run(n1/2, n2*n3, input, output)
  diff.each { |elem|
   raise  "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
(0..8).each{ |unroll|
  k = ConvolutionGenerator::Synthesis(FILTER,7,unroll,true)
#  k.print if unroll == 0
  k.build({:FC => 'fortran',:CC => 'gcc',:CFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
  stats = k.run(n1/2, n2*n3, input, output)
  stats = k.run(n1/2, n2*n3, input, output)
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{32*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
}
