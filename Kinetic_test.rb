require "./Kinetic.rb"

filt = ["-6.924474940639200152025730585882e-18",
        " 2.70800493626319438269856689037647576e-13",   # ---------------------
        "-5.813879830282540547959250667e-11",
        "-1.05857055496741470373494132287e-8",
        "-3.7230763047369275848791496973044e-7",
        " 2.0904234952920365957922889447361e-6",
        "-0.2398228524507599670405555359023135e-4",
        " 0.45167920287502235349480037639758496e-3",
        "-0.409765689342633823899327051188315485e-2",
        " 0.02207029188482255523789911295638968409e0",
        "-0.0822663999742123340987663521e0",
        " 0.2371780582153805636239247476e0",
        "-0.6156141465570069496314853949e0",
        " 2.2191465938911163898794546405e0",
        "-3.5536922899131901941296809374e0",
        " 2.2191465938911163898794546405e0",
        "-0.6156141465570069496314853949e0",
        " 0.2371780582153805636239247476e0",
        "-0.0822663999742123340987663521e0",
        " 0.02207029188482255523789911295638968409e0",
        "-0.409765689342633823899327051188315485e-2",
        " 0.45167920287502235349480037639758496e-3",
        "-0.2398228524507599670405555359023135e-4",
        " 2.0904234952920365957922889447361e-6",
        "-3.7230763047369275848791496973044e-7",
        "-1.05857055496741470373494132287e-8",
        "-5.813879830282540547959250667e-11",
        " 2.70800493626319438269856689037647576e-13",   # ---------------------
        "-6.924474940639200152025730585882e-18"]

n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1,n2,n3).random
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
hgrid = NArray.float(3)
hgrid[0] = 0.5
hgrid[1] = 0.6
hgrid[2] = 0.7
kstrten_ref = NArray.float(3)
kstrten = NArray.float(3)
epsilon = 10e-15
ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
k = ConvolutionGenerator::kinetic_per_ref
stats = k.run(n1, n2, n3, hgrid, input, output_ref, 0.5)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
k = ConvolutionGenerator::kinetic_per_ref_optim
k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
stats = k.run(n1-1, n2-1, n3-1, hgrid, input, output, 0.5)
stats = k.run(n1-1, n2-1, n3-1, hgrid, input, output, 0.5)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
diff = (output_ref - output).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

(1..14).each{ |unroll|
  k = ConvolutionGenerator::kinetic(filt,14,unroll)
  #k.print
  #k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fbounds-check",:LDFLAGS => "-lgfortran"})
  k.build({:FC => 'gfortran',:CC => 'gcc',:FCFLAGS => "-O2 -fopenmp",:LDFLAGS => "-fopenmp"})
  begin
    stats = k.run(n1, n2, n3, hgrid, input, output, 0.5)
    stats = k.run(n1, n2, n3, hgrid, input, output, 0.5)
  rescue Exception => e
    puts e.inspect
  end
  puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
  diff = (output_ref - output).abs
  diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
  }
}
k = ConvolutionGenerator::kinetic_per_ref_optim_ekin
stats = k.run(n1, n2, n3, hgrid, input, output_ref, kstrten_ref)
k = ConvolutionGenerator::kinetic(filt,14,1,true)
k.print
stats = k.run(n1, n2, n3, hgrid, input, output, kstrten)
diff = (output_ref - output).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
diff = (kstrten_ref - kstrten)
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
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
#
