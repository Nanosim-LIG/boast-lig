require "./Kinetic-2.rb"

FILTER= ["-6.924474940639200152025730585882e-18",
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

conv_filter = BOAST::ConvolutionFilter::new('kin',FILTER,14)

n1 = 124
n2 = 132
n3 = 130
input = NArray.float(n1,n2,n3).random
input_y = NArray.float(n1,n2,n3)
input_y = 0.5*input
output_y = NArray.float(n1,n2,n3)
output_x = NArray.float(n1,n2,n3)
work1 = NArray.float(n1,n2,n3)
work2 = NArray.float(n1,n2,n3)
output_ref = NArray.float(n1,n2,n3)
output = NArray.float(n1,n2,n3)
hgrid = NArray.float(3)
hgrid[0] = 0.5
hgrid[1] = 0.6
hgrid[2] = 0.7
scal=NArray.float(3)
(0..2).each{ |ind| scal[ind] = -0.5 / (hgrid[ind]*hgrid[ind])   }
kstrten_ref = NArray.float(3)
kstrten = NArray.float(3)
epsilon = 10e-13
k = BOAST::kinetic_per_ref
stats = k.run(n1, n2, n3, hgrid, input, output_ref, 0.5)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
k = BOAST::kinetic_per_ref_optim
k.build(:openmp => true)
#k.build({:FC => 'ifort',:CC => 'icc',:FCFLAGS => "-O2 -axSSE4.2 -openmp",:LDFLAGS => "-openmp"})
stats = k.run(n1-1, n2-1, n3-1, hgrid, input, output, 0.5)
stats = k.run(n1-1, n2-1, n3-1, hgrid, input, output, 0.5)
stats = k.run(n1-1, n2-1, n3-1, hgrid, input, output, 0.5)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
diff = (output_ref - output).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
output = NArray.float(n1,n2,n3).random
n = NArray.int(3)
n[0] = n1
n[1] = n2
n[2] = n3

k = BOAST::kineticG(conv_filter)
#k.print
k.build(:openmp => true)

bc = NArray.int(3)
bc[0] = BOAST::BC::PERIODIC
bc[1] = BOAST::BC::PERIODIC
bc[2] = BOAST::BC::PERIODIC

repeat = 5

begin
  stats_a = []
  repeat.times { 
    stats_a.push  k.run(3, n, bc, input, output, scal, 0.5)
  }

rescue Exception => e
  puts e.inspect
end
stats_a.sort_by! { |a| a[:duration] }
  
stats = stats_a.first
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n1*n2*n3 / (stats[:duration]*1.0e9)} GFlops"
diff = (output_ref - output).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
