require './Poisson.rb'

#for a common case, the 3 filters are the same in all directions. This is the middle line divided by hgrids(x)

n01 = 200
n02 = 200
n03 = 200
nord = 16

filter= NArray.float(nord+1, nord+1)


hgrids = NArray.float(3)
hgrids[0] = 0.1
hgrids[1] = 0.1
hgrids[2] = 0.1
geocode = 3
epsilon = 10e-11
$options = {}
$parser = OptionParser::new do |opts|
  opts.on("-o", "--output", "output code") {
    $options[:output] = true
  }
  opts.on("-f", "--fast", "use less optimization flags to run faster") {
    $options[:fast] = true
  }
  opts.on("-i", "--input-file FILE", "use input convolutions already optimized") { |file|
    $options[:input_file] = file
  }
  opts.parse!
end

n = NArray.int(3)
n[0] = n01
n[1] = n02
n[2] = n03

bc = NArray.int(3)

bc[0] = BC::PERIODIC
bc[1] = BC::PERIODIC
bc[2] = BC::PERIODIC
if(geocode == 1) then
bc[0] = BC::NPERIODIC
bc[2] = BC::NPERIODIC
end
if(geocode != 1) then
bc[1] = BC::NPERIODIC
end


generate_filter.run(nord, filter)

u = NArray.float(n01,n02,n03,3).random
du_ref = NArray.float(n01,n02,n03)
du = NArray.float(n01,n02,n03)
du_3D = NArray.float(n01,n02,n03)
du_boast = NArray.float(n01,n02,n03)
cc_ref = NArray.float(n01,n02,n03)

k = div_u_i_ref
stats = k.run(geocode, n01, n02, n03, u, du_ref, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

k = div_u_i_opt
stats = k.run(geocode, n01, n02, n03, u, du, nord, hgrids, filter,cc_ref)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

diff = (du_ref - du).abs
diff.each { |elem|
  raise "Warning: residue opt too big: #{elem}" if elem > epsilon
}



#3*1D version with boast convolutions : works

k = div_u_i_1d
stats = k.run(geocode, n01, n02, n03, u, du_3D,1, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
stats = k.run(geocode, n01, n02, n03, u, du_3D,2, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
stats = k.run(geocode, n01, n02, n03, u, du_3D,3, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"


diff = (du_ref - du_3D).abs
diff.each { |elem|
  raise "Warning: residue opt too big: #{elem}" if elem > epsilon
}




#single kernel that calls the 3 convolutions

#compute and optimize outside of the kernel, as we will reuse these convolutions later
conv_filter = PoissonFilter::new('poisson',filter.to_a,nord)


if $options[:input_file] then
  puts "Reading already optimized kernels from #{$options[:input_file]}"
  kconv = Poisson_brocker_from_file($options[:input_file],n01,n02,n03)
  puts "Compiling"
  kconv.build(:openmp => true)
  puts "Compilation done."
else

  if $options[:fast] then
  optims = GenericOptimization::new(:unroll_range => [5,5], :mod_arr_test => false,:tt_arr_test => false,:unrolled_dim_index_test => false, :unroll_inner_test =>false, :dimensions => [n01,n02,n03])
  else
  optims = GenericOptimization::new(:unroll_range => 5, :mod_arr_test => true,:tt_arr_test => true,:unrolled_dim_index_test => true, :unroll_inner_test =>true, :dimensions => [n01,n02,n03])
  end

  kconv = Poisson_brocker(optims,n01,n02,n03)
  kconv.build(:openmp => true)
end

cc = NArray.float(n01,n02,n03)

#kconv = Poisson_conv(conv_filter, optims)
#kconv.build(:openmp => true)

#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts kconv
#  }

k2 = div_u_i(n01,n02,n03,kconv)
k2.build(:openmp => true)


stats_a = []
stats_a.push k2.run(geocode, n01,n02,n03,u,du_boast,nord,hgrids,cc)

diff = (du_ref - du_boast).abs
diff.each { |elem|
    raise "Warning: residue too big: #{elem}" if elem > epsilon
}

diff = (cc_ref - cc).abs
diff.each { |elem|
    raise "Warning: residue too big for cc: #{elem}" if elem > epsilon
}
#div_u_i_boast(geocode, n01, n02, n03, u, du_3D3, nord, hgrids)
#stats = k2.run(geocode, n01,n02,n03,hgrids,u,u1,u2,du_boast)
repeat = 5
begin
  repeat.times { |i|
  stats_a.push k2.run(geocode, n01,n02,n03,u,du_boast,nord,hgrids,cc)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k2.procedure.name}: #{stats[:duration]*1.0e3} #{k2.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"

cc=nil
u = nil
du_ref = nil
du = nil
du_3D = nil
du_boast = nil


#test for the fssnord3dmatnabla3varde2_lg function

u = NArray.float(n01,n02,n03).random
du_ref = NArray.float(n01,n02,n03,3)
du2_ref = NArray.float(n01,n02,n03)
du_boast = NArray.float(n01,n02,n03,3)
du2_boast = NArray.float(n01,n02,n03)

k = nabla_u_and_square_ref
stats = k.run(geocode, n01, n02, n03, u, du_ref,du2_ref , nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

#annoyingly, too much computation to get a good precision here ...
epsilon = 10e-6

k3 = nabla_u_and_square(n01,n02,n03,kconv)

k3.build(:openmp => true)
#du0=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 0]
#du1=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 1]
#du2=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 2]

stats_a = []
stats_a.push k3.run(geocode, n01,n02,n03,u,du_boast,du2_boast,nord,hgrids)

diff = (du2_ref - du2_boast).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

repeat = 5
begin
  repeat.times { |i|
  stats_a.push k3.run(geocode, n01,n02,n03,u,du_boast,du2_boast, nord, hgrids)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k3.procedure.name}: #{stats[:duration]*1.0e3} #{k3.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"


u=nil
du_boast=nil
du2_boast=nil
du_ref=nil
du2_ref=nil

#now the fssnord3dmatnabla_lg variant


u = NArray.float(n01,n02,n03).random
dlogeps = NArray.float(3,n01,n02,n03).random
rhopol_ref = NArray.float(n01,n02,n03)
rhopol_boast = NArray.float(n01,n02,n03)
rhores2_ref = 0.0
rhores2_boast = 0.0
eta = 0.5
du_boast = NArray.float(n01,n02,n03,3)
#du2_boast = NArray.float(n01,n02,n03)

k = update_rhopol_ref
stats = k.run(geocode, n01, n02, n03, u, nord, hgrids,eta,dlogeps,rhopol_ref,rhores2_ref,  filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"


k4 = update_rhopol(n01,n02,n03,kconv)


k4.build(:openmp => true)
#du0=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 0]
#du1=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 1]
#du2=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 2]

stats_a = []
stats_a.push k4.run(geocode, n01,n02,n03,u, nord,hgrids,eta,dlogeps,rhopol_boast,rhores2_boast)

diff = (rhopol_ref - rhopol_boast).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}


repeat = 5
begin
  repeat.times { |i|
  stats_a.push k4.run(geocode, n01,n02,n03,u,nord, hgrids,eta,dlogeps,rhopol_boast,rhores2_boast)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k4.procedure.name}: #{stats[:duration]*1.0e3} #{k4.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"

u = nil
dlogeps = nil
rhopol_ref = nil
rhopol_boast = nil


#now the fssnord3dmatnabla3var_lg variant
epsilon = 10e-11

u = NArray.float(n01,n02,n03).random
du_boast = NArray.float(n01,n02,n03,3)
du_ref = NArray.float(n01,n02,n03,3)

k = nabla_u_ref
stats = k.run(geocode, n01, n02, n03, u, du_ref, nord, hgrids,filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"


k5 = nabla_u(n01,n02,n03,kconv)


k5.build(:openmp => true)

stats_a = []
stats_a.push k5.run(geocode, n01,n02,n03,u,du_boast,nord,hgrids)

diff = (du_ref - du_boast).abs
diff.each { |elem|
  raise "Warning: residue 0 too big: #{elem}" if elem > epsilon
}

repeat = 5
begin
  repeat.times { |i|
  stats_a.push k5.run(geocode, n01,n02,n03,u,du_boast,nord,hgrids)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k5.procedure.name}: #{stats[:duration]*1.0e3} #{k5.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"


u = nil
du_boast = nil
du_ref = nil



#now the nabla_u_epsilon variant


u = NArray.float(n01,n02,n03).random
du_boast = NArray.float(n01,n02,n03,3)
du_ref = NArray.float(n01,n02,n03,3)
du_3D = NArray.float(n01,n02,n03).random

k = nabla_u_epsilon_ref
stats = k.run(geocode, n01, n02, n03, u, du_ref,nord, hgrids,du_3D, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"


k6 = nabla_u_epsilon(n01,n02,n03,kconv)


k6.build(:openmp => true)

stats_a = []
stats_a.push k6.run(geocode, n01,n02,n03,u,du_boast, nord,hgrids,du_3D)

epsilon=10e-11
diff = (du_ref - du_boast).abs
diff.each { |elem|
  raise "Warning: residue 0 too big: #{elem}" if elem > epsilon
}

repeat = 5
begin
  repeat.times { |i|
  stats_a.push k6.run(geocode, n01,n02,n03,u,du_boast, nord,hgrids,du_3D)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k6.procedure.name}: #{stats[:duration]*1.0e3} #{k5.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"


u = nil
du_boast = nil
du_ref = nil
du_3D = nil

#test for the nabla_u_square function

u = NArray.float(n01,n02,n03).random
du2_ref = NArray.float(n01,n02,n03)
du_boast = NArray.float(n01,n02,n03,3)
du2_boast = NArray.float(n01,n02,n03)

k = nabla_u_square_ref
stats = k.run(geocode, n01, n02, n03, u, du2_ref , nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

#annoyingly, too much computation to get a good precision here ...
epsilon = 10e-6


#cheat to use and square ?
k7 = nabla_u_square(n01,n02,n03,kconv)

k7.build(:openmp => true)
#du0=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 0]
#du1=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 1]
#du2=du_boast[0..(n01-1), 0..(n02-1), 0..(n03-1), 2]

stats_a = []
stats_a.push k7.run(geocode, n01,n02,n03,u,du2_boast,nord,hgrids)

diff = (du2_ref - du2_boast).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

repeat = 5
begin
  repeat.times { |i|
  stats_a.push k7.run(geocode, n01,n02,n03,u,du2_boast,nord,hgrids)
}
rescue Exception => e
  puts e.inspect
end

stats_a.sort_by! { |a| a[:duration] }
stats = stats_a.first

puts "#{k7.procedure.name}: #{stats[:duration]*1.0e3} #{k3.cost(0,n, bc) / (stats[:duration]*1.0e9)} GFlops"


u=nil
du_boast=nil
du2_boast=nil
du2_ref=nil


if $options[:output] then
  suffix = ".c" if BOAST::get_lang == BOAST::C
  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
  File::open("poisson#{suffix}","w") { |f|
    f.puts k2
    f.puts k3
    f.puts k4
    f.puts k5
    f.puts k6
    f.puts k7
  }
end

#test for the fssnord3dmatnabla_lg function

