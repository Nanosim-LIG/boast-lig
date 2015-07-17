require './Poisson.rb'

#for a common case, the 3 filters are the same in all directions. This is the middle line divided by hgrids(x)

n01 = 200
n02 = 200
n03 = 200
nord = 16

filter= NArray.float(nord, nord)

u = NArray.float(n01,n02,n03,3).random
du_ref = NArray.float(n01,n02,n03)
du = NArray.float(n01,n02,n03)
du_3D = NArray.float(n01,n02,n03)
du_3D2 = NArray.float(n01,n02,n03)
du_3D3 = NArray.float(n01,n02,n03)
du_3D4 = NArray.float(n01,n02,n03)
hgrids = NArray.float(3)
hgrids[0] = 0.1
hgrids[1] = 0.1
hgrids[2] = 0.1
geocode = 2
epsilon = 10e-13
$options = { }
$parser = OptionParser::new do |opts|
  opts.on("-o", "--output", "output code") {
    $options[:output] = true
  }
  opts.parse!
end


generate_filter.run(nord, filter)

k = fssnord3dmatnabla3var_lg_ref
stats = k.run(geocode, n01, n02, n03, u, du_ref, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

k = fssnord3dmatnabla3var_lg_opt
stats = k.run(geocode, n01, n02, n03, u, du, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

diff = (du_ref - du).abs
diff.each { |elem|
  raise "Warning: residue opt too big: #{elem}" if elem > epsilon
}


k = fssnord3dmatnabla3var_lg_1d
stats = k.run(geocode, n01, n02, n03, u, du_3D,1, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
stats = k.run(geocode, n01, n02, n03, u, du_3D,2, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
stats = k.run(geocode, n01, n02, n03, u, du_3D,3, nord, hgrids, filter)
puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"


diff = (du_ref - du_3D).abs
diff.each { |elem|
  raise "Warning: residue 3d too big: #{elem}" if elem > epsilon
}


#k = fssnord3dmatnabla3var_lg_3d(fssnord3dmatnabla3var_lg_1d)
#stats = k.run(geocode, n01, n02, n03, u, du_3D4,1, nord, hgrids)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

#fails with `block in get_sub_kernels': undefined method `rewind' for nil:NilClass


#3*1D version with boast convolutions : work

n = NArray.int(3)
n[0] = n01
n[1] = n02
n[2] = n03
FILTER_x=FILTER.collect { |n| n / hgrids[0] }
conv_filter = ConvolutionFilter::new('poisson',FILTER_x,8)
optims = GenericOptimization::new(:unroll_range => 1, :mod_arr_test => true,:tt_arr_test => true, :dimensions => [n01,n02,n03])

k2 = Poisson_conv(conv_filter, optims)
k2.build(:openmp => true)
bc = NArray.int(3)
bc[0] = BC::PERIODIC
bc[1] = BC::PERIODIC
bc[2] = BC::PERIODIC

stats1 = k2.run(3, 0, n, BC::PERIODIC, u, du_3D)
stats2 = k2.run(3, 1, n, BC::PERIODIC, u[0..(n01-1), 0..(n02-1), 0..(n03-1), 1], du_3D2)
stats3 = k2.run(3, 2, n, BC::PERIODIC, u[0..(n01-1), 0..(n02-1), 0..(n03-1), 2], du_3D3)

du_3D3= du_3D+du_3D2+du_3D3


diff = (du_ref - du_3D3).abs
diff.each { |elem|
  raise "Warning: residue 3*1D too big: #{elem}" if elem > epsilon
}

#attempt with a single kernel that calls the 3 convolutions -> fails


k2 = fssnord3dmatnabla3var_lg_full(n01,n02,n03,hgrids)
k2.build(:openmp => true)
#fssnord3dmatnabla3var_lg_full(geocode, n01, n02, n03, u, du_3D3, nord, hgrids)
k2.run(n01,n02,n03,u,du_3D4)
#stats = k.run(geocode, n01, n02, n03, u, du_3D3, nord, hgrids)

#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"

#stats_a.sort_by! { |a| a[:duration] }
#stats = stats_a.first

#m.times { |i|
#  diff = (du_ref - du_3D2[0..(n01-1), 0..(n02-1), 0..(n03-1), i]).abs
#  diff.each { |elem|
#    raise "Warning: residue too big: #{elem}" if elem > epsilon
#  }
#}
#puts "{k2.procedure.name}: #{stats1[:duration]*1.0e3} #{k2.cost(0,n, bc) / (stats1[:duration]*1.0e9)} GFlops"
#puts "{k2.procedure.name}: #{stats2[:duration]*1.0e3} #{k2.cost(1,n, bc) / (stats2[:duration]*1.0e9)} GFlops"
#puts "{k2.procedure.name}: #{stats3[:duration]*1.0e3} #{k2.cost(2,n, bc) / (stats3[:duration]*1.0e9)} GFlops"

##stats = k2.run(geocode, n01, n02, n03, u, du_3D2, nord, hgrids)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
##let's keep these 2 for now'
#stats = k.run(geocode, n01, n02, n03, u, du_3D2,2, nord, hgrids)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
#stats = k.run(geocode, n01, n02, n03, u, du_3D2,3, nord, hgrids)
#puts "#{k.procedure.name}: #{stats[:duration]*1.0e3} #{3*59*n01*n02*n03 / (stats[:duration]*1.0e9)} GFlops"
n=0
diff = (du_ref - du_3D3).abs
diff.each { |elem|
  n=n+1 if elem > epsilon
}
puts "#{n} errors"


if $options[:output] then
  suffix = ".c" if BOAST::get_lang == BOAST::C
  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
  File::open("poisson#{suffix}","w") { |f|
    f.puts k
  }
end
#diff = (du_ref - du).abs
#diff.each { |elem|
#  raise "Warning: residue too big: #{elem}" if elem > epsilon
#}


#repeat = 1
#begin
#  stats_a = []
#  repeat.times {
#    stats_a.push k.run(geocode, n01, n02, n03, u, du, nord, hgrids)
#  }
#rescue Exception => e
#  puts e.inspect
#end

