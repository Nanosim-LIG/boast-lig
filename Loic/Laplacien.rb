require 'BOAST'
require 'narray_ffi'
include BOAST

set_array_start(0)
set_default_real_size(4)
set_default_align(32)

def laplace_kernel_ref(opts = {})
  default = { :laplace_order => 1 }
  default.update(opts)
  laplace_order = default[:laplace_order]
  n1 = Int("n1", :dir => :in)
  n2 = Int("n2", :dir => :in)
  n3 = Int("n3", :dir => :in)
  a = Real("a", :restrict => true, :dir => :in, :dim => [Dim(n1), Dim(n2), Dim(n3)])
  b = Real("b", :restrict => true, :dir => :out, :dim => [Dim(n1), Dim(n2), Dim(n3)])


  p = Procedure("laplace_#{laplace_order}", [n1, n2, n3, a, b]) {
    i = Int(:i)
    j = Int(:j)
    k = Int(:k)
    decl i,j,k
    pr For(i, laplace_order, n1 -1 - laplace_order) {
      pr For(j, laplace_order, n2 -1 - laplace_order) {
        pr For(k, laplace_order, n3 -1 - laplace_order) {
          pr b[i,j,k] === a[i,j,k] * 6.0
          (-laplace_order..laplace_order).each { |l|
            pr b[i,j,k] === b[i,j,k] - a[i+l,j,k] unless l == 0
            pr b[i,j,k] === b[i,j,k] - a[i,j+l,k] unless l == 0
            pr b[i,j,k] === b[i,j,k] - a[i,j,k+l] unless l == 0
          }
        }
      }
    }
  }
  return p.ckernel
 
end
def laplace_kernel(opts = {})
  default = { :laplace_order => 1, :unroll => 1 }
  default.update(opts)
  laplace_order = default[:laplace_order]
  unroll = default[:unroll]
  n1 = Int("n1", :dir => :in)
  n2 = Int("n2", :dir => :in)
  n3 = Int("n3", :dir => :in)
  a = Real("a", :restrict => true, :dir => :in, :dim => [Dim(n1), Dim(n2), Dim(n3)])
  b = Real("b", :restrict => true, :dir => :out, :dim => [Dim(n1), Dim(n2), Dim(n3)])


  p = Procedure("laplace_#{laplace_order}", [n1, n2, n3, a, b]) {
    i = Int(:i)
    j = Int(:j)
    k = Int(:k)
    decl i,j,k
    inner_block = lambda { |o|
      pr b[i+o,j,k] === a[i+o,j,k] * 6.0
      (-laplace_order..laplace_order).each { |l|
        pr b[i+o,j,k] === b[i+o,j,k] - a[i+o+l,j,k] unless l == 0
        pr b[i+o,j,k] === b[i+o,j,k] - a[i+o,j+l,k] unless l == 0
        pr b[i+o,j,k] === b[i+o,j,k] - a[i+o,j,k+l] unless l == 0
      }
    }

    pr For(i, laplace_order, n1 -1 - laplace_order - (unroll - 1), :step => unroll) {
      pr For(j, laplace_order, n2 -1 - laplace_order) {
        pr For(k, laplace_order, n3 -1 - laplace_order) {
          unroll.times(&inner_block) 
        }
      }
    }
    pr For(i, i+(unroll-1), n1 -1 - laplace_order) {
      pr For(j, laplace_order, n2 -1 - laplace_order) {
        pr For(k, laplace_order, n3 -1 - laplace_order) {
          1.times(&inner_block)
        }
      }
    }
  }
  return p.ckernel
 
end


dim = 128
a = ANArray.sfloat(32, dim, dim, dim).random!
b = ANArray.sfloat(32, dim, dim, dim)
b_ref = ANArray.sfloat(32, dim, dim, dim)
p_ref = laplace_kernel_ref
#(1...dim-1).each { |i|
#  (1...dim-1).each { |j|
#    (1...dim-1).each { |k|
#      b_ref[i,j,k] = a[i,j,k]*6.0
#      (-1..1).each  { |l|
#        b_ref[i,j,k] = b_ref[i,j,k] - a[i+l,j,k] unless l == 0
#        b_ref[i,j,k] = b_ref[i,j,k] - a[i,j+l,k] unless l == 0
#        b_ref[i,j,k] = b_ref[i,j,k] - a[i,j,k+l] unless l == 0
#      }
#    }
#  }
#}
p_ref.run(dim, dim, dim, a, b_ref)
epsilon = 10e-7
set_lang(C)
opt_space = OptimizationSpace::new(:unroll => 1..24)
optimizer = BruteForceOptimizer::new( opt_space )

puts optimizer.optimize { |opts|
  laplace_order = 1
  p = laplace_kernel(opts)
  puts opts
  p.run(dim, dim, dim, a, b)
  stats = 7.times.collect { p.run(dim, dim, dim, a, b) }
  stats.sort! { |s1, s2| s1[:duration] <=> s2[:duration] }
  diff = b - b_ref
  puts stats.first
  puts "Perf: #{7*dim*dim*dim/(stats.first[:duration]*1000000000.0)}"
  diff[laplace_order..-1-laplace_order,laplace_order..-1-laplace_order,laplace_order..-1-laplace_order] .each { |d|
    raise "Error #{d}!" if d > epsilon
  }
  stats.first[:duration]
}
