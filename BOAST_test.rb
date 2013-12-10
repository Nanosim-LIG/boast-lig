require './BOAST.rb'
require 'rubygems'
require 'narray'
module BOAST
  def BOAST::vector_add
    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang
    function_name = "vector_add"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    a = Variable::new("a",Real,{:direction => :in, :dimension => [ Dimension::new(0,n-1)] })
    b = Variable::new("b",Real,{:direction => :in, :dimension => [ Dimension::new(0,n-1)] })
    c = Variable::new("c",Real,{:direction => :out, :dimension => [ Dimension::new(0,n-1)] })
    i = Variable::new("i",Int,{:signed => false})
    ig = Variable::new("ig", Sizet)
    if kernel.lang == BOAST::CL then
      @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    end
    if(BOAST::get_lang == BOAST::CUDA) then
      p = Procedure::new(function_name, [n,a,b,c])
      @@output.print <<EOF
__global__ void vector_add(unsigned int n, const double *a, const double *b, double *c) {
   unsigned int ig = blockDim.x * blockIdx.x + threadIdx.x;
   c[ig] = a[ig] + b[ig];
}
EOF
    else
      p = Procedure::new(function_name, [n,a,b,c]) {
        if(BOAST::get_lang == BOAST::CL) then
          ig.decl
          (ig === FuncCall::new( "get_global_id", 0)).print
          (c[ig] === a[ig] + b[ig]).print
        else
          i.decl
          f = For::new(i,0,n-1) {
            (c[i] === a[i] + b[i]).print
          }
          f.print
        end
      }
      p.print
    end
    kernel.procedure = p
    return kernel
  end
end

def rndup( val, div)
  return (val%div) == 0 ? val : val + div - (val%div)
end


n = 1024*1024
a = NArray.float(n).random
b = NArray.float(n).random
c = NArray.float(n)
c_ref = NArray.float(n)

epsilon = 10e-15

BOAST::set_lang( BOAST::FORTRAN )
puts "FORTRAN"
k = BOAST::vector_add
puts k.print
k.run(n,a,b,c_ref)
BOAST::set_lang( BOAST::C )
puts "C"
c.random
k = BOAST::vector_add
puts k.print
k.run(n,a,b,c)
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
BOAST::set_lang( BOAST::CL )
puts "CL"
c.random
k = BOAST::vector_add
puts k.print
k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
BOAST::set_lang( BOAST::CUDA )
puts "CUDA"
c.random
k = BOAST::vector_add
puts k.print
k.build(:LDFLAGS => " -L/usr/local/cuda-5.5.22/lib64")
k.run(n, a, b, c, :block_number => [rndup(n,32)/32, 1,1], :block_size => [32,1,1] )
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
