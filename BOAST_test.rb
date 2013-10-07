require './BOAST.rb'
require 'rubygems'
require 'narray'
module ConvolutionGenerator
  def ConvolutionGenerator::vector_add
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "vector_add"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    a = Variable::new("a",Real,{:direction => :in, :dimension => [ Dimension::new(n)] })
    b = Variable::new("b",Real,{:direction => :in, :dimension => [ Dimension::new(n)] })
    c = Variable::new("c",Real,{:direction => :out, :dimension => [ Dimension::new(n)] })
    i = Variable::new("i",Int,{:signed => false})
    p = Procedure::new(function_name, [n,a,b,c]) {
      if(ConvolutionGenerator::get_lang == ConvolutionGenerator::CUDA) then
      else
        i.decl
        f = For::new(i,1,n) {
           (c[i] === a[i] + b[i]).print
        }
        f.print
      end
    }
    p.print
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
ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
puts "FORTRAN"
k = ConvolutionGenerator::vector_add
k.run(n,a,b,c_ref)
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
puts "C"
k = ConvolutionGenerator::vector_add
k.run(n,a,b,c)
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
ConvolutionGenerator::set_lang( ConvolutionGenerator::CL )
puts "CL"
k = ConvolutionGenerator::vector_add
k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
