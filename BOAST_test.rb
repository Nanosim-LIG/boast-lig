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
    a = Variable::new("a",Real,{:direction => :in, :dimension => [ Dimension::new(0,n-1)] })
    b = Variable::new("b",Real,{:direction => :in, :dimension => [ Dimension::new(0,n-1)] })
    c = Variable::new("c",Real,{:direction => :out, :dimension => [ Dimension::new(0,n-1)] })
    i = Variable::new("i",Int,{:signed => false})
    ig = Variable::new("ig", Sizet)
    if kernel.lang == ConvolutionGenerator::CL and ConvolutionGenerator::get_default_real_size == 8 then
      $output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    end
    p = Procedure::new(function_name, [n,a,b,c]) {
      if(ConvolutionGenerator::get_lang == ConvolutionGenerator::CUDA) then
      elsif(ConvolutionGenerator::get_lang == ConvolutionGenerator::CL) then
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

ConvolutionGenerator::set_default_real_size(8)
epsilon = 10e-15

ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
puts "FORTRAN"
k = ConvolutionGenerator::vector_add
puts k.print
k.run(n,a,b,c_ref)
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
puts "C"
k = ConvolutionGenerator::vector_add
puts k.print
k.run(n,a,b,c)
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
ConvolutionGenerator::set_lang( ConvolutionGenerator::CL )
puts "CL"
k = ConvolutionGenerator::vector_add
puts k.print
k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}

ConvolutionGenerator::set_default_real_size(4)
epsilon = 10e-5

ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
puts "FORTRAN"
k = ConvolutionGenerator::vector_add
k.run(n,a,b,c_ref)
puts k.print
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
puts "C"
k = ConvolutionGenerator::vector_add
puts k.print
k.run(n,a,b,c)
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}
ConvolutionGenerator::set_lang( ConvolutionGenerator::CL )
puts "CL"
k = ConvolutionGenerator::vector_add
puts k.print
k.run(n, a, b, c, :global_work_size => [rndup(n,32), 1,1], :local_work_size => [32,1,1] )
diff = (c_ref - c).abs
diff.each { |elem|
  raise "Warning: residue too big: #{elem}" if elem > epsilon
}


