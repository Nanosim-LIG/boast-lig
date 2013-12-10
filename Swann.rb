require './BOAST.rb'
require 'rubygems'
require 'narray'

module BOAST
  def BOAST::swann_kernel(sec, var_idx)
    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang
    function_name = "swann_#{sec}_#{var_idx}"
    ncolidx = Variable::new("ncolidx",Int, :constant => (150000*16*16+150000*17))
    colidx = Variable::new("colidx",Int,{:direction => :in, :dimension => [ Dimension::new(0,ncolidx-1)] })
    a = Variable::new("a",Real,{:direction => :in, :dimension => [ Dimension::new(0,ncolidx-1)] })
    np = Variable::new("np",Int, :constant => 150003)
    p = Variable::new("p",Real,{:direction => :in, :dimension => [ Dimension::new(0,np)] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    sum = Variable::new("sum",Real)
    vars = [colidx,a,p]
    sw = Procedure::new(function_name, [colidx,a,p],[ncolidx,np]) {
      i.decl
      j.decl
      sum.decl
      @@output.puts "#pragma statement cache_subsector_size #{sec} #{12-sec}"
      @@output.puts "#pragma statement cache_subsector_assign #{vars[var_idx]}"
      For::new(i,1,np) {
        (sum === 0.0).print
        For::new(j,1,np) {
          (sum === sum + a[j] * p[colidx[j]]).print
        }.print
      }.print
    }
    sw.print
    kernel.procedure = sw
    return kernel
  end
end

BOAST::set_lang(BOAST::C)
(1..11).each { |sec|
  (0..2).each { |var_idx|
     k = BOAST::swann_kernel(sec,var_idx)
     k.print
  }
}
