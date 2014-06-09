require "BOAST"
require 'narray'
require "./GenericConvolution.rb" 
module BOAST
  def BOAST::magicfilter_per_ref( invert = false, free = false )
    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::FORTRAN)
    kernel = CKernel::new
    kernel.lang = BOAST::FORTRAN
    lowfil=-8
    lupfil=7
    function_name = "magicfilter"
    if not free then
      function_name += "_per"
    else
      function_name += "_free"
    end
    function_name += "_t" if invert
    function_name += "_ref"
    n = BOAST::Int("n",:dir => :in, :signed => false)
    ndat = BOAST::Int("ndat",:dir => :in, :signed => false)

    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    if free then
      if invert then
        dim_in_min = lowfil
        dim_in_max = n-1+lupfil
      else
        dim_out_min = -lupfil
        dim_out_max = n-1-lowfil
      end
    end

    x = BOAST::Real("x",:dir => :in, 
                    :dim => [ BOAST:Dim(dim_in_min, dim_in_max), BOAST:Dim(ndat) ])
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print(File::read("magicfilter_refs.f90"))
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end

  def BOAST::MagicFilter(filt, center, unroll, invert, free=false )
    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang
    function_name = "magicfilter"
    if free then
      function_name += "_free"
    else 
      function_name += "_per"
    end
    if invert then
      function_name += "_inv"
    end
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in})
    ndat = Variable::new("ndat",Int,{:direction => :in})
    if invert then
      lowfil = Variable::new("lowfil",Int,{:constant => 1-center})
      upfil = Variable::new("upfil",Int,{:constant => filt.length-center})
    else
      lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
      upfil = Variable::new("upfil",Int,{:constant => center-1})
    end

    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    if free then
      if invert then
        dim_in_min = lowfil
        dim_in_max = n-1+upfil
      else
        dim_out_min = -upfil
        dim_out_max = n-1-lowfil
      end
    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    tt = [Variable::new("tt1",Real)]
    2.upto(unroll) { |index|
      tt.push(Variable::new("tt#{index}",Real))
    }
    if invert then filt = filt.reverse end 
    arr = ConstArray::new(filt,Real)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil,upfil) ]})

    if BOAST::get_lang == C then
      @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      @@output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      @@output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end  
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil]) {
      i.decl
      j.decl
      k.decl
      l.decl
      tt.each { |t|
        t.decl
      }
  
      fil.decl
  
      if unroll>1 then 
        for1 = For::new(j,1,ndat-(unroll-1), unroll)
      else
        for1 = For::new(j,1,ndat)
      end
      for1.print
        for2 = For::new(i, dim_out_min, dim_out_max) {
          tt.each{ |t|
            (t === 0.0).print
          }
          if free and not invert
            for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) )
          else
            for3 = For::new(l, lowfil, upfil)
          end
          for3.print
            if not free then
              (k === FuncCall::new( "modulo", i+l, n)).print
            else
              (k === i+l).print
            end
            tt.each_index{ |index|
              (tt[index] === tt[index] + x[k,j+index]*fil[l]).print
            }
          for3.close
          tt.each_index{ |index|
            (y[j+index,i] === tt[index]).print
        }
      }.print
      for1.close
      
      if unroll>1 then
        for1 = For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          for2 = For::new(i,dim_out_min,dim_out_max) {
            (tt[0] === 0.0).print
            if free and not invert then
              for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) ) {
                (k === i+l).print
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.print
            else
              for3 = For::new(l, lowfil, upfil) {
                if not free then
                  (k === FuncCall::new( "modulo", i+l, n)).print
                else
                  (k === i+l).print
                end
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.print
            end
            (y[j,i] === tt[0]).print
          }.print
        }.print
      end
    }
    p.print
    kernel.procedure = p
    return kernel
  end

  def BOAST::MF3d(filt, center, unroll, free=[BOAST::BC::PERIODIC]*3)
    kernel = CKernel::new
    BOAST::set_output( kernel.code )
    kernel.lang = BOAST::get_lang

    if BOAST::get_lang == C then
      @@output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      @@output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      @@output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end
    
    conv_filter = ConvolutionFilter::new('magfilt',filt,center)
    
    conv_operation = ConvolutionOperator::new(conv_filter,3,free,:work => true)

    optim = ConvolutionOptimization::new(conv_operation,:use_mod => true,:unroll => unroll, :transpose => 1,:tt_arr => true)

    p, subops= conv_operation.procedure(optim)
    subops[0].print
    #subops.each{ | ops| ops.print}
    p.print
    
    kernel.procedure = p
    return kernel
  end
end
