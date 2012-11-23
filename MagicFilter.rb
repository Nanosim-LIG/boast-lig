require "./BOAST.rb"
module ConvolutionGenerator
  def ConvolutionGenerator::MagicFilter(filt, center, unroll, invert, free=false )
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
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

    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
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
    arr = ConstArray::new(filt)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil,upfil) ]})
  
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
      for2 = For::new(i,dim_out_min, dim_out_max) {
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
end

FILTER = [ "8.4334247333529341094733325815816e-7",
       "-0.1290557201342060969516786758559028e-4",
       "0.8762984476210559564689161894116397e-4",
       "-0.30158038132690463167163703826169879e-3",
       "0.174723713672993903449447812749852942e-2",
       "-0.942047030201080385922711540948195075e-2",
       "0.2373821463724942397566389712597274535e-1",
       "0.612625895831207982195380597e-1",
       "0.9940415697834003993178616713",
       "-0.604895289196983516002834636e-1",
       "-0.2103025160930381434955489412839065067e-1",
       "0.1337263414854794752733423467013220997e-1",
       "-0.344128144493493857280881509686821861e-2",
       "0.49443227688689919192282259476750972e-3",
       "-0.5185986881173432922848639136911487e-4",
       "2.72734492911979659657715313017228e-6" ]



ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
k = ConvolutionGenerator::MagicFilter(FILTER,8,0,false)
k.print
k.build
#ConvolutionGenerator::MagicFilter(FILTER,8,5,true).print
#ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true).print
#ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true).print
ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false).print
#ConvolutionGenerator::MagicFilter(FILTER,8,5,true).print
#ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true).print
k = ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true)
k.print
k.build

