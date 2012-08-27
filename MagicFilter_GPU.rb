require './Algorithm.rb'

module ConvolutionGenerator
  def ConvolutionGenerator::MagicFilter_GPU(filt, center, size_n, max_work_item=256, local_mem_size=16384 )
    array_start = $array_start
    $array_start = 0

    function_name = "magicfilter_n#{size_n}"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
    upfil = Variable::new("upfil",Int,{:constant => center-1})
    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(ndat), Dimension::new(dim_in_min, dim_in_max) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(dim_out_min, dim_out_max), Dimension::new(ndat) ] })
    arr = ConstArray::new(filt)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(filt.length) ]})
    
    wgs = 8
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil], {:reqd_work_group_size => [wgs,wgs,1]}) {
      buff_length = size_n+filt.length-1
      buff_length = buff_length.modulo(16) == 0 ? buff_length : buff_length + 16 - buff_length.modulo(16)
      buff_length += 1
      tmp_buff = Variable::new("tmp_buff", Real, { :local => true, :dimension => [Dimension::new( buff_length), Dimension::new(8)] } )
      tmp_buff.decl
      jg = Variable::new("jg", Sizet)
      jg.decl
      (jg === FuncCall::new( "get_global_id", 1)).print
      i2 = Variable::new("i2", Sizet) #, {:constant =>FuncCall::new( "get_local_id", 0)})
      i2.decl
      (i2 === FuncCall::new( "get_local_id", 0)).print
      j2 = Variable::new("j2", Sizet) #, {:constant =>FuncCall::new( "get_local_id", 1)})
      j2.decl
      (j2 === FuncCall::new( "get_local_id", 1)).print
      jgt = Variable::new("jgt", Sizet, {:signed => true})
      jgt.decl
      igt = Variable::new("igt", Sizet, {:signed => true})
      igt.decl
      igtf = Variable::new("igtf", Sizet, {:signed => true})
      igtf.decl
      (jgt === FuncCall::new( "get_group_id", 1)).print
      (jg === (Ternary::new( jgt == FuncCall::new( "get_num_groups", 1) - 1 , jg + ndat - FuncCall::new( "get_global_size", 1) , jg))).print
      replace_constants = $replace_constants
      (igt === j2 + lowfil).print
      (jgt === jg - j2 + i2).print
      #load into buffers
      k = Variable::new("k", Int)
      k.decl
      f1 = For::new(k, 0, size_n+filt.length-1-1-(wgs-1), wgs){
        (igtf === Ternary::new(igt < 0, igt + size_n, Ternary::new(igt >= n, igt - size_n, igt))).print
        (tmp_buff[j2+k,i2] === x[jgt, igtf]).print
        (igt === igt + wgs).print
      }.unroll
      if (size_n+filt.length-1).modulo(wgs) != 0 then
        (igt === j2 + lowfil + size_n+filt.length-1-wgs).print
        (igtf === Ternary::new(igt < 0, igt + size_n, Ternary::new(igt >= n, igt - size_n, igt))).print
        (tmp_buff[j2+(size_n+filt.length-1-wgs),i2] === x[jgt, igtf]).print
      end
      FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print
      tt = Variable::new( "tt", Real)
      tt.decl
      l = Variable::new( "l",Int)
      l.decl
      f1 = For::new(k, 0, size_n-1-(wgs-1), wgs) {
        (tt === 0.0).print
        f2 = For::new(l, 0, filt.length-1) {
          (tt === tt + tmp_buff[k+l+i2, j2]*fil[l]).print
        }.unroll
        ( y[i2 + k, jg] == tt ).print
      }.unroll
      if (size_n.modulo(wgs) != 0) then
        (i2 === i2 + size_n - wgs).print
        (tt === 0.0).print
        f2 = For::new(l, 0, filt.length-1) {
          (tt === tt + tmp_buff[i2+l, j2]*fil[l]).print
        }.unroll
        ( y[i2 , jg] == tt ).print
      end
    }.print
    $array_start = array_start
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



ConvolutionGenerator::set_lang( ConvolutionGenerator::OpenCL )
ConvolutionGenerator::MagicFilter_GPU(FILTER,8,31)
