require 'BOAST'

module BOAST
  def BOAST::magicfilter_GPU_per_ref
    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::CL)
    kernel = CKernel::new
    kernel.lang = BOAST::CL
    function_name = "magicfilter_per_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y], [],{:reqd_work_group_size => [16,16,1]})
    kernel.code.print <<EOF
#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#elif defined (cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64: enable
#endif
#define FILTER_WIDTH 16
#define FILT0   8.4334247333529341094733325815816e-7
#define FILT1  -0.1290557201342060969516786758559028e-4
#define FILT2   0.8762984476210559564689161894116397e-4
#define FILT3  -0.30158038132690463167163703826169879e-3
#define FILT4   0.174723713672993903449447812749852942e-2
#define FILT5  -0.942047030201080385922711540948195075e-2
#define FILT6   0.2373821463724942397566389712597274535e-1
#define FILT7   0.612625895831207982195380597e-1
#define FILT8   0.9940415697834003993178616713
#define FILT9  -0.604895289196983516002834636e-1
#define FILT10 -0.2103025160930381434955489412839065067e-1
#define FILT11  0.1337263414854794752733423467013220997e-1
#define FILT12 -0.344128144493493857280881509686821861e-2
#define FILT13  0.49443227688689919192282259476750972e-3
#define FILT14 -0.5185986881173432922848639136911487e-4
#define FILT15  2.72734492911979659657715313017228e-6
#define filter(tmp) double tt = 0.0;tt = mad(*tmp++, FILT0, tt);tt = mad(*tmp++, FILT1, tt);tt = mad(*tmp++, FILT2, tt);tt = mad(*tmp++, FILT3, tt);tt = mad(*tmp++, FILT4, tt);tt = mad(*tmp++, FILT5, tt);tt = mad(*tmp++, FILT6, tt);tt = mad(*tmp++, FILT7, tt);tt = mad(*tmp++, FILT8, tt);tt = mad(*tmp++, FILT9, tt);tt = mad(*tmp++, FILT10, tt);tt = mad(*tmp++, FILT11, tt);tt = mad(*tmp++, FILT12, tt);tt = mad(*tmp++, FILT13, tt);tt = mad(*tmp++, FILT14, tt);tt = mad(*tmp++, FILT15, tt);
//n is supposed to be greater or equal than get_local_size(0)
//this filter is for periodic boundary conditions
__kernel __attribute__((reqd_work_group_size(16,16, 1))) void #{function_name}(uint n, uint ndat, __global const double * restrict psi, __global double * restrict out){
__local double tmp1[FILTER_WIDTH*(2*FILTER_WIDTH+1)];
__local double *tmp = &tmp1[0];
//get our position in the local work group
size_t ig = get_global_id(0);
size_t jg = get_global_id(1);
//get our position in the result matrix
const size_t i2 = get_local_id(0);
const size_t j2 = get_local_id(1);
//get our group number
ptrdiff_t igt = get_group_id(0);
ptrdiff_t jgt = get_group_id(1);
//if data are ill dimentioned border blocks recomputes part of the data
jg  = jgt == get_num_groups(1) - 1 ? jg - ( get_global_size(1) - ndat ) : jg;
ig  = igt == get_num_groups(0) - 1 ? ig - ( get_global_size(0) - n ) : ig;
//transpose indexes in the work group in order to read transposed data
igt = ig - i2 + j2 - FILTER_WIDTH/2;
jgt = jg - j2 + i2;
//if we are on the outside, select a border element to load, wrapping around
//we will be loading 2 elements each
if ( igt < 0 ) 
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + ( n + igt ) * ndat];
else 
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2] = psi[jgt + igt * ndat];
igt += FILTER_WIDTH;
if ( igt >= n ) 
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt + ( igt - n ) * ndat];
else
  tmp[i2 * (2 * FILTER_WIDTH + 1) + j2 + FILTER_WIDTH] = psi[jgt +  igt * ndat];
//rest position in the buffer to first element involved in the convolution
tmp += j2*(2*FILTER_WIDTH+1) + i2;
//wait for buffer to be full
barrier(CLK_LOCAL_MEM_FENCE);
//apply filter
filter(tmp);
//store the result
out[(jg*n+ig)]=tt;
};
EOF
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end

  def BOAST::magicfilter_GPU(filt, center, size_n, max_work_item=256, local_mem_size=16384 )
    array_start = @@array_start
    @@array_start = 0

    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::CL)
    kernel = CKernel::new
    kernel.lang = BOAST::CL
    BOAST::set_output( kernel.code )

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
    @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
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
      (igt === j2 + lowfil).print
      (jgt === jg - j2 + i2).print
      #load into buffers
      k = Variable::new("k", Int)
      k.decl

      (igtf === Ternary::new(igt < 0, igt + size_n, Ternary::new(igt >= n, igt - size_n, igt))).print
      (tmp_buff[j2+0,i2] === x[jgt, igtf]).print
      (igt === igt + wgs).print
      f1 = For::new(k, filt.length/2, size_n+filt.length/2-1-1-(wgs-1), wgs){
        (igtf === igt).print
        (tmp_buff[j2+k,i2] === x[jgt, igtf]).print
        (igt === igt + wgs).print
      }.unroll
      (igtf === Ternary::new(igt < 0, igt + size_n, Ternary::new(igt >= n, igt - size_n, igt))).print
      (tmp_buff[j2+(size_n+filt.length-1-1-(wgs-1)-(size_n+filt.length-1-1-(wgs-1)).modulo(wgs)),i2] === x[jgt, igtf]).print
      (igt === igt + wgs).print

      if (size_n+filt.length-1).modulo(wgs) != 0 then
        (igt === j2 + lowfil + size_n+filt.length-1-wgs).print
        (igtf === Ternary::new(igt < 0, igt + size_n, Ternary::new(igt >= n, igt - size_n, igt))).print
        (tmp_buff[j2+(size_n+filt.length-1-wgs),i2] === x[jgt, igtf]).print
      end
      FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print
      fil.decl
      tt = Variable::new( "tt", Real)
      tt.decl
      l = Variable::new( "l",Int)
      l.decl
      f1 = For::new(k, 0, size_n-1-(wgs-1), wgs) {
        (tt === 0.0).print
        f2 = For::new(l, 0, filt.length-1) {
          (tt === tt + tmp_buff[k+l+i2, j2]*fil[l]).print
        }.unroll
        ( y[i2 + k, jg] === tt ).print
      }.unroll
      if (size_n.modulo(wgs) != 0) then
        (i2 === i2 + size_n - wgs).print
        (tt === 0.0).print
        f2 = For::new(l, 0, filt.length-1) {
          (tt === tt + tmp_buff[i2+l, j2]*fil[l]).print
        }.unroll
        ( y[i2 , jg] === tt ).print
      end
    }
    p.print
    @@array_start = array_start
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end

  def BOAST::magicfilter_GPU_next(filt, center, size_n, max_work_item=256, local_mem_size=16384 )
    array_start = @@array_start
    @@array_start = 0

    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::CL)
    kernel = CKernel::new
    kernel.lang = BOAST::CL
    BOAST::set_output( kernel.code )

    function_name = "magicfilter_n#{size_n}"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
    upfil = Variable::new("upfil",Int,{:constant => center-1})
    dim_in_min = 0
    dim_in_max = size_n-1
    dim_out_min = 0
    dim_out_max = size_n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(ndat), Dimension::new(dim_in_min, dim_in_max) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(dim_out_min, dim_out_max), Dimension::new(ndat) ] })
    arr = ConstArray::new(filt)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil, upfil) ]})
    
    wgs = 16
    @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil], {:reqd_work_group_size => [wgs,wgs,1]}) {
      buff_length = size_n
      buff_length = buff_length.modulo(16) == 0 ? buff_length : buff_length + 16 - buff_length.modulo(16)
      buff_length += 1
      tmp_buff = Variable::new("tmp_buff", Real, { :local => true, :dimension => [Dimension::new( buff_length), Dimension::new(wgs)] } )
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
      (jgt === FuncCall::new( "get_group_id", 1)).print
      (jg === (Ternary::new( jgt == FuncCall::new( "get_num_groups", 1) - 1 , jg + ndat - FuncCall::new( "get_global_size", 1) , jg))).print
      (jgt === jg - j2 + i2).print
      #load into buffers
      k = Variable::new("k", Int)

      if size_n > wgs then
        (igt === j2 + size_n - wgs).print
        (tmp_buff[igt,i2] === x[jgt, igt]).print
      end
      (igt === j2).print
      (tmp_buff[igt,i2] === x[jgt, igt]).print
      (igt === igt + wgs).print
      For::new(k, wgs, size_n-wgs-1, wgs){
        (tmp_buff[igt,i2] === x[jgt, igt]).print
        (igt === igt + wgs).print
      }.unroll

      FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print

      tt = Variable::new( "tt", Real)
      tt.decl
      l = Variable::new( "l",Int)
      ig = Variable::new("ig", Sizet, {:signed => true})
      ig.decl
      (ig === i2 - wgs/2).print

      (ig === Ternary::new( ig < 0 , ig + size_n , ig ) ).print
        
      (tt === 0.0).print
      For::new( l, lowfil, -1) {
        (igt === ig+l).print
        (igt === Ternary::new(igt < 0, igt + size_n, igt) ).print 
        (tt === tt + tmp_buff[igt, j2]*fil[l]).print
      }.unroll
      For::new( l, 0, upfil) {
        (igt === ig+l).print
        (igt === Ternary::new( igt >= size_n, igt -size_n, igt) ).print 
        (tt === tt + tmp_buff[igt, j2]*fil[l]).print
      }.unroll
      (y[ig, jg] === tt ).print
      
      (ig === i2 + wgs/2).print
      ((wgs/2)..(size_n-wgs/2-wgs)).step(wgs) { |k|
#        if(k+wgs/2 <= size_n-wgs-1) then
#          (igt === j2 + k+wgs/2).print
#          (tmp_buff[igt,i2] === x[jgt, igt]).print
#          FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print
#        end
        (igt === ig + lowfil).print
        (tt === 0.0).print
        f2 = For::new(l, lowfil, upfil) {
          (tt === tt + tmp_buff[igt.inc, j2]*fil[l]).print
        }.unroll
        ( y[ig, jg] === tt ).print
        (ig === ig + wgs).print
      }

      if (size_n.modulo(wgs) != 0) then
        (ig === i2 + size_n - wgs - wgs/2).print
        (tt === 0.0).print
        if (size_n < wgs+filt.length) then
          (ig === Ternary::new( ig < 0 , ig + size_n , ig ) ).print
          For::new( l, lowfil, -1) {
            (igt === ig+l).print
            (igt === Ternary::new(igt < 0, igt + size_n, igt) ).print 
            (tt === tt + tmp_buff[igt, j2]*fil[l]).print
          }.unroll
          For::new( l, 0, upfil) {
            (igt === ig+l).print
            (igt === Ternary::new( igt >= size_n, igt - size_n, igt) ).print 
            (tt === tt + tmp_buff[igt, j2]*fil[l]).print
          }.unroll
        else
          (igt === ig + lowfil).print
          f2 = For::new(l, lowfil, upfil) {
            (tt === tt + tmp_buff[igt.inc, j2]*fil[l]).print
          }.unroll
        end
        ( y[ig , jg] === tt ).print
      end
    }
    p.print
    @@array_start = array_start
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end

  def BOAST::magicfilter_GPU_next_next(filt, center, size_n, max_work_item=256, local_mem_size=16384 )
    array_start = @@array_start
    @@array_start = 0

    lang = BOAST::get_lang
    BOAST::set_lang(BOAST::CL)
    kernel = CKernel::new
    kernel.lang = BOAST::CL
    BOAST::set_output( kernel.code )

    function_name = "magicfilter_n#{size_n}"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
    upfil = Variable::new("upfil",Int,{:constant => center-1})
    dim_in_min = 0
    dim_in_max = size_n-1
    dim_out_min = 0
    dim_out_max = size_n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(ndat), Dimension::new(dim_in_min, dim_in_max) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(dim_out_min, dim_out_max), Dimension::new(ndat) ] })
    arr = ConstArray::new(filt)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil, upfil) ]})
    
    wgs = 16
    @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil], {:reqd_work_group_size => [wgs,wgs,1]}) {
      buff_length = size_n
      buff_length = buff_length.modulo(16) == 0 ? buff_length : buff_length + 16 - buff_length.modulo(16)
      buff_length += 1
      tmp_buff = Variable::new("tmp_buff", Real, { :local => true, :dimension => [Dimension::new( buff_length), Dimension::new(wgs)] } )
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
      (jgt === FuncCall::new( "get_group_id", 1)).print
      (jg === (Ternary::new( jgt == FuncCall::new( "get_num_groups", 1) - 1 , jg + ndat - FuncCall::new( "get_global_size", 1) , jg))).print
      (jgt === jg - j2 + i2).print
      #load into buffers
      k = Variable::new("k", Int)

      if size_n > wgs then
        (igt === j2 + size_n - wgs).print
        (tmp_buff[igt,i2] === x[jgt, igt]).print
      end
      (igt === j2).print
      (tmp_buff[igt,i2] === x[jgt, igt]).print
      (igt === igt + wgs).print
#      For::new(k, wgs, size_n-wgs-1, wgs){
#        (tmp_buff[igt,i2] === x[jgt, igt]).print
#        (igt === igt + wgs).print
#      }.unroll

      FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print

      tt = Variable::new( "tt", Real)
      tt.decl
      l = Variable::new( "l",Int)
      ig = Variable::new("ig", Sizet, {:signed => true})
      ig.decl
      (ig === i2 - wgs/2).print

      (ig === Ternary::new( ig < 0 , ig + size_n , ig ) ).print
        
      (tt === 0.0).print
      For::new( l, lowfil, -1) {
        (igt === ig+l).print
        (igt === Ternary::new(igt < 0, igt + size_n, igt) ).print 
        (tt === tt + tmp_buff[igt, j2]*fil[l]).print
      }.unroll
      For::new( l, 0, upfil) {
        (igt === ig+l).print
        (igt === Ternary::new( igt >= size_n, igt -size_n, igt) ).print 
        (tt === tt + tmp_buff[igt, j2]*fil[l]).print
      }.unroll
      (y[ig, jg] === tt ).print
      
      (ig === i2 + wgs/2).print
      ((wgs/2)..(size_n-wgs/2-wgs)).step(wgs) { |k|
        if(k+wgs/2 <= size_n-wgs-1) then
          (igt === j2 + k+wgs/2).print
          (tmp_buff[igt,i2] === x[jgt, igt]).print
          FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print
        end
        (igt === ig + lowfil).print
        (tt === 0.0).print
        f2 = For::new(l, lowfil, upfil) {
          (tt === tt + tmp_buff[igt.inc, j2]*fil[l]).print
        }.unroll
        ( y[ig, jg] === tt ).print
        (ig === ig + wgs).print
      }

      if (size_n.modulo(wgs) != 0) then
        (ig === i2 + size_n - wgs - wgs/2).print
        (tt === 0.0).print
        if (size_n < wgs+filt.length) then
          (ig === Ternary::new( ig < 0 , ig + size_n , ig ) ).print
          For::new( l, lowfil, -1) {
            (igt === ig+l).print
            (igt === Ternary::new(igt < 0, igt + size_n, igt) ).print 
            (tt === tt + tmp_buff[igt, j2]*fil[l]).print
          }.unroll
          For::new( l, 0, upfil) {
            (igt === ig+l).print
            (igt === Ternary::new( igt >= size_n, igt - size_n, igt) ).print 
            (tt === tt + tmp_buff[igt, j2]*fil[l]).print
          }.unroll
        else
          (igt === ig + lowfil).print
          f2 = For::new(l, lowfil, upfil) {
            (tt === tt + tmp_buff[igt.inc, j2]*fil[l]).print
          }.unroll
        end
        ( y[ig , jg] === tt ).print
      end
    }
    p.print
    @@array_start = array_start
    kernel.procedure = p
    BOAST::set_lang(lang)
    return kernel
  end



end


