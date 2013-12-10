module BOAST
  def BOAST::get_maximum_scalar_kernel
    old_array_start = @@array_start
    @@array_start = 0
    kernel = CKernel::new
    function_name = "get_maximum_scalar_kernel"
    size = Variable::new("size",Int,{:direction => :in})
    array = Variable::new("array", Real,{:direction => :in, :dimension => [ Dimension::new(size)]})
    d_max = Variable::new("d_max", Real,{:direction => :out, :dimension => [ Dimension::new]})
    blocksize_transfer =  Variable::new("blocksize_transfer", Int, :constant => 256)
    if kernel.lang == BOAST::CL and BOAST::get_default_real_size == 8 then
      @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
    end
    p = Procedure::new(function_name, [array, size, d_max], [blocksize_transfer])
    if(BOAST::get_lang == BOAST::CUDA) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(BOAST::get_lang == BOAST::CL) then
      p.decl
      sdata = Variable::new("sdata", Real, { :local => true, :dimension => [Dimension::new(blocksize_transfer)] } )
      tid = Variable::new("tid", Int)
      bx = Variable::new("bx", Int)
      i = Variable::new("i", Int)
      s = Variable::new("s", Int)
      sdata.decl
      tid.decl
      bx.decl
      i.decl
      s.decl
      (tid === FuncCall::new("get_local_id",0)).print
      (bx === FuncCall::new("get_group_id",1)*FuncCall::new("get_num_groups",0) + FuncCall::new("get_group_id",0)).print
      (i === tid + bx*FuncCall::new("get_local_size",0)).print;
      (sdata[tid] === Ternary::new( i < size, FuncCall::new("fabs", array[i]), 0.0)).print;
      FuncCall::new( "barrier","CLK_LOCAL_MEM_FENCE").print
      (s === Expression::new("/",FuncCall::new("get_local_size",0),2)).print
      While::new(s > 0) {
        If::new(tid < s) {
          If::new( sdata[tid] < sdata[tid + s] ) {
            (sdata[tid] === sdata[tid + s]).print
          }.print
        }.print
        (s === Expression::new(">>",s,1)).print
      }.print
      If::new(tid == 0) {
        (d_max[bx] === sdata[0]).print
      }.print
      p.close
    else
      raise "Unsupported language!"
    end
    kernel.procedure = p
    @@array_start = old_array_start
    return kernel
  end
end
