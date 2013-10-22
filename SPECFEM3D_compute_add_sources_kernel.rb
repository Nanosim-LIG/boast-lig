module ConvolutionGenerator
  def ConvolutionGenerator::compute_add_sources_kernel
    old_array_start = $array_start
    $array_start = 0
    kernel = CKernel::new
    ConvolutionGenerator::set_output( kernel.code )
    kernel.lang = ConvolutionGenerator::get_lang
    function_name = "compute_add_sources_kernel"
    accel = Variable::new("accel", Real,{:direction => :out, :dimension => [ Dimension::new ]})
    ibool = Variable::new("ibool",Int,{:direction => :in, :dimension => [ Dimension::new ]})
    sourcearrays = Variable::new("sourcearrays", Real,{:direction => :in, :dimension => [ Dimension::new ]})
    stf_pre_compute = Variable::new("stf_pre_compute", Real,{:size => 8, :direction => :in, :dimension => [ Dimension::new ]})
    myrank = Variable::new("myrank",Int,{:direction => :in})
    islice_selected_source = Variable::new("islice_selected_source",Int,{:direction => :in, :dimension => [ Dimension::new ]})
    ispec_selected_source = Variable::new("ispec_selected_source",Int,{:direction => :in, :dimension => [ Dimension::new ]})
    nsources = Variable::new("nsources",Int,{:direction => :in})

    ndim =  Variable::new("NDIM", Int, :constant => 3)
    ngllx =  Variable::new("NGLLX", Int, :constant => 5)
    if kernel.lang == ConvolutionGenerator::CL and ConvolutionGenerator::get_default_real_size == 8 then
      $output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      $output.puts "#pragma OPENCL EXTENSION cl_khr_int64_base_atomics: enable"
    end
    p = Procedure::new(function_name, [accel,ibool,sourcearrays,stf_pre_compute,myrank,islice_selected_source,ispec_selected_source,nsources], [ndim,ngllx])
    if(ConvolutionGenerator::get_lang == ConvolutionGenerator::CUDA) then
      $output.print File::read("specfem3D/#{function_name}.cu")
    elsif(ConvolutionGenerator::get_lang == ConvolutionGenerator::CL) then
      type_f = Real::new.decl
      if ConvolutionGenerator::get_default_real_size == 8 then
        type_i = "unsigned long int"
        cmpx_name = "atom_cmpxchg"
      else
        type_i = "unsigned int"
        cmpx_name = "atomic_cmpxchg"
      end
      $output.print <<EOF
static inline void atomicAdd_f(volatile __global float *source, const float val) {
  union {
    #{type_i} iVal;
    #{type_f} fVal;
  } res, orig;
  do {
    orig.fVal = *source;
    res.fVal = orig.fVal + val;
  } while (#{cmpx_name}((volatile __global #{type_i} *)source, orig.iVal, res.iVal) != orig.iVal);
}
#define INDEX4(xsize,ysize,zsize,x,y,z,i) x + xsize*(y + ysize*(z + zsize*i))
#define INDEX5(xsize,ysize,zsize,isize,x,y,z,i,j) x + xsize*(y + ysize*(z + zsize*(i + isize*(j))))
EOF
      p.decl
      ispec = Variable::new("ispec", Int)
      iglob = Variable::new("iglob", Int)
      stf = Variable::new("stf", Real)
      isource = Variable::new("isource", Int)
      i = Variable::new("i", Int)
      j = Variable::new("j", Int)
      k = Variable::new("k", Int)
      ispec.decl
      iglob.decl
      stf.decl
      isource.decl
      i.decl
      j.decl
      k.decl
      (i === FuncCall::new("get_local_id",0)).print
      (j === FuncCall::new("get_local_id",1)).print
      (k === FuncCall::new("get_local_id",2)).print
      (isource === FuncCall::new("get_group_id",0) + FuncCall::new("get_num_groups",0)*FuncCall::new("get_group_id",1)).print

      If::new(isource < nsources) {
        If::new(myrank == islice_selected_source[isource]) {
          (ispec === ispec_selected_source[isource] - 1).print
          (stf === stf_pre_compute[isource]).print
          (iglob === ibool[FuncCall::new("INDEX4",ngllx,ngllx,ngllx,i,j,k,ispec)] - 1).print
          (0..2).each { |indx|
            (FuncCall::new("atomicAdd_f",accel+iglob*3+indx, sourcearrays[FuncCall::new("INDEX5",ndim,ngllx,ngllx,ngllx,indx,i,j,k,isource)]*stf)).print
          }
        }.print
      }.print

      p.close
    else
      raise "Unsupported language!"
    end
    kernel.procedure = p
    $array_start = old_array_start
    return kernel
  end
end
