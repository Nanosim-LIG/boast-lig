module BOAST
  @@ocl_cuda_dim_assoc = { 0 => "x", 1 => "y", 2 => "z" }


  def BOAST::get_work_dim
    if @@lang == CL then
      return FuncCall::new("get_work_dim")
    else
      raise "Unsupported language!"
    end
  end
  
  def BOAST::get_global_size(dim)
    if @@lang == CL then
      return FuncCall::new("get_global_size",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "gridDim", d)*Expression::new(".", "blockDim", d)
    else
      raise "Unsupported language!"
    end
  end

  def BOAST::get_global_id(dim)
    if @@lang == CL then
      return FuncCall::new("get_global_id",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "threadIdx", d)+Expression::new(".", "blockIdx", d)*Expression::new(".", "blockDim", d)
    else
      raise "Unsupported language!"
    end
  end

  def BOAST::get_local_size(dim)
    if @@lang == CL then
      return FuncCall::new("get_local_size",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "blockDim", d)
    else
      raise "Unsupported language!"
    end
  end

  def BOAST::get_local_id(dim)
    if @@lang == CL then
      return FuncCall::new("get_local_id",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "threadIdx", d)
    else
      raise "Unsupported language!"
    end
  end
  
  def BOAST::get_num_groups(dim)
    if @@lang == CL then
      return FuncCall::new("get_num_groups",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "gridDim", d)
    else
      raise "Unsupported language!"
    end
  end

  def BOAST::get_group_id(dim)
    if @@lang == CL then
      return FuncCall::new("get_group_id",dim)
    elsif @@lang == CUDA then
      d = @@ocl_cuda_dim_assoc[dim]
      raise "Unsupported dimension!" if not d
      return Expression::new(".", "blockIdx", d)
    else
      raise "Unsupported language!"
    end
  end
  
end
