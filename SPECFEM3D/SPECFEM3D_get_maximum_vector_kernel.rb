module BOAST
  def BOAST::get_maximum_vector_kernel(ref = true)
    push_env( :array_start => 0 )
    kernel = CKernel::new
    function_name = "get_maximum_vector_kernel"
    size =               Int( "size",               :dir => :in )
    array =              Real("array",              :dir => :in,  :dimension => [ Dim(size)] )
    d_max =              Real("d_max",              :dir => :out, :dimension => [ Dim()] )
    blocksize_transfer = Int( "blocksize_transfer", :const => 256 )
    p = Procedure::new(function_name, [array, size, d_max], [blocksize_transfer])
    if(get_lang == CUDA and ref) then
      @@output.print File::read("specfem3D/#{function_name}.cu")
    elsif(get_lang == CUDA or get_lang == CL) then
      if( get_lang == CL and get_default_real_size == 8) then
        @@output.puts "#pragma OPENCL EXTENSION cl_khr_fp64: enable"
      end
      decl p
      sdata = Real("sdata", :local => true, :dimension => [Dim(blocksize_transfer)] )
      tid =   Int( "tid")
      bx =    Int( "bx")
      i =     Int( "i")
      s =     Int( "s")
      decl sdata
      decl tid
      decl bx
      decl i
      decl s
      print tid === get_local_id(0)
      print bx === get_group_id(1)*get_num_groups(0) + get_group_id(0)
      print i === tid + bx*get_local_size(0)
      print sdata[tid] === Ternary( i < size, sqrt(array[i*3+0]*array[i*3+0]+array[i*3+1]*array[i*3+1]+array[i*3+2]*array[i*3+2]), 0.0)
      print barrier(:local)
      print s === get_local_size(0)/2
      print While(s > 0) {
        print If(tid < s) {
          print If( sdata[tid] < sdata[tid + s] ) {
            print sdata[tid] === sdata[tid + s]
          }
        }
        print s === Expression(">>",s,1)
        print barrier(:local)
      }
      print If(tid == 0) {
        print d_max[bx] === sdata[0]
      }
      close p
    else
      raise "Unsupported language!"
    end
    pop_env( :array_start )
    kernel.procedure = p
    return kernel
  end
end
