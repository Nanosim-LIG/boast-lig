module BOAST

  def self.filter_boastruct(filt,center)
    arr = ConstArray::new(filt,Real)
    fil=Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((-center),(filt.length - center -1)) ]})
    lowfil = Variable::new("lowfil",Int,{:constant => -center})
    upfil = Variable::new("upfil",Int,{:constant => filt.length - center -1})
    return [fil,lowfil,upfil]
  end

  #returns the indices of the output according to which of the directions is unrolled
  def self.output_index_unroll(unrolling_dim, i_in,unroll_index)
    i_out=(0...i_in.length).collect { |indx| unrolling_dim == indx ? i_in[unrolling_dim] + (unroll_index) : i_in[indx]}
    return i_out
  end
  def self.output_index_k(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index)
    i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
    return (0...i_in.length).collect { |indx| processed_dim == indx ? lconv_index +i_in[processed_dim] : i_out[indx]}
  end
  def self.output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array)
    i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
    return (0...i_in.length).collect { |indx| processed_dim == indx ? wrapping_array[lconv_index +i_in[processed_dim]] : i_out[indx]}
  end
  def self.output_index_k_mod(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,ndim_processed)
    i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
    return (0...i_in.length).collect { |indx| processed_dim == indx ?  lconv_index + i_in[processed_dim] - ((i_in[processed_dim]+lconv_index +  ndim_processed * 2 )/ndim_processed - 2) *ndim_processed  : i_out[indx]}
  end

  #returns the indices of the output array according to the starting point in the input and of the
  ## processed dimension as well as th eposition in the convolution
  def self.output_index(unrolling_dim, i_in,unroll_index,processed_dim=nil,lconv_index=nil,
                        ndim_processed=nil,wrapping_array=nil)
    if ndim_processed then
      i_out=output_index_k_mod(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,ndim_processed)
    elsif wrapping_array then
      i_out=output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array)
    elsif processed_dim then
      i_out=output_index_k(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index)
    else
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
    end
    return i_out
  end

  def self.for_conv(nobc,i_in,l,t,processed_dim,unro,init,free,c,scal,x,y,eks,filter,dims,mods,nrotate)
    t.each_index { |ind|
      i_out = output_index(unro, i_in, ind)
      #WARNING: the eks conditional here can be relaxed
      print t[ind] === ((init and not eks) ? c * x[*i_out] / scal : 0.0)
    }
    if (free and not nobc) then
      loop_start=max(-i_in[processed_dim], filter.lowfil)
      loop_end=min(filter.upfil, dims[processed_dim] - 1 - i_in[processed_dim])
    elsif
      loop_start=filter.lowfil
      loop_end=filter.upfil
    end
    For( l,loop_start,loop_end) {
      t.each_index { |ind|
        if free or nobc then
          i_out = output_index(unro, i_in, ind,processed_dim,l)
        elsif mods then
          i_out = output_index(unro, i_in, ind,processed_dim,l,nil,mods) 
        else
          i_out = output_index(unro, i_in, ind,processed_dim,l,dims[processed_dim])
        end
        print t[ind] === t[ind] + x[*i_out]*filter.fil[l]
      }
    }.unroll#.print#

    t.each_index { |ind|
      i_out = output_index(unro, i_in, ind)
      print t[ind] === t[ind] * scal
      print eks === eks + t[ind] * x[*i_out] if eks
      print y[*i_out.rotate(nrotate)] === (init ? t[ind] : y[*i_out] + t[ind] )
    }
  end

  class ConvolutionFilter
    # List of the floating point values of the convolution filter
    attr_reader :fil_array
    # central point of the filter
    attr_reader :center
    # BOAST variables
    # Filter array (to be used on BOAST functions)
    attr_reader :fil
    # extremes of the filter, calculated via its central point (BOAST object)
    attr_reader :lowfil, :upfil

    def initialize(filt,center)
      arr = ConstArray::new(filt,Real)
      @fil=Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((-center),(filt.length - center -1)) ]})
      @lowfil = Variable::new("lowfil",Int,{:constant => -center})
      @upfil = Variable::new("upfil",Int,{:constant => filt.length - center -1})
      @fil_array=filt
      @center=center
    end

  end

  def self.conv_lines(filter,free,dims,iters,l,t,processed_dim,unro,nrotate,init,c,scal,x,y,eks,mods)
    #to be checked: grow and shrink operations have to be applied there!
    For::new(iters[processed_dim], 0, -filter.lowfil) {
      for_conv(false,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
               dims,mods,nrotate)
    }.print
    For::new(iters[processed_dim], -filter.lowfil+1, dims[processed_dim] - 1 - filter.upfil) {
      for_conv(true,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
               dims,mods,nrotate)
    }.print
    #to be checked: grow and shrink operations have to be applied there!
    For::new(iters[processed_dim], dims[processed_dim] - filter.upfil, dims[processed_dim]-1) {
      for_conv(false,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
               dims,mods,nrotate)
    }.print
  end

  #returns the starting and ending points of the convolutions according to unrolling and unrolled dimension
  def self.startendpoints(dim,unroll,unrolling_length,in_reliq)
    istart= (in_reliq and unroll) ? (dim/unrolling_length)*unrolling_length : 0
    iend  = (unroll and not in_reliq) ? dim-unrolling_length : dim-1
    istep = (unroll and not in_reliq) ? unrolling_length : 1 
    return [istart,iend,istep]
  end

  class ConvolutionOperator
    # Class describing to the convolution operator. Should contain the filter values and its central point
    attr_reader :filter
    # Dimensions of the problem. Can be multi-dimensional
    attr_reader :n
    # Boundary conditions of the problem Integer values: -1 (shrink), 0 (periodic), 1 (grow) - to be implemented in for_conv
    attr_reader :bc
    # Constants (in matricial sense) to which the output arrays has to be initialized, y <- alpha* Conv * x + beta * x
    attr_reader :beta, :a
    # Matrix where the results of the convolution reductions can be written
    attr_reader :eks

    def initialize(filter,n,bc,beta,a,eks)
      @filter = filter
      @n = n
      @bc = bc
      @beta = beta
      @a = a
      @eks = eks
    end
    
    def dim_indexes(processed_dim, transpose)
      raise "Invalid processed_dim: #{processed_dim}" if processed_dim >= @n.length or processed_dim < 0 
      if transpose then
        if processed_dim == 0 then
          return [1, 0]
        else
          return [0, 1]
        end
      else
        if processed_dim == 0 then
          return [1, 0]
        elsif processed_dim == @n.length - 1
          return [0, 1]
        else
          return [2,0,1]
        end
      end
    end
  end

  def self.convolution1d(filter,dims,free,iters,l,t,dim_indexes,nrotate,init,c,scal,x,y,eks,mods,unro,unrolling_length)
    convgen= lambda { |dim,t,reliq|
      ises0 = startendpoints(dims[dim[0]],unro == dim[0],unrolling_length,reliq)
      @@output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
      For::new(iters[dim[0]], *ises0 ) {
        if dim.length == 3 then
          ises1 = startendpoints(dims[dim[1]],unro == dim[1],unrolling_length,reliq)
          For::new(iters[dim[1]], *ises1) {
            conv_lines(filter,free,dims,iters,l,t,dim[-1],unro,nrotate,init,c,scal,x,y,eks,mods)
          }.print
        else
          conv_lines(filter,free,dims,iters,l,t,dim[-1],unro,nrotate,init,c,scal,x,y,eks,mods)
        end
      }.print
      @@output.print("!$omp end do\n") if BOAST::get_lang == BOAST::FORTRAN
    }
    #first without the reliq
    convgen.call(dim_indexes,t[0..unrolling_length-1],false)
    #then with the reliq but only if the unrolling patterns need it
    convgen.call(dim_indexes,t[0..0],true) if (unrolling_length > 1)
  end


end


