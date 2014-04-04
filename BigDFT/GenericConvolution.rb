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

  def self.for_conv(nobc,i_in,l,t,processed_dim,unro,init,free,c,scal,x,y,eks,fil,lowfil,upfil,dims,mods)
    t.each_index { |ind|
      i_out = output_index(unro, i_in, ind)
      print t[ind] === ((init and not eks) ? c * x[*i_out] / scal : 0.0)
    }
    if (free and not nobc) then
      loop_start=max(-i_in[processed_dim], lowfil)
      loop_end=min(upfil, dims[processed_dim] - 1 - i_in[processed_dim])
    elsif
      loop_start=lowfil
      loop_end=upfil
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
        print t[ind] === t[ind] + x[*i_out]*fil[l]
      }
    }.unroll#.print#

    t.each_index { |ind|
      i_out = output_index(unro, i_in, ind)
      print t[ind] === t[ind] * scal
      print eks === eks + t[ind] * x[*i_out] if eks
      print y[*i_out] === (init ? t[ind] : y[*i_out] + t[ind] )
    }
  end
end


