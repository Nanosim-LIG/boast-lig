module BOAST

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
    # name of the filter
    attr_reader :name

    def initialize(name,filt,center)
      arr = ConstArray::new(filt,Real)
      @fil=BOAST::Real("fil",:constant => arr,:dim => [ BOAST::Dim((-center),(filt.length - center -1)) ])
      @lowfil = BOAST::Int("lowfil",:constant => -center)
      @upfil = BOAST::Int("upfil",:constant => filt.length - center -1)
      @fil_array=filt
      @center=center
      @name=name
    end

  end

  class ConvolutionOperator1d
    # Convolution filter
    attr_reader :filter
    # dimensions
    attr_reader :dims, :dim_n
    # input array, unchanged on exit
    attr_reader :in
    # output array, reduced on exit
    attr_reader :out
    # out <- [out +] alpha conv(in) + beta in
    # reduce the array or not
    attr_reader :reduce
    # constants of the convolution
    attr_reader :beta, :alpha
    # value of <in | conv (in)>
    attr_reader :dotp
    # order of the dimensions to be convolved (las one is the processed dim)
    attr_reader :dim_indexes
    # variables of the procedure
    attr_reader :vars
    # options
    attr_reader :option
    # Creates new 1d convolution
    # 
    # ==== Attributes
    # 
    # * +filter+ - ConvolutionFilter object corresponding to the operations to be applied on data
    # * +bc+ Boundary conditions: control the way in which the convolution has to be applied.
    #        Typical values are 0 : periodic BC, the size of input and output arrays are identical
    #                           1 : Free BC, grow: the size of the output array is equal to the 
    #                                              size of input array plus the one of the filter
    #                          -1 : Free BC, shrink: the size of the input array is equal to the one
    #                                        of output array plus the filter
    #                     Given a convolution and its inverse, in free BC -1 is the inverse of 1
    #                      but not viceversa as the number of point treated is lower.
    #                         10:  Free BC, the size of input and output arrays are identical
    #                              (loss of information)
    #                       
    # * +options+ - Hash table of allowed options (see options descritpion)
    #
    # ==== Options
    #
    # * +:alpha+ -  Convolution constant of the operation out <- [out +] alpha conv(in) + beta in
    # * +:beta+ -  Convolution constant of the operation out <- [out +] alpha conv(in) + beta in
    def initialize(filter,bc,transpose,dim_indexes,init,options={})
      @filter = filter
      @bc = bc
      @transpose = transpose
      @dim_indexes = dim_indexes
      @dim_n = BOAST::Int("n",:dir =>:in)
      @dims = [@dim_n]
      if (dim_indexes.length == 3) then
        @dims = [BOAST::Int("ndat1",:dir =>:in)] + @dims + [BOAST::Int("ndat2",:dir =>:in)]
      elsif dim_indexes.last == 0
        @dims = @dims + [BOAST::Int("ndat",:dir =>:in)]
      else
        @dims = [BOAST::Int("ndat",:dir =>:in)] + @dims
      end
      @vars = @dims.dup
      dimx =  @dims.collect{ |dim|
        BOAST::Dim(0, dim-1)
      }
      if transpose !=0  then
        dimy = dimx.rotate(1)
      else
        dimy= dimx
      end
      @vars.push @in = BOAST::Real("x",:dir => :in, :dim => dimx, :restrict => true)
      @vars.push @out = BOAST::Real("y",:dir => :out, :dim => dimy, :restrict => true)
      @vars.push @alpha = BOAST::Real("alpha",:dir => :in) if options[:alpha]
      @vars.push @beta = BOAST::Real("beta",:dir => :in) if options[:beta] and init
      @vars.push @dotp = BOAST::Real("dotp",:dir => :out) if options[:dotp]
      @init = init
      @options = options
    end

    def procedure(unroll,unrolled_dim,use_mod)
      function_name = @filter.name + "_" + ( (@bc != 0 and @bc != 10) ? "f" : "p") + "_#{@dim_indexes.join('')}" +
        "_u#{unroll}_#{unrolled_dim}_#{use_mod}"

      l = BOAST::Int("l")
      tt = (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("tt#{index}") }
      iters =  (1..@dims.length).collect{ |index| BOAST::Int("i#{index}")}
      
      if use_mod then
        mods=BOAST::Int("mod_arr", :allocate => true, :dim => [BOAST::Dim(@filter.lowfil,@dim_n+@filter.upfil)])
      else
        mods=nil
      end
      return BOAST::Procedure(function_name,vars,[@filter.lowfil,@filter.upfil]){
        BOAST::decl @filter.fil
        BOAST::decl *iters
        BOAST::decl l
        BOAST::decl *tt
        if use_mod then
          BOAST::decl mods 
          BOAST::print BOAST::For(l, @filter.lowfil, @dim_n + @filter.upfil) {
            BOAST::print mods[l] === BOAST::modulo(l, @dim_n)
          } 
        end
        BOAST::print @dotp === 0.0 if @options[:dotp]
        if BOAST::get_lang == BOAST::FORTRAN then
          BOAST::get_output.print("!$omp parallel default(shared)&\n")
          BOAST::get_output.print("!$omp reduction(+:#{dotp})&\n") if @options[:dotp]
          BOAST::get_output.print("!$omp private(#{iters.join(",")},#{l})&\n")
          BOAST::get_output.print("!$omp private(#{tt.join(",")})\n")
        elsif BOAST::get_lang == BOAST::C then
          BOAST::get_output.print("#pragma omp parallel default(shared) #{@options[:dotp] ? "reduction(+:#{dotp})" : ""} private(#{iters.join(",")},#{l},#{tt.join(",")})\n")
        end

        convolution1d(@filter,@dims,@bc,iters,l,tt,@dim_indexes,@transpose,@init,
                            @alpha,@beta,@in,@out,@dotp,mods,unrolled_dim,unroll)

        BOAST::get_output.print("!$omp end parallel\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("#pragma omp end parrallel\n")  if BOAST::get_lang == BOAST::C
      }
    end

    #here follows the internal operations for the convolution 1d
    def convolution1d(filter,dims,free,iters,l,t,dim_indexes,nrotate,init,scal,c,x,y,eks,mods,unro,unrolling_length)
      convgen= lambda { |dim,t,reliq|
        ises0 = startendpoints(dims[dim[0]],unro == dim[0],unrolling_length,reliq)
        BOAST::get_output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("#pragma omp for\n")
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
        BOAST::get_output.print("!$omp end do\n") if BOAST::get_lang == BOAST::FORTRAN
      }
      #first without the reliq
      convgen.call(dim_indexes,t[0..unrolling_length-1],false)
      #then with the reliq but only if the unrolling patterns need it
      convgen.call(dim_indexes,t[0..0],true) if (unrolling_length > 1)
    end

    #returns the starting and ending points of the convolutions according to unrolling and unrolled dimension
    def startendpoints(dim,unroll,unrolling_length,in_reliq)
      istart= (in_reliq and unroll) ? (dim/unrolling_length)*unrolling_length : 0
      iend  = (unroll and not in_reliq) ? dim-unrolling_length : dim-1
      istep = (unroll and not in_reliq) ? unrolling_length : 1 
      return [istart,iend,istep]
    end

    def conv_lines(filter,free,dims,iters,l,t,processed_dim,unro,nrotate,init,c,scal,x,y,eks,mods)
      #to be checked: grow and shrink operations have to be applied there!
      BOAST::For(iters[processed_dim], 0, -filter.lowfil) {
        for_conv(false,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
                 dims,mods,nrotate)
      }.print
      BOAST::For(iters[processed_dim], -filter.lowfil+1, dims[processed_dim] - 1 - filter.upfil) {
        for_conv(true,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
                 dims,mods,nrotate)
      }.print
      #to be checked: grow and shrink operations have to be applied there!
      BOAST::For(iters[processed_dim], dims[processed_dim] - filter.upfil, dims[processed_dim]-1) {
        for_conv(false,iters,l,t, processed_dim, unro, init,free,c,scal,x,y,eks,filter,
                 dims,mods,nrotate)
      }.print
    end

    def for_conv(nobc,i_in,l,t,processed_dim,unro,init,free,c,scal,x,y,eks,filter,dims,mods,nrotate)
      t.each_index { |ind|
        i_out = output_index(unro, i_in, ind)
        #WARNING: the eks conditional here can be relaxed
        BOAST::print t[ind] === ((init and not eks) ? c * x[*i_out] / scal : 0.0)
      }
      if ((free !=0 and free != 10) and not nobc) then
        loop_start=max(-i_in[processed_dim], filter.lowfil)
        loop_end=min(filter.upfil, dims[processed_dim] - 1 - i_in[processed_dim])
      elsif
        loop_start=filter.lowfil
        loop_end=filter.upfil
      end
      BOAST::For( l,loop_start,loop_end) {
        t.each_index { |ind|
          if (free != 0 and free !=10) or nobc then
            i_out = output_index(unro, i_in, ind,processed_dim,l)
          elsif mods then
            i_out = output_index(unro, i_in, ind,processed_dim,l,nil,mods) 
          else
            i_out = output_index(unro, i_in, ind,processed_dim,l,dims[processed_dim])
          end
          BOAST::print t[ind] === t[ind] + x[*i_out]*filter.fil[l]
        }
      }.unroll#.print#

      t.each_index { |ind|
        i_out = output_index(unro, i_in, ind)
        BOAST::print t[ind] === t[ind] * scal if scal
        BOAST::print eks === eks + t[ind] * x[*i_out] if eks
        if nrotate != 0 then
          BOAST::print y[*i_out.rotate(nrotate)] === t[ind]
        else
          BOAST::print y[*i_out.rotate(nrotate)] === 
            (init ? t[ind] : y[*i_out.rotate(nrotate)] + t[ind] )
        end
          
      }
    end

    #returns the indices of the output array according to the starting point in the input and of the
    ## processed dimension as well as th eposition in the convolution
    def output_index(unrolling_dim, i_in,unroll_index,processed_dim=nil,lconv_index=nil,
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

    def filter_boastruct(filt,center)
      arr = ConstArray::new(filt,Real)
      fil=Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new((-center),(filt.length - center -1)) ]})
      lowfil = Variable::new("lowfil",Int,{:constant => -center})
      upfil = Variable::new("upfil",Int,{:constant => filt.length - center -1})
      return [fil,lowfil,upfil]
    end

    #returns the indices of the output according to which of the directions is unrolled
    def output_index_unroll(unrolling_dim, i_in,unroll_index)
      i_out=(0...i_in.length).collect { |indx| unrolling_dim == indx ? i_in[unrolling_dim] + (unroll_index) : i_in[indx]}
      return i_out
    end
    def output_index_k(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      return (0...i_in.length).collect { |indx| processed_dim == indx ? lconv_index +i_in[processed_dim] : i_out[indx]}
    end
    def output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      return (0...i_in.length).collect { |indx| processed_dim == indx ? wrapping_array[lconv_index +i_in[processed_dim]] : i_out[indx]}
    end
    def output_index_k_mod(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,ndim_processed)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      return (0...i_in.length).collect { |indx| processed_dim == indx ?  lconv_index + i_in[processed_dim] - ((i_in[processed_dim]+lconv_index +  ndim_processed * 2 )/ndim_processed - 2) *ndim_processed  : i_out[indx]}
    end


  end


  class ConvolutionOperator
    # Class describing to the convolution operator. Should contain the filter values and its central point
    attr_reader :filter
    # Dimensions of the problem. Can be multi-dimensional
    attr_reader :ndim
    # Boundary conditions of the problem Integer values: -1 (shrink), 0 (periodic), 1 (grow) - to be implemented in for_conv
    attr_reader :bc
    # Constants (in matricial sense) to which the output arrays has to be initialized, 
    # y <- alpha* Conv * x + beta * x
    attr_reader :beta, :alpha
    # Matrix where the results of the convolution reductions can be written
    attr_reader :eks
    # Array of boast dimensions
    attr_reader :dims
    # data in (x) and data out (y), also with work arrays (work)
    attr_reader :x, :y, :work
    # Array of subconvolutions steps
    attr_reader :subops
    # Variables of the convolution operations
    attr_reader :vars

    def initialize(filter,ndim,bc,options={})
      @filter = filter
      @ndim = ndim
      @bc = bc
      @dims = (0...@ndim).collect{ |indx| 
        BOAST::Int("n#{indx+1}",{:direction => :in, :signed => false})
      }
      @eks = []
      @vars = @dims.dup

      #each of the dimension should be adapted to the bc of the system according to filter lengths
      dimx = @dims.collect{ |dim|
        BOAST::Dim(0, dim-1)
      }
      dimy = @dims.collect{ |dim|
        BOAST::Dim(0, dim-1)
      }
      dimw = @dims.collect{ |dim|
        BOAST::Dim(0, dim-1)
      }
      #shoud put input in the case of convolutions with work array
      @x = BOAST::Real("x",:dir => :in, :dim => dimx, :restrict => true)
      @y = BOAST::Real("y",:dir => :out, :dim => dimy, :restrict => true)
      @vars += [@x, @y]
      @vars.push @work = BOAST::Real("work",:dir => :inout, :dim => dimw) if options[:work]
      @vars.push @alpha = BOAST::Real("scal",:dir => :in,:dim => [ BOAST::Dim(0, @ndim-1)]) if options[:alpha]
      @vars.push @beta = BOAST::Real("c",:dir => :in) if options[:beta]
      @vars.push @eks = BOAST::Real("kstrten",:dir => :out,:dim =>[ BOAST::Dim(0, @ndim-1)]) if options[:eks]
      #save options for future reference
      @options=options
    end
    
    def dim_indexes(processed_dim, transpose)
      raise "Invalid processed_dim: #{processed_dim}" if processed_dim >= @ndim or processed_dim < 0 
      if transpose == 1 then
        return [1, 0]
      elsif transpose == -1
        return [0, 1]
      else
        if processed_dim == 0 then
          return [1, 0]
        elsif processed_dim == @ndim - 1
          return [0, 1]
        else
          return [2,0,1]
        end
      end
    end

    # return BOAST Procedure corresponding to the application of the multidimensional convolution
    def procedure(optimization)
      function_name = @filter.name
      fr = @bc.collect { |e| ( (e != 0 and e != 10) ? "f" : "p") }
      function_name += "_" + fr.join("")
      unroll, unroll_dims = optimization.unroll
      if unroll.max > 0 then
        function_name += "_u#{unroll.join("_")}"
      end

      transpose = optimization.transpose
      processed_dims = optimization.dim_order

      subops_dims = self.convolution_steps_dims(processed_dims,transpose)

      subops_data = self.convolution_steps_data
      #set of 1 chaining of convolutions
      subops= (0...@ndim).collect{ |ind|
        ConvolutionOperator1d::new(@filter, @bc[ind], transpose,
                                   self.dim_indexes(ind,transpose),(ind==0 and @options[:beta]),
                                   @options)
      }
      procs = []
      subops.each_with_index{ |subop,ind|
        procs.push subop.procedure(unroll[ind],unroll_dims[ind],optimization.use_mod[ind])
      }
      p= BOAST::Procedure(function_name,@vars){
        procs.each_with_index{ |sub,ind|
          vars = subops_dims[ind]
          vars += subops_data[ind]
          vars += [@alpha[ind]] if @options[:alpha]
          vars += [@beta] if @options[:beta] and ind==0
          vars += [@eks[ind]] if @options[:eks]
          BOAST::print sub.call(*vars)
        }
      }
      return [ p, procs ]
    end
    def convolution_steps_data()
      ndim = @ndim
      if @work then 
        if ndim%2 then
          out = [ [ @x, @y] ]
          ndim -= 1
        else
          out = [ [ @x, @work], [ @work, @y ] ]
          ndim -= 2
        end
        while ndim > 0 do
          out += [ [ @y, @work], [ @work, @y ] ]
          ndim -= 2
        end
        return out
      else
        out = [ [ @x, @y] ]*ndim
      end
    end
    def convolution_steps_dims(processed_dims,transpose)
      dims_tmp=@dims.dup
      return processed_dims.collect{ |ind| 
        if transpose != 0 then
          dims_tmp2 = dims_tmp.dup
          n = dims_tmp2.delete_at(ind)
          ndat = eval '"#{dims_tmp2.join("*")}"'
          if transpose == 1 then
            res = [n,ndat]
          else
            res = [ndat,n]
          end
        else
          n = dims_tmp[ind]
          ndat1 = eval '"#{dims_tmp[0...ind].join("*")}"'
          ndat2 = eval '"#{dims_tmp[ind+1..-1].join("*")}"'

          res=[]
          res += [ndat1] if ndat1 != ""
          res += [n]
          res += [ndat2] if ndat2 != ""
        end
        #here we should apply the correction to the dimension depending on bc
        #count where we are and which dimensions have alredy been treated
        bc=@bc[ind]
        dims_tmp[ind]=n
        #end of correction
        res
      }
    end


  end

  class ConvolutionOptimization
    # apply transposition paradigm or not in different directions: +1,0,-1
    attr_reader :transpose
    # order of the treated dimensions
    attr_reader :dim_order
    # 2-uple containing the unrolling lengths and the unrolling dims
    attr_reader :unroll
    # use the mod_arr strategy for the convolutions
    attr_reader :use_mod
    def initialize(convolution,options)

      ndim = convolution.dims.length

      @transpose = 0
      @transpose = options[:transpose] if options[:transpose]
      
      @use_mod = convolution.dims.collect { false }
      convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc) } if options[:use_mod]

      @dim_order=(0...ndim).collect{|i| i}
      @dim_order.reverse!  if @transpose == -1 
      @dim_order = options[:dim_order] if options[:dim_order] and @transpose == 0

      if @transpose ==1 then
        unrolled_dim = [ 1 ] * ndim
      elsif @transpose == -1 then
        unrolled_dim = [ 0 ] * ndim
      else
        hdim = ndim / 2
        unrolled_dim = (1..ndim).collect { |i| i < hdim ? 2 : 0 }
        unrolled_dim = [1] + unrolled_dim 
        if options[:unrolled_dim] then
          raise 'Incoherent dimensions specified in unrolling options: #{ndim}, #{options[:unrolled_dim].length}' if options[:unrolled_dim].length != ndim
          unrolled_dim = options[:unrolled_dim]
        end
      end
      if options[:unroll] then
        unro = [options[:unroll]].flatten
        if unro.length == 1 then
          unroll_lengths = unro * ndim
        elsif unro.length == ndim then
          unroll_lengths = unro
        else
          raise 'Incoherent dimensions specified in unrolling options: #{ndim}, #{unro.length}'
        end
      else
          unroll_lengths = [1] * ndim
      end
      @unroll = [ unroll_lengths, unrolled_dim ]
    end
  end


end


