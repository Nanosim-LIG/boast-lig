module BOAST

  class ConvolutionFilter
    # List of the floating point values of the convolution filter
    attr_reader :fil_array
    # central point of the filter
    attr_reader :center
    # BOAST variables
    # Filter array (to be used on BOAST functions)
    attr_reader :fil
    # extremes of the filter, calculated via its central point (integers)
    attr_reader :lowfil_val, :upfil_val
    # extremes of the filter, calculated via its central point (BOAST object)
    attr_reader :lowfil, :upfil
    # name of the filter
    attr_reader :name

    def initialize(name,filt,center)
      arr = ConstArray::new(filt,Real)
      @fil=BOAST::Real("fil",:constant => arr,:dim => [ BOAST::Dim((-center),(filt.length - center -1)) ])
      @lowfil_val = -center
      @upfil_val = filt.length - center -1
      @lowfil = BOAST::Int("lowfil",:constant => @lowfil_val)
      @upfil = BOAST::Int("upfil",:constant => @upfil_val)
      @fil_array=filt
      @center=center
      @name=name
    end

  end #class ConvolutionFilter

  # determine the BC of the convolutions
  #        Typical values are 0 : periodic BC, the size of input and output arrays are identical
  #                           1 : Free BC, grow: the size of the output array is equal to the 
  #                                              size of input array plus the one of the filter
  #                          -1 : Free BC, shrink: the size of the input array is equal to the one
  #                                        of output array plus the filter
  #                     Given a convolution and its inverse, in free BC -1 is the inverse of 1
  #                      but not viceversa as the number of point treated is lower.
  #                         10:  Free BC, the size of input and output arrays are identical
  #                              (loss of information)
  class BoundaryConditions
    # conditions names
    PERIODIC = 0
    GROW = 1
    SHRINK = -1
    # determine if the boundary condition is free or not
    attr_reader :free 
    # name of the boundary condition, used for the name of the routines
    attr_reader :name
    # determine if the convolution is a grow-like type (cannot be true if shrink is true)
    attr_reader :grow
    # determine if the convolution is a shrink-like type (cannot be true if grow is true)
    attr_reader :shrink
    # original id of the boundary conditions
    attr_reader :id
    def initialize(ibc)
      @id     = ibc
      @free   = (ibc != PERIODIC)
      @grow   = (ibc == GROW)
      @shrink = (ibc == SHRINK)
      if not @free then
        @name = 'p'
      else
        @name = 'f'
        if @grow then
          @name += 'g'
        elsif @shrink then
          @name += 's'
        end
      end
    end
  end #class BoundaryConditions

  BC = BoundaryConditions

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
      #growed dimension, to be used either for extremes or for mod_arr
      @dim_ngs = BOAST::Dim(@filter.lowfil,@dim_n+@filter.upfil-1)
      #dimensions corresponding to the output of a grow operation
      dim_nsg = BOAST::Dim(-@filter.upfil,@dim_n-@filter.lowfil-1)
      @dims = [@dim_n]
      if (dim_indexes.length == 3) then
        @dims = [BOAST::Int("ndat1",:dir =>:in)] + @dims + [BOAST::Int("ndat2",:dir =>:in)]
        #shrinked and growed dimensions to be used for non-periodic bcs
        dimss = [ BOAST::Dim(0, @dims[0] -1) , @dim_ngs, BOAST::Dim(0, @dims[-1] -1) ] 
        dimsg = [ BOAST::Dim(0, @dims[0] -1) , dim_nsg, BOAST::Dim(0, @dims[-1] -1) ] 
      elsif dim_indexes.last == 0
        @dims = @dims + [BOAST::Int("ndat",:dir =>:in)]
        dimss = [ @dim_ngs, BOAST::Dim(0, @dims[-1] -1) ] 
        dimsg = [ dim_nsg, BOAST::Dim(0, @dims[-1] -1) ] 
      else
        @dims = [BOAST::Int("ndat",:dir =>:in)] + @dims
        dimss = [ BOAST::Dim(0, @dims[0] -1) , @dim_ngs ] 
        dimsg = [ BOAST::Dim(0, @dims[0] -1) , dim_nsg ] 
      end
      @vars = @dims.dup
      # dimension of the problem
      dim_data = @dims.collect{ |dim|
        BOAST::Dim(0, dim-1)
      } 
      if bc.shrink then
        dimx = dimss
      else
        dimx = dim_data
      end
      if bc.grow then
        dimy = dimsg
      else
        dimy = dim_data
      end
      if transpose !=0  then
        dimy = dimy.rotate(1)
      end
      @vars.push @in = BOAST::Real("x",:dir => :in, :dim => dimx, :restrict => true)
      @vars.push @out = BOAST::Real("y",:dir => :out, :dim => dimy, :restrict => true)
      @vars.push @alpha = BOAST::Real("alpha",:dir => :in) if options[:alpha]
      @vars.push @beta = BOAST::Real("beta",:dir => :in) if options[:beta] and init
      @vars.push @dotp = BOAST::Real("dotp",:dir => :out) if options[:dotp]
      @init = init
      @options = options
    end

    def procedure(unroll,unrolled_dim,use_mod,tt_arr)
      function_name = @filter.name + "_" + @bc.name + "_#{@dim_indexes.join('')}" +
        "_u#{unroll}_#{unrolled_dim}_#{use_mod}"

      l = BOAST::Int("l")
      #try to modify tt scalars into arrays of size unroll
      if tt_arr then
        tt = BOAST::Real("tt", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true)
      else
        tt = (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("tt#{index}") }
      end
      iters =  (1..@dims.length).collect{ |index| BOAST::Int("i#{index}")}
      
      if use_mod then
        # the mod_arr behaves as a shrink operation
        #mods=BOAST::Int("mod_arr", :allocate => true, :dim => [@dim_ngs])
        mods=BOAST::Int("mod_arr", :allocate => true, 
                        :dim => [BOAST::Dim(@filter.lowfil - @filter.upfil, @filter.upfil - @filter.lowfil - 1)])
      else
        mods=nil
      end
      return BOAST::Procedure(function_name,vars,[@filter.lowfil,@filter.upfil]){
        BOAST::decl @filter.fil
        BOAST::decl *iters
        BOAST::decl l
        if tt_arr then
          BOAST::decl tt
        else
          BOAST::decl *tt
        end
        if use_mod then
          BOAST::decl mods 
          #BOAST::print BOAST::For(l, @filter.lowfil, @dim_n -1 + @filter.upfil) {
          BOAST::print BOAST::For(l, @filter.lowfil - @filter.upfil, @filter.upfil - @filter.lowfil - 1) {
            BOAST::print mods[l] === BOAST::modulo(l, @dim_n)
          }
        end
        BOAST::print @dotp === 0.0 if @options[:dotp]
        if BOAST::get_lang == BOAST::FORTRAN then
          BOAST::get_output.print("!$omp parallel default(shared)&\n")
          BOAST::get_output.print("!$omp reduction(+:#{dotp})&\n") if @options[:dotp]
          BOAST::get_output.print("!$omp private(#{iters.join(",")},#{l})&\n")
          if tt_arr then
            BOAST::get_output.print("!$omp private(#{tt})\n")
          else
            BOAST::get_output.print("!$omp private(#{tt.join(",")})\n")
          end
        elsif BOAST::get_lang == BOAST::C then
          if tt_arr then
            BOAST::get_output.print("#pragma omp parallel default(shared) #{@options[:dotp] ? "reduction(+:#{dotp})" : ""} private(#{iters.join(",")},#{l},#{tt})\n{\n")
          else
            BOAST::get_output.print("#pragma omp parallel default(shared) #{@options[:dotp] ? "reduction(+:#{dotp})" : ""} private(#{iters.join(",")},#{l},#{tt.join(",")})\n{\n")
          end
        end

        convolution1d(iters,l,tt,mods,unrolled_dim,unroll)

        BOAST::get_output.print("!$omp end parallel\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("}\n")  if BOAST::get_lang == BOAST::C
      }
    end

    #here follows the internal operations for the convolution 1d
    def convolution1d(iters,l,t,mods,unro,unrolling_length)
      convgen= lambda { |dim,t,tlen,reliq|
        ises0 = startendpoints(@dims[dim[0]],unro == dim[0],unrolling_length,reliq)
        BOAST::get_output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("#pragma omp for\n") if BOAST::get_lang == BOAST::C
        For::new(iters[dim[0]], *ises0 ) {
          if dim.length == 3 then
            ises1 = startendpoints(@dims[dim[1]],unro == dim[1],unrolling_length,reliq)
            For::new(iters[dim[1]], *ises1) {
              conv_lines(iters,l,t,tlen,dim[-1],unro,mods)
            }.print
          else
            conv_lines(iters,l,t,tlen,dim[-1],unro,mods)
          end
        }.print
        BOAST::get_output.print("!$omp end do\n") if BOAST::get_lang == BOAST::FORTRAN
      }
      #first without the reliq
      convgen.call(@dim_indexes,t,unrolling_length,false)
      #then with the reliq but only if the unrolling patterns need it
      convgen.call(@dim_indexes,t,1,true) if (unrolling_length > 1)
    end

    #returns the starting and ending points of the convolutions according to unrolling and unrolled dimension
    def startendpoints(dim,unroll,unrolling_length,in_reliq)
      istart= (in_reliq and unroll) ? (dim/unrolling_length)*unrolling_length : 0
      iend  = (unroll and not in_reliq) ? dim-unrolling_length : dim-1
      istep = (unroll and not in_reliq) ? unrolling_length : 1 
      return [istart,iend,istep]
    end

    def conv_lines(iters,l,t,tlen,processed_dim,unro,mods)
      # the only difference holds for the grow operation
      if @bc.grow then
        line_start = -@filter.upfil
        line_end = @dims[processed_dim]-1-@filter.lowfil
      else
        #for the periodic or constant operation the loop is between 0 and n-1
        line_start = 0
        line_end = @dims[processed_dim]-1
      end
      # the shrink operation contains the central part only
      if @bc.shrink then
        BOAST::For(iters[processed_dim], line_start, line_end) {
          for_conv(0,iters,l,t,tlen,processed_dim, unro,mods)
        }.print
      else
        BOAST::For(iters[processed_dim], line_start, -@filter.lowfil-1) {
          for_conv(-1,iters,l,t,tlen,processed_dim, unro,mods)
        }.print
        BOAST::For(iters[processed_dim], -@filter.lowfil, @dims[processed_dim]-1 - @filter.upfil) {
          for_conv(0,iters,l,t,tlen,processed_dim, unro,mods)
        }.print
        BOAST::For(iters[processed_dim], @dims[processed_dim] - @filter.upfil, line_end) {
          for_conv(1,iters,l,t,tlen,processed_dim, unro,mods)
        }.print
      end
    end

    def for_conv(nobc,i_in,l,t,tlen,processed_dim,unro,mods)
      #t.each_index { |ind|
      (0...tlen).each{ |ind|
        i_out = output_index(unro, i_in, ind)
        #WARNING: the eks conditional here can be relaxed
        BOAST::print t[ind] === ((@init and not @dotp) ? @beta * @in[*i_out] / @alpha : 0.0)
      }
      if ( @bc.free and nobc != 0) then
        loop_start = BOAST::max(-i_in[processed_dim], @filter.lowfil)
        loop_end   = BOAST::min(@filter.upfil, @dims[processed_dim] - 1 - i_in[processed_dim])
      elsif
        loop_start=@filter.lowfil
        loop_end=@filter.upfil
      end
      BOAST::For( l,loop_start,loop_end) {
        (0...tlen).each{ |ind|
         #t.each_index { |ind|
          if @bc.free or (nobc == 0) then
            i_out = output_index(unro, i_in, ind,processed_dim,l)
          elsif mods then
            i_out = output_index(unro, i_in, ind,processed_dim,l,nil,mods,nobc) 
          else
            i_out = output_index(unro, i_in, ind,processed_dim,l,@dims[processed_dim])
          end
          BOAST::print t[ind] === t[ind] + @in[*i_out]*@filter.fil[l]
        }
      }.unroll#

      #t.each_index { |ind|
      (0...tlen).each{ |ind|
        i_out = output_index(unro, i_in, ind)
        BOAST::print t[ind] === t[ind] * @alpha if @alpha
        BOAST::print @dotp === @dotp + t[ind] * @in[*i_out] if @dotp
        if @transpose != 0 then
          BOAST::print @out[*i_out.rotate(@transpose)] === t[ind]
        else
          BOAST::print @out[*i_out] === 
            (@init ? t[ind] : @out[*i_out] + t[ind] )
        end
          
      }
    end

    #returns the indices of the output array according to the starting point in the input and of the
    ## processed dimension as well as th eposition in the convolution
    def output_index(unrolling_dim, i_in,unroll_index,processed_dim=nil,lconv_index=nil,
                          ndim_processed=nil,wrapping_array=nil,side=nil)
      if ndim_processed then
        i_out=output_index_k_mod(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,ndim_processed)
      elsif wrapping_array then
        i_out=output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array,side)
      elsif processed_dim then
        i_out=output_index_k(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index)
      else
        i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      end
      return i_out
    end

    #returns the indices of the output according to which of the directions is unrolled
    def output_index_unroll(unrolling_dim, i_in,unroll_index)
      i_out=(0...i_in.length).collect { |indx| unrolling_dim == indx ? i_in[unrolling_dim] + (unroll_index) : i_in[indx]}
      return i_out
    end
    # index of the convolution in the internal region, k=i+l, otherwise the index in the unrolling dimension
    def output_index_k(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      return (0...i_in.length).collect { |indx| processed_dim == indx ? lconv_index +i_in[processed_dim] : i_out[indx]}
    end
    # index in the external region wrapped around (periodic BC), thanks to the presence of the wrapping_array
    # if the side is the left one (-1), the recipe is the usual, otherwise (+1) the recipe is subtracted
    # the the value of the processed dim. In this way the size of mod_arr is only dependent by the size of the 
    # filter which makes life easier for the compiler
    def output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array,side)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      if (side == 1) then
        return (0...i_in.length).collect { |indx| 
          processed_dim == indx ? wrapping_array[lconv_index +i_in[processed_dim] - @dims[processed_dim]] : i_out[indx]}
      else
        return (0...i_in.length).collect { |indx| processed_dim == indx ? wrapping_array[lconv_index +i_in[processed_dim]] : i_out[indx]}
      end
    end
    # index in the external region wrapped around (periodic BC), where the wrapping is given by the integer division
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
      @bc = bc.collect{ |ibc| BoundaryConditions::new(ibc)}
      @dims = (0...@ndim).collect{ |indx| 
        BOAST::Int("n#{indx+1}",{:direction => :in, :signed => false})
      }
      @eks = []
      @vars = @dims.dup

      #each of the dimension should be adapted to the bc of the system according to filter lengths
      dimx = []
      dimy = []
      dimw = []
      #the sizes of the work array are the maximum between inpu and output
      #but this can be reduced by considering the order of the subops
      #however this order is only fixed in the optimization part
      #therefore the value of the array work might be updated
      #before calling the procedure
      @dims.each_with_index { |dim,ind|
        dimin, dimout, dimwork = dim_from_bc(dim,@bc[ind])
        dimx.push(dimin)
        dimy.push(dimout)
        dimw.push(dimwork)
      }
      #dimx = @dims.collect{ |dim|
      #  BOAST::Dim(0, dim-1)
      #}
      #dimy = @dims.collect{ |dim|
      #  BOAST::Dim(0, dim-1)
      #}
      #dimw = @dims.collect{ |dim|
      #  BOAST::Dim(0, dim-1)
      #}

      #should put input in the case of convolutions with work array
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

    #define the boast dimensions of the input, output and work array as a function of the 1d bc
    def dim_from_bc(dim,bc)
      dim_data = BOAST::Dim(0, dim-1)
      if bc.grow then
        dim_nsg = BOAST::Dim(-@filter.upfil,dim - 1 - @filter.lowfil) 
        dim_arrays = [ dim_data, dim_nsg , dim_nsg]
      elsif bc.shrink then
        dim_ngs = BOAST::Dim(@filter.lowfil,dim - 1 + @filter.upfil)
        dim_arrays = [ dim_ngs, dim_data, dim_ngs ]
      else
        dim_arrays = [ dim_data, dim_data, dim_data ]
      end
      return dim_arrays
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
      function_name = @filter.name + "_" +
        @bc.collect { |e| e.name }.join("")
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
                                   self.dim_indexes(ind,transpose),
                                   (ind==0 and @options[:beta]),
                                   @options)
      }
      procs = []
      subops.each_with_index{ |subop,ind|
        procs.push subop.procedure(unroll[ind],unroll_dims[ind],optimization.use_mod[ind],optimization.tt_arr[ind])
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
    #this procedure returns the values of the dimensions of the 
    #data problem that has to be treated
    #the values of the ndat has to be coherent with the
    #boundary conditions 
    def convolution_steps_dims(processed_dims,transpose)
      # first copy the data dimension as reference
      dims_data = []
      # then consider the starting point of the data
      dims_actual = []
      @dims.each_with_index { |d,ind|
        dims_data.push(d)
        if @bc[ind].shrink then
          dims_actual.push(d-@filter.lowfil+@filter.upfil)
        else
          dims_actual.push(d)
        end
      }
      return processed_dims.collect{ |ind| 
        n = dims_data[ind]
        if transpose != 0 then
          dims_tmp2 = dims_actual.dup
          dims_tmp2.delete_at(ind)
          ndat = eval '"(#{dims_tmp2.join(")*(")})"'
          if transpose == 1 then
            res = [n,ndat]
          else
            res = [ndat,n]
          end
        else
          #n = dims_data[ind]
          ndat1 = eval '"(#{dims_actual[0...ind].join(")*(")})"'
          ndat2 = eval '"(#{dims_actual[ind+1..-1].join(")*(")})"'

          res=[]
          res += [ndat1] if ndat1 != "()"
          res += [n]
          res += [ndat2] if ndat2 != "()"
        end
        #here we should apply the correction to the dimension depending on bc
        #count where we are and which dimensions have already been treated
        if @bc[ind].grow then
          dims_actual[ind]=n+@filter.upfil-@filter.lowfil
        else
          dims_actual[ind]=n
        end
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
    # use the tt_arr strategy for the convolutions (temporary variables are arrays and not scalars)
    attr_reader :tt_arr
    def initialize(convolution,options)

      ndim = convolution.dims.length

      @transpose = 0
      @transpose = options[:transpose] if options[:transpose]
      
      @use_mod = convolution.dims.collect { false }
      convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc.free) } if options[:use_mod]

      @tt_arr = convolution.dims.collect { false }
      if options[:tt_arr] then
        ttopt=[options[:tt_arr]].flatten
        if ttopt.length == 1 then
          @tt_arr = convolution.dims.collect { ttopt[0] }
        elsif ttopt.length == ndim then
          @tt_arr = ttopt
        else
          raise 'Incoherent dimensions specified in tt_arr options: #{ndim}, #{ttopt.length}'
        end
      end
        
      convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc.free) } 

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


