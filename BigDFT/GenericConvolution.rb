module BOAST

  class ConvolutionFilter
    # List of the floating point values of the convolution filter
    attr_reader :fil_array
    # central point of the filter
    attr_reader :center
    attr_reader :length
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
      @upfil_val = filt.length - center - 1
      @lowfil = BOAST::Int("lowfil",:constant => @lowfil_val)
      @upfil = BOAST::Int("upfil",:constant => @upfil_val)
      @fil_array = filt
      @center = center
      @name = name
      @length = filt.length
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


  #handle the choice of the best kernels in a given architecture
  class GenericOptimization

    attr_reader :repeat
    attr_reader :dimensions
    
    class DataRange
      def initialize(start,stop,step)
        @range = start..stop
        @step = step
      end
      def each(&block)
        return @range.step(@step,&block) 
      end
    end

    class ParamSpace
      attr_accessor :space
      def initialize(space={})
        @space=space
      end
      def size
        return @space.size
      end
      def points
        pts=[]
        space2 = @space.dup
        dimension,value = space2.shift 
        space2=ParamSpace::new(space2)
        value.each{ |val| 
          pts.push({dimension.to_s.chomp("_range").to_sym => val})
        }
        if space2.size == 0 then
          return pts
        else
          pts3=[]
          pts.each{ |p1| 
            space2.each { |p2| 
              pts3.push(p1.dup.update(p2))
            }
          }
          return pts3
        end
      end
      def each(&block)
        return self.points.each(&block)
      end
    end

#    def initialize(unroll_range=1,mod_arr_test=true,tt_arr_test=true,
#                   unrolled_dim_index_test=false)
    def initialize(options = {})
      unroll_range = 1
      unroll_range = options[:unroll_range] if options[:unroll_range]
      mod_arr_test = true
      mod_arr_test = options[:mod_arr_test] if not options[:mod_arr_test].nil?
      tt_arr_test = false
      tt_arr_test = options[:tt_arr_test] if not options[:tt_arr_test].nil?
      unrolled_dim_index_test = false
      unrolled_dim_index_test = options[:unrolled_dim_index_test] if not options[:unrolled_dim_index_test].nil?
      unroll_inner_test = false
      unroll_inner_test = options[:unroll_inner_test] if not options[:unroll_inner_test].nil?
      @repeat = 3
      @repeat = options[:repeat] if options[:repeat]
      @dimensions = [124,132,130]
      @dimensions = [options[:dimensions]].flatten if options[:dimensions]

      unrl_rng=[unroll_range].flatten
      if unrl_rng.length == 2 then
        unrl_rng=[*unrl_rng,1]
      elsif unrl_rng.length == 1 then
        unrl_rng=[1,*unrl_rng,1]
      end
      space={}
      space[:unroll_range] = DataRange::new(*unrl_rng[0..2])
      space[:mod_arr_range] = [true]
      space[:mod_arr_range] = [true,false] if mod_arr_test
      space[:tt_arr_range] = [false]
      space[:tt_arr_range] = [true,false] if tt_arr_test
      space[:unrolled_dim_index_range] = [0]
      space[:unrolled_dim_index_range] = [0,1] if unrolled_dim_index_test
      space[:unroll_inner_range] = [true]
      space[:unroll_inner_range] = [true, false] if unroll_inner_test
      @space=ParamSpace::new(space)
    end
    
    def each(&block)
      return @space.each(&block)
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
    attr_reader :base_name
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
    # * +:ld+ - leading dimensions enable
    def initialize(filter,bc,transpose,dim_indexes,init,options={})
      @filter = filter
      @bc = bc
      @transpose = transpose
      @dim_indexes = dim_indexes
      @ld = options[:ld]
      @dim_n = BOAST::Int("n",:dir =>:in)
      if @ld then
        @ld_in = BOAST::Int("ld_in",:dir =>:in)
        @ld_out = BOAST::Int("ld_out",:dir =>:in)
      else
        @ld_in = @dim_n
        @ld_out = @dim_n
      end
      @dims = [@dim_n]
      @dims_in = [@ld_in]
      @dims_out = [@ld_out]
      if (dim_indexes.length == 3) then
        ndat1 = BOAST::Int("ndat1",:dir =>:in)
        ndat2 = BOAST::Int("ndat2",:dir =>:in)
        @dims     = [ndat1] + @dims     + [ndat2]
        @dims_in  = [ndat1] + @dims_in  + [ndat2]
        @dims_out = [ndat1] + @dims_out + [ndat2]
      elsif dim_indexes.last == 0
        ndat = BOAST::Int("ndat",:dir =>:in)
        @dims     = @dims     + [ndat]
        @dims_in  = @dims_in  + [ndat]
        @dims_out = @dims_out + [ndat]
      else
        ndat = BOAST::Int("ndat",:dir =>:in)
        @dims     = [ndat] + @dims
        @dims_in  = [ndat] + @dims_in
        @dims_out = [ndat] + @dims_out
      end
      @vars = @dims.dup
      if @ld then
        @vars.push( @ld_in )
        @vars.push( @ld_out )
      end
      #growed dimension, to be used either for extremes or for mod_arr
      @dim_ngs = BOAST::Dim( @filter.lowfil, @dim_n + @filter.upfil  - 1)
      #dimensions corresponding to the output of a grow operation
      @dim_nsg = BOAST::Dim(-@filter.upfil,  @dim_n - @filter.lowfil - 1)
      # dimension of the problem
      dim_data = @dims.collect{ |dim|
        BOAST::Dim(0, dim - 1)
      }
      dimx = @dims_in.collect{ |dim|
        if not dim.name.match("ndat") and bc.shrink and not @ld then
          @dim_ngs
        elsif not dim.name.match("ndat") and bc.shrink
          BOAST::Dim(@filter.lowfil, dim + @filter.lowfil - 1)
        else
          BOAST::Dim(0, dim - 1)
        end
      }
      dimy = @dims_out.collect{ |dim|
        if not dim.name.match("ndat") and bc.grow and not @ld then
          @dim_nsg
        elsif not dim.name.match("ndat") and bc.grow then
          BOAST::Dim( -@filter.upfil, dim - @filter.lowfil - 1)
        else
          BOAST::Dim(0, dim - 1)
        end
      }
      if transpose !=0  then
        dimy = dimy.rotate(transpose)
      end
      @vars.push @in = BOAST::Real("x",:dir => :in, :dim => dimx, :restrict => true)
      @vars.push @out = BOAST::Real("y",:dir => :out, :dim => dimy, :restrict => true)
      @vars.push @alpha = BOAST::Real("alpha",:dir => :in) if options[:alpha]
      @vars.push @beta = BOAST::Real("beta",:dir => :in) if options[:beta] and init
      @vars.push @dotp = BOAST::Real("dotp",:dir => :out) if options[:dotp]
      @init = init
      @accumulate = false
      @accumulate = options[:accumulate] if options[:accumulate]
      @options = options
      @base_name = @filter.name + "_" + @bc.name + "_#{@dim_indexes.join('')}"
      @base_name += "_alpha" if @alpha
      @base_name += "_beta" if @beta
      @base_name += "_dotp" if @dotp
      @base_name += "_acc" if @accumulate
      @base_name += "_ld" if @ld
    end

    def params(dim, index=0)
      vars=[]
      varsin=[]
      varsout=[]
      nd = { "n" => dim[index], "ndat" => 1, "ndat1" => 1, "ndat2" => 1 }
      if @dims.length==2 then
        dim.each_index { |indx|
          nd["ndat"] *= dim[indx] if indx != index
        }
      else
        dim.each_index { |indx|
          nd["ndat1"] *= dim[indx] if indx < index
          nd["ndat2"] *= dim[indx] if indx > index
        }
      end
      @dims.each { |dim|

        vars.push(nd[dim.name])
        if dim.name == "n" then
          if @bc.grow then
            varsin.push(nd["n"])
            varsout.push(nd["n"] + @filter.length - 1)
          elsif @bc.shrink
            varsin.push(nd["n"] + @filter.length - 1)
            varsout.push(nd["n"])
          else
            varsin.push(nd["n"])
            varsout.push(nd["n"])
          end
        else
          varsin.push(nd[dim.name])
          varsout.push(nd[dim.name])
        end
      }
      varsout.rotate!(@transpose)
      #input and output arrays
      if @ld then
        if @bc.grow then
          vars.push(nd["n"])
          vars.push(nd["n"] + @filter.length - 1)
        elsif @bc.shrink
          vars.push(nd["n"] + @filter.length - 1)
          vars.push(nd["n"])
        else
          vars.push(nd["n"])
          vars.push(nd["n"])
        end
      end
      vars.push(NArray.float(*varsin).random)
      vars.push(NArray.float(*varsout))
      #accessory scalars
      nscal=0
      nscal+=1 if @alpha
      nscal+=1 if @beta
      nscal.times{vars.push(0.5)}
      vars.push(NArray.float(1).random) if @dotp
      return vars
    end

    def optimize(opt_space)
      opt_space=GenericOptimization::new if not opt_space
      t_best=Float::INFINITY
      p_best = nil
      opt_space.each{ |optim|
        next if optim[:unrolled_dim_index] == 1 and @dims.length < 3
        next if optim[:mod_arr] and @bc.free
        #puts optim
        kernel = CKernel::new
        BOAST::set_output( kernel.code )
        kernel.lang = BOAST::get_lang
        if BOAST::get_lang == C then
          BOAST::get_output.print "inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
          BOAST::get_output.print "inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
          BOAST::get_output.print "inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
        end

        p = self.procedure(optim)
        BOAST::print p
        kernel.procedure = p
        kernel.build(:openmp => true)
        #kernel.print if @bc.free
        t_mean = 0
        dimensions = opt_space.dimensions
        par = nil
        if dimensions.length < @dims.length then
          dimensions += [dimensions[0]]*(@dims.length-dimensions.length)
        end
	dimensions.length.times { |indx|
          stats_a = []
          par = self.params(dimensions, indx)
          #puts par.inspect
          opt_space.repeat.times {
            stats_a.push kernel.run(*par)
          }
          stats_a.sort_by! { |a| a[:duration] }
          stats = stats_a.first
          if BOAST::get_verbose then
            puts "#{indx} - [#{par[0...@dims.length].join(", ")}] - #{kernel.procedure.name}: #{stats[:duration]*1.0e3} us #{self.cost(*par[0...@dims.length]) / (stats[:duration]*1.0e9)} GFlops"
            puts optim
          end
          t_mean += stats[:duration]
        }
        t_mean /= dimensions.length
        puts "#{kernel.procedure.name}: #{t_mean*1.0e3} us #{self.cost(*par[0...@dims.length]) / (t_mean*1.0e9)} GFlops"
        if t_best > t_mean then
          t_best = t_mean
          p_best = p
        end
      }
      return p_best
    end

    def cost(*dimens)
      n = dimens[@dim_indexes.last]
      ndat = 1
      @dim_indexes[0...-1].each { |indx|
        ndat *= dimens[indx]
      }
      return n * ( 2 * @filter.length ) * ndat
    end

    def procedure(options={})
      #(unroll, unrolled_dim, use_mod, tt_arr)
      #default values
      unroll = 1
      mod_arr = true
      tt_arr = false
      unroll_inner = true
      unroll = options[:unroll] if options[:unroll]
      mod_arr = options[:mod_arr] if not options[:mod_arr].nil?
      tt_arr = options[:tt_arr] if not options[:tt_arr].nil?
      unroll_inner = options[:unroll_inner] if not options[:unroll_inner].nil?

      unrolled_dim=@dim_indexes[0]
      unrolled_dim=@dim_indexes[options[:unrolled_dim_index]] if @dim_indexes.length > 2 and options[:unrolled_dim_index]
     
      mod_arr = false if @bc.free

      function_name = @base_name + 
        "_u#{unroll}_#{unrolled_dim}_#{mod_arr}_#{tt_arr}_#{unroll_inner}"

      l = BOAST::Int("l")
      #try to modify tt scalars into arrays of size unroll
      if tt_arr then
        tt = BOAST::Real("tt", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true)
      else
        tt = (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("tt#{index}") }
      end
      iters =  (1..@dims.length).collect{ |index| BOAST::Int("i#{index}")}
      
      if mod_arr then
        # the mod_arr behaves as a shrink operation
        #mods=BOAST::Int("mod_arr", :allocate => true, :dim => [@dim_ngs])
        mods=BOAST::Int("mod_arr", :allocate => true, 
                        :dim => [BOAST::Dim(@filter.lowfil - @filter.upfil, @filter.upfil - @filter.lowfil - 1)])
      else
        mods=nil
      end
      return BOAST::Procedure(function_name, vars, [@filter.lowfil, @filter.upfil] ){
        BOAST::decl @filter.fil
        BOAST::decl *iters
        BOAST::decl l
        if tt_arr then
          BOAST::decl tt
        else
          BOAST::decl *tt
        end
        if mod_arr then
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

        convolution1d(iters, l, tt, mods, unrolled_dim, unroll, unroll_inner)

        BOAST::get_output.print("!$omp end parallel\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("}\n")  if BOAST::get_lang == BOAST::C
      }
    end

    #here follows the internal operations for the convolution 1d
    def convolution1d(iters, l, t, mods, unro, unrolling_length, unroll_inner)
      convgen= lambda { |dim,t,tlen,reliq|
        ises0 = startendpoints(@dims[dim[0]], unro == dim[0], unrolling_length, reliq)
        BOAST::get_output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
        BOAST::get_output.print("#pragma omp for\n") if BOAST::get_lang == BOAST::C
        For::new(iters[dim[0]], *ises0 ) {
          if dim.length == 3 then
            ises1 = startendpoints(@dims[dim[1]], unro == dim[1], unrolling_length, reliq)
            For::new(iters[dim[1]], *ises1) {
              conv_lines(iters, l, t, tlen, dim[-1], unro, mods, unroll_inner)
            }.print
          else
            conv_lines(iters, l, t, tlen, dim[-1], unro, mods, unroll_inner)
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

    def conv_lines(iters, l, t, tlen, processed_dim, unro, mods, unroll_inner)
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
          for_conv(:center, iters, l, t, tlen, processed_dim, unro, mods, unroll_inner)
        }.print
      else
        BOAST::For(iters[processed_dim], line_start, -@filter.lowfil-1) {
          for_conv(:begin, iters, l, t, tlen, processed_dim, unro, mods, unroll_inner)
        }.print
        BOAST::For(iters[processed_dim], -@filter.lowfil, @dims[processed_dim]-1 - @filter.upfil) {
          for_conv(:center, iters, l, t, tlen, processed_dim, unro, mods, unroll_inner)
        }.print
        BOAST::For(iters[processed_dim], @dims[processed_dim] - @filter.upfil, line_end) {
          for_conv(:end, iters, l, t, tlen, processed_dim, unro, mods, unroll_inner)
        }.print
      end
    end

    def for_conv(side, i_in, l, t, tlen, processed_dim, unro, mods, unroll_inner)
      #t.each_index { |ind|
      (0...tlen).each{ |ind|
        i_out = output_index(unro, i_in, ind)
        #WARNING: the eks conditional here can be relaxed
        BOAST::print t[ind] === ((@init and not @dotp) ? @beta * @in[*i_out] / @alpha : 0.0)
      }
      if ( @bc.free and side != :center) then
        loop_start = BOAST::max(-i_in[processed_dim], @filter.lowfil)
        loop_end   = BOAST::min(@filter.upfil, @dims[processed_dim] - 1 - i_in[processed_dim])
      elsif
        loop_start=@filter.lowfil
        loop_end=@filter.upfil
      end
      f = BOAST::For( l,loop_start,loop_end) {
        (0...tlen).each{ |ind|
         #t.each_index { |ind|
          if @bc.free or (side == :center) then
            i_out = output_index(unro, i_in, ind,processed_dim,l)
          elsif mods then
            i_out = output_index(unro, i_in, ind,processed_dim,l,nil,mods,side) 
          else
            i_out = output_index(unro, i_in, ind,processed_dim,l,@dims[processed_dim])
          end
          BOAST::print t[ind] === t[ind] + @in[*i_out]*@filter.fil[l]
        }
      }
      if unroll_inner then
        f.unroll
      else
        f.print
      end

      #t.each_index { |ind|
      (0...tlen).each{ |ind|
        i_out = output_index(unro, i_in, ind)
        BOAST::print t[ind] === t[ind] * @alpha if @alpha
        BOAST::print @dotp === @dotp + t[ind] * @in[*i_out] if @dotp
        if not @accumulate then
          BOAST::print @out[*i_out.rotate(@transpose)] === t[ind]
        else
          BOAST::print @out[*i_out.rotate(@transpose)] === 
            (@init ? t[ind] : @out[*i_out.rotate(@transpose)] + t[ind] )
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
    # if the side is :begin, the recipe is the usual, otherwise (:end) the recipe is subtracted
    # the the value of the processed dim. In this way the size of mod_arr is only dependent by the size of the 
    # filter which makes life easier for the compiler
    def output_index_k_mod_arr(unrolling_dim, i_in,unroll_index,processed_dim,lconv_index,wrapping_array,side)
      i_out=output_index_unroll(unrolling_dim, i_in,unroll_index)
      if (side == :end) then
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

  class GenericConvolutionOperator
    attr_accessor :procs
    attr_reader :needed_subops
    def initialize(filter,options={})
      @filter = filter
      @ld = options[:ld]
      @m = options[:m]

      @vars = []
      @vars.push @ndim  = BOAST::Int( "d",  :dir => :in )
      @vars.push @dims  = BOAST::Int( "n",  :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      @vars.push @bc    = BOAST::Int( "bc", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      if @ld then
        @vars.push @ld_in  = BOAST::Int( "ld_in", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
        @vars.push @ld_out = BOAST::Int( "ld_out", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      end
      @vars.push @m     = BOAST::Int( "m", :dir => :in ) if @m
      @vars.push @x     = BOAST::Real("x",  :dir => :in,  :restrict => true, :dim => [ BOAST::Dim() ] )
      @vars.push @y     = BOAST::Real("y",  :dir => :out, :restrict => true, :dim => [ BOAST::Dim() ] )
      @vars.push @w1 = BOAST::Real("w1", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
      @vars.push @w2 = BOAST::Real("w2", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
      @vars.push @alpha = BOAST::Real("alpha",:dir => :in,:dim => [ BOAST::Dim(0, @ndim - 1)]) if options[:alpha]
      @vars.push @beta = BOAST::Real("beta",:dir => :in) if options[:beta]
      @vars.push @eks = BOAST::Real("eks",:dir => :out,:dim =>[ BOAST::Dim(0, @ndim - 1)]) if options[:eks]

      @transpose = 0
      @transpose = options[:transpose] if options[:transpose]
      @options=options
      @procs = {}
      @needed_subops = {}
      [BOAST::BC::PERIODIC, BOAST::BC::GROW, BOAST::BC::SHRINK].each { |bc|
        dim_indexes_a = []
        if @transpose == 0 then
          dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
        elsif @transpose == 1
          dim_indexes_a = [ [1, 0] ]
        elsif @transpose == -1
          dim_indexes_a = [ [0, 1] ]
        end
        dim_indexes_a.each{ |dim_indexes|
          p = ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc), @transpose, dim_indexes, false, @options)
          @needed_subops[p.base_name] = p
          @procs[p.base_name] = p.procedure
          if @options[:beta] then
            p = ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc), @transpose, dim_indexes, true, @options)
            @needed_subops[p.base_name] = p
            @procs[p.base_name] = p.procedure
          end
        }
      }
    end

    def optimize(opt_space=nil)
      @needed_subops.each { |name, subop|
        @procs[name] = subop.optimize(opt_space)
      }
    end

    def cost(n, bc, m = 1)
      dims_actual = []
      compute_ni_ndat = lambda { |indx|
        indx = n.length - indx - 1 if @transpose == -1
        ndat = 1
        (0..(n.length - 1)).each { |i|
          ndat *= n[i] if i != indx
        }
        ni_ndat = [ n[indx], ndat ]
        ni_ndat.reverse! if @transpose == -1 or (@transpose == 0 and indx = n.length - 1 )
        indexes = [ 1, 0]
        indexes.reverse! if @transpose == -1 or (@transpose == 0 and indx = n.length - 1 )
        return [ni_ndat, indexes]
      }
      compute_ndat_ni_ndat2 = lambda { |indx|
        ndat = 1
        ndat2 = 1
        (0..(n.length - 1)).each { |i|
          ndat *= n[i] if i < indx
          ndat2 *= n[i] if i > indx
        }
        ni = n[indx]
        ndat_ni_ndat2 = [ndat, ni, ndat2]
        indexes = [2, 0, 1]
        return [ndat_ni_ndat2, indexes]
      }
      (0...n.length).each { |indx|
        if bc[indx] == BOAST::BC::SHRINK then
          dims_actual[indx] = n[indx] + @filter.length - 1
        else
          dims_actual[indx] = n[indx]
        end
      }
      if n.length == 1 then
        d = [ n[0], 1 ]
        d_indexes = [ 1, 0 ]
        dims.reverse! if @transpose == -1
        dim_indexes.reverse! if @transpose == -1
        cost = ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc[0]), @transpose, d_indexes, (true and @options[:beta]), @options).cost( *d )
      else
        cost = 0
        dims, dim_indexes = compute_ni_ndat.call(0)
        cost += ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc[0]), @transpose, dim_indexes, (true and @options[:beta]), @options).cost( *dims )
        dims_left = n.length - 1
        while dims_left > 1 do
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call(n.length-dims_left)
          else
            dims, dim_indexes = compute_ni_ndat.call(n.length-dims_left)
          end
          cost += ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc[0]), @transpose, dim_indexes, false, @options).cost( *dims )
          dims_left -= 1
        end
        dims, dim_indexes = compute_ni_ndat.call(n.length-1)
        cost += ConvolutionOperator1d::new(@filter, BOAST::BC::new(bc[0]), @transpose, dim_indexes, false, @options).cost( *dims )
      end
      return cost * m
    end

    def procedure
      function_name = @filter.name
      function_name += "_ld" if @ld
      function_name += "_m" if @m
      p = BOAST::Procedure(function_name,@vars) {
        dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, @ndim - 1) ] ) if BOAST::get_lang == FORTRAN
        dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, 16) ] ) if BOAST::get_lang == C
        dims_left   = BOAST::Int "dims_left"
        ni = BOAST::Int "ni"
        ndat = BOAST::Int "ndat"
        ndat2 = BOAST::Int "ndat2"
        ndat_tot_in = BOAST::Int "nti" if @m
        ndat_tot_out = BOAST::Int "nto" if @m
        i = BOAST::Int "i"
        j = BOAST::Int "j"
        BOAST::decl i, j, dims_actual, dims_left, ni, ndat, ndat2
        BOAST::decl ndat_tot_in, ndat_tot_out if @m
        dims = []
        dim_indexes = []
        BOAST::print dims_left === @ndim
        BOAST::print BOAST::For( i, 0, @ndim - 1 ) {
          if @ld then
            BOAST::print dims_actual[i] === @ld_in[i]
          else
            BOAST::print BOAST::If(@bc[i] == BOAST::BC::SHRINK, lambda {
              BOAST::print dims_actual[i] === @dims[i] + @filter.length - 1
            }, lambda {
              BOAST::print dims_actual[i] === @dims[i]
            })
          end
        }
        compute_ni_ndat = lambda { |indx|
          indx = @ndim - 1 - indx if @transpose == -1
          BOAST::print ni === @dims[indx]
          BOAST::print ndat === 1
          BOAST::print BOAST::For(j, 0, @ndim - 1) {
            BOAST::print BOAST::If( j != indx ) {
              BOAST::print ndat === ndat * dims_actual[j]
            }
          }
          d = [ ni, ndat ]
          d_indexes = [ 1, 0]
          d.reverse! if @transpose == -1
          d_indexes.reverse! if @transpose == -1
          return [ d, d_indexes ]
        }
        compute_ndat_ni_ndat2 = lambda { |indx|
          BOAST::print ni === @dims[indx]
          BOAST::print ndat === 1
          BOAST::print ndat2 === 1
          BOAST::print BOAST::For(j, 0, @ndim - 1) {
            BOAST::print BOAST::If( j < indx, lambda {
              BOAST::print ndat === ndat * dims_actual[j]
            }, j > indx, lambda {
              BOAST::print ndat2 === ndat2 * dims_actual[j]
            })
          }
          return [ [ndat, ni, ndat2], [2, 0, 1] ]
        }

        print_call = lambda { |indx, init, datas, multi_conv|
          vars = dims
          if multi_conv then
            if dim_indexes.length == 2 then
              BOAST::print ndat_tot_in === dims[dim_indexes[0]]
            else
              BOAST::print ndat_tot_in === dims[dim_indexes[0]] * dims[dim_indexes[1]]
            end
            BOAST::print ndat_tot_out === ndat_tot_in
            BOAST::print ndat_tot_in === ndat_tot_in * dims_actual[indx]
            BOAST::print ndat_tot_out === ndat_tot_out * @ld_out[indx] if @ld
            f = BOAST::For(j, 0, @m-1)
          end
          if @ld then
            vars.push @ld_in[indx]
            vars.push @ld_out[indx]
          end
          indx = @ndim - 1 - indx if @transpose == -1
          vars2 = []
          vars2.push( @alpha[indx] ) if @options[:alpha]
          vars2.push( @beta ) if (init and @options[:beta])
          vars2.push( @eks[indx] ) if @options[:eks]
          dats = []
          BOAST::print BOAST::Case( @bc[indx], BOAST::BC::PERIODIC, lambda {
            procname = ConvolutionOperator1d::new(@filter, BOAST::BC::new(BOAST::BC::PERIODIC), @transpose, dim_indexes, (init and @options[:beta]), @options).base_name

            if multi_conv then
              BOAST::print ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              f.print
            else
              dats = datas.dup
            end
            BOAST::print @procs[procname].call( *vars, *dats, *vars2 )
            f.close if multi_conv
            BOAST::print dims_actual[indx] === @ld_out[indx]  if @ld

          }, BOAST::BC::GROW, lambda {
            procname = ConvolutionOperator1d::new(@filter, BOAST::BC::new(BOAST::BC::GROW),     @transpose, dim_indexes, (init and @options[:beta]), @options).base_name

            if multi_conv then
              BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 1 ) if not @ld
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              f.print
            else
              dats = datas.dup
            end
            BOAST::print @procs[procname].call( *vars, *dats, *vars2 )
            f.close if multi_conv
            if @ld then
              BOAST::print dims_actual[indx] === @ld_out[indx]  if @ld
            else
              BOAST::print dims_actual[indx] === dims_actual[indx] + @filter.length - 1
            end
          }, BOAST::BC::SHRINK, lambda {
            procname = ConvolutionOperator1d::new(@filter, BOAST::BC::new(BOAST::BC::SHRINK),    @transpose, dim_indexes, (init and @options[:beta]), @options).base_name

            if multi_conv then
              BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 1 )  if not @ld
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              f.print
            else
              dats = datas.dup
            end
            BOAST::print @procs[procname].call( *vars, *dats, *vars2 )
            f.close if multi_conv
            if @ld then
              BOAST::print dims_actual[indx] === @ld_out[indx]  if @ld
            else
              BOAST::print dims_actual[indx] === dims_actual[indx] - @filter.length + 1
            end
          })
        }

        BOAST::print BOAST::If( @ndim == 1 , lambda {
          conv_number = 1
          conv_number = @m if @m
          dims = [ @dims[0], conv_number ]
          dim_indexes = [ 1, 0 ]
          dims.reverse! if @transpose == -1
          dim_indexes.reverse! if @transpose == -1
          datas = [ @x, @y ]
          print_call.call( 0, true, datas, false )
        }, lambda {
          dims, dim_indexes = compute_ni_ndat.call( 0 )
          datas = [ @x, @w1 ]
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( 0, true, datas, @m )
          BOAST::print dims_left === dims_left - 1
          BOAST::print i === 1
          BOAST::print BOAST::While( dims_left > 2 ) {
            if @transpose == 0 then
              dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
            else
              dims, dim_indexes = compute_ni_ndat.call( i )
            end
            datas = [ @w1, @w2 ]
            datas = [ @x, @y ] if not @options[:work]
            print_call.call( i, false, datas, @m )
            BOAST::print i === i + 1
            if @transpose == 0 then
              dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
            else
              dims, dim_indexes = compute_ni_ndat.call( i )
            end
            datas = [ @w2, @w1 ]
            datas = [ @x, @y ] if not @options[:work]
            print_call.call( i, false, datas, @m )
            BOAST::print i === i + 1
            BOAST::print dims_left === dims_left - 2
          }
          BOAST::print BOAST::If( dims_left == 2, lambda {
            if @transpose == 0 then
              dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
            else
              dims, dim_indexes = compute_ni_ndat.call( i )
            end
            datas = [ @w1, @w2 ]
            datas = [ @x, @y ] if not @options[:work]
            print_call.call( i, false, datas, @m )
            BOAST::print i === i + 1
            dims, dim_indexes = compute_ni_ndat.call( i )
            if @transpose == 0 then
              dims.reverse!
              dim_indexes.reverse!
            end
            datas = [ @w2, @y ]
            datas = [ @x, @y ] if not @options[:work]
            print_call.call( i, false, datas, @m )
          }, lambda {
            dims, dim_indexes = compute_ni_ndat.call( i )
            if @transpose == 0 then
              dims.reverse! 
              dim_indexes.reverse!
            end
            datas = [ @w1, @y ]
            datas = [ @x, @y ] if not @options[:work]
            print_call.call( i, false, datas, @m )
          })
        })
      }
      return [ p, @procs ]
    end
  end

  class ConvolutionOperator
    # Class describing to the convolution operator. Should contain the filter valus and its central point
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
        BOAST::Int("n#{indx+1}",:direction => :in)
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

      ndim = convolution.ndim

      @transpose = 0
      @transpose = options[:transpose] if options[:transpose]
      
      @use_mod =  [ false ] * ndim
      convolution.bc.each_with_index { |bc,ind| @use_mod[ind] = (not bc.free) } if options[:use_mod]

      @tt_arr = convolution.dims.collect { false }
      if options[:tt_arr] then
        ttopt=[options[:tt_arr]].flatten
        if ttopt.length == 1 then
          @tt_arr = [ttopt[0]] * ndim
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


