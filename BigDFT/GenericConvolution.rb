require 'BOAST'
require 'narray'

BOAST::register_funccall("modulo")
BOAST::register_funccall("min")
BOAST::register_funccall("max")

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

  def initialize(name, filt, center)
    @fil_array = filt.dup
    arr = BOAST::ConstArray::new(@fil_array,BOAST::Real)
    @fil = BOAST::Real("#{name}_fil",:constant => arr,:dim => [ BOAST::Dim((-center),(@fil_array.length - center -1)) ])
    @lowfil_val = -center
    @upfil_val = @fil_array.length - center - 1
    @lowfil = BOAST::Int("lowfil",:constant => @lowfil_val)
    @upfil = BOAST::Int("upfil",:constant => @upfil_val)
    @center = center
    @name = name
    @length = @fil_array.length
  end

end #class ConvolutionFilter

class WaveletFilter
  attr_reader :low, :high
  attr_reader :low_even, :low_odd
  attr_reader :high_even, :high_odd
  attr_reader :low_reverse_even, :low_reverse_odd
  attr_reader :high_reverse_even, :high_reverse_odd
  attr_reader :name
  attr_reader :center
  attr_reader :fil_array
  attr_reader :length
  def initialize(name, filt1)
    @fil_array = filt1.dup
    @center = @fil_array.length/2
    @center -= @center%2
    @low = ConvolutionFilter::new(name+"_l", @fil_array, @center)
    filt2 = []
    @fil_array.each_with_index { |e,i|
      if i % 2 == 0 then
        filt2.push(e)
      else
        filt2.push(-e)
      end
    }
    filt2.reverse!
    @high = ConvolutionFilter::new(name+"_h", filt2, @center)

    center_half = @center / 2

    filt_3 = @fil_array.values_at(*(0..(@fil_array.length-1)).step(2).collect)
    @low_even = ConvolutionFilter::new(name+"_le", filt_3, center_half)

    filt_4 = @fil_array.values_at(*(1..(@fil_array.length-1)).step(2).collect)
    @low_odd = ConvolutionFilter::new(name+"_lo", filt_4, center_half)

    filt_5 = filt2.values_at(*(0..(filt2.length-1)).step(2).collect)
    @high_even = ConvolutionFilter::new(name+"_he", filt_5, center_half)

    filt_6 = filt2.values_at(*(1..(filt2.length-1)).step(2).collect)
    @high_odd = ConvolutionFilter::new(name+"_ho", filt_6, center_half)

    center_half = (filt1.length - @center - 1)/2

    filt_r_3 = @fil_array.reverse.values_at(*(0..(@fil_array.length-1)).step(2).collect)
    @low_reverse_even = ConvolutionFilter::new(name+"_lre", filt_r_3, center_half)

    filt_r_4 = @fil_array.reverse.values_at(*(1..(@fil_array.length-1)).step(2).collect)
    @low_reverse_odd = ConvolutionFilter::new(name+"_lro", filt_r_4, center_half)

    filt_r_5 = filt2.reverse.values_at(*(0..(filt2.length-1)).step(2).collect)
    @high_reverse_even = ConvolutionFilter::new(name+"_hre", filt_r_5, center_half)

    filt_r_6 = filt2.reverse.values_at(*(1..(filt2.length-1)).step(2).collect)
    @high_reverse_odd = ConvolutionFilter::new(name+"_hro", filt_r_6, center_half)

    @length = @fil_array.length
    @name = name
    @default_wavelet=:decompose
  end

  def default_wavelet=(wavelet)
    @default_wavelet = wavelet
  end

  def lowfil(wavelet=nil)
    wavelet = @default_wavelet if not wavelet
    if wavelet == :decompose then
      return @low_even.lowfil
    else
      return @low_reverse_even.lowfil
    end
  end

  def upfil(wavelet=nil)
    wavelet = @default_wavelet if not wavelet
    if wavelet == :decompose then
      return @low_even.upfil
    else
      return @low_reverse_even.upfil
    end
  end
end

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
  FREE = 2
  CONDITIONS = [PERIODIC, GROW, SHRINK]
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
    @discard = (idc == FREE)
    if not @free then
      @name = 'p'
    else
      @name = 'f'
      if @grow then
        @name += 'g'
      elsif @shrink then
        @name += 's'
      elsif @discard then
        @name += 'd'
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
  # output array
  # wavelet:     y <- a wt_fil (x) in + [ a_y * y ]
  # convolution: y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ]
  attr_reader :y
  # reduce the array or not
  attr_reader :reduce
  # constants of the convolution
  attr_reader :a_x, :a, :a_y
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
  # * +:a+ - constant of y <- a wt_fil (x) in + [ a_y * y ] or y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 1
  # * +:a_y+ - constant of y <- a wt_fil (x) in + [ a_y * y ] or y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 0
  # * +:a_x+ - constant of y <- a    fil (x) in + [ a_y * y ] + [ a_x * in ], default 0
  # * +:ld+ - leading dimensions enable
  # * +:wavelet+ - specify a wavelet operation, :decompose or :recompose
  def initialize(filter, bc, transpose, dim_indexes, options={})
    @filter = filter
    @bc = bc
    @transpose = transpose
    @dim_indexes = dim_indexes
    @ld = options[:ld]
    @kinetic = options[:kinetic]
    @wavelet = options[:wavelet]
    if @wavelet then
      @filter.default_wavelet = @wavelet
    end

    self.compute_dims

    @vars = @dims.dup
    if @ld then
      @vars.push( @nx )
      @vars.push( @ny )
    end

    self.compute_inner_loop_boundaries

    dimx, dimy = self.compute_dimx_dimy
    #@init = init

    @vars.push @x = BOAST::Real("x",:dir => :in, :dim => dimx, :restrict => true)
    if @kinetic and @kinetic != :inplace and not @options[:zero_out] then
      @vars.push @x2 = BOAST::Real("x2",:dir => :in, :dim => dimx, :restrict => true)
    end
    @vars.push @y = BOAST::Real("y",:dir => :out, :dim => dimy, :restrict => true)
    if @kinetic and @transpose != 0 then
      @vars.push @y2 =  BOAST::Real("y2", :dir => :out, :dim => dimy, :restrict => true)
    end
    @vars.push @a = BOAST::Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = BOAST::Real("a_x",:dir => :in) if options[:a_x] #and init
    if options[:a_y] then
      if options[:a_y] == 1 then
        @accumulate = true
      else
        @vars.push @a_y = BOAST::Real("a_y",:dir => :in) if options[:a_y]
      end
    end
    @vars.push @dotp = BOAST::Real("dotp",:dir => :out) if options[:dotp]
    @options = options
    @base_name = ""
    if @wavelet then
      if @wavelet == :decompose then
        @base_name += "dwt_"
      else
        @base_name += "idwt_"
      end
    end
    @base_name += @filter.name + "_" + @bc.name + "_#{@dim_indexes.join('')}"
    @base_name += "_ascal" if @a
    @base_name += "_ain" if @a_x
    @base_name += "_ay" if @a_y
    @base_name += "_acc" if @accumulate
    @base_name += "_dotp" if @dotp
    @base_name += "_ld" if @ld
  end

  def compute_dims
    @dim_n = BOAST::Int("n",:dir =>:in)
    if @ld then
      @nx = BOAST::Int("nx",:dir =>:in)
      @ny = BOAST::Int("ny",:dir =>:in)
    else
      @nx = @dim_n
      @ny = @dim_n
    end
    @dims = [@dim_n]
    @dims_in = [@nx]
    @dims_out = [@ny]
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

    if @wavelet then
      if @wavelet == :decompose then
        @dim_ngs = [ BOAST::Dim(0, 1), BOAST::Dim( @filter.lowfil, @dim_n + @filter.upfil - 1 ) ]
        @dim_nsg = [ BOAST::Dim( -@filter.upfil, @dim_n - @filter.lowfil - 1 ), BOAST::Dim(0, 1) ] 
      else
        @dim_ngs = [ BOAST::Dim( @filter.lowfil, @dim_n + @filter.upfil - 1 ), BOAST::Dim(0, 1) ]
        @dim_nsg = [ BOAST::Dim(0, 1), BOAST::Dim( -@filter.upfil, @dim_n - @filter.lowfil - 1 ) ]
      end
    else
      #growed dimension, to be used either for extremes or for mod_arr
      @dim_ngs = BOAST::Dim( @filter.lowfil, @dim_n + @filter.upfil  - 1)
      #dimensions corresponding to the output of a grow operation
      @dim_nsg = BOAST::Dim(-@filter.upfil,  @dim_n - @filter.lowfil - 1)
    end

  end

  def compute_dimx_dimy
    dimx = @dims_in.collect{ |dim|
      if not dim.name.match("ndat") and @bc.shrink and not @ld then
        @dim_ngs
      elsif not dim.name.match("ndat") and @bc.shrink
        if @wavelet then
          if @wavelet == :decompose then
            [ BOAST::Dim(0, 1), BOAST::Dim( @filter.lowfil, dim + @filter.lowfil - 1 ) ]
          else
            [ BOAST::Dim( @filter.lowfil, dim + @filter.lowfil - 1 ), BOAST::Dim(0, 1) ]
          end
        else
          BOAST::Dim(@filter.lowfil, dim + @filter.lowfil - 1)
        end
      elsif not dim.name.match("ndat")
        if @wavelet then
          if @wavelet == :decompose then
            [ BOAST::Dim(0, 1), BOAST::Dim(0, dim - 1) ]
          else
            [ BOAST::Dim(0, dim - 1), BOAST::Dim(0, 1) ]
          end
        else
          BOAST::Dim(0, dim - 1)
        end
      else
        BOAST::Dim(0, dim - 1)
      end
    }
    dimx.flatten!
    dimy = @dims_out.collect{ |dim|
      if not dim.name.match("ndat") and @bc.grow and not @ld then
        @dim_nsg
      elsif not dim.name.match("ndat") and @bc.grow then
        if @wavelet then
          if @wavelet == :decompose then
            [ BOAST::Dim( -@filter.upfil, dim - @filter.upfil - 1 ), BOAST::Dim(0, 1) ]
          else
            [ BOAST::Dim(0, 1), BOAST::Dim( -@filter.upfil, dim - @filter.upfil - 1 ) ]
          end
        else
          BOAST::Dim( -@filter.upfil, dim - @filter.upfil - 1)
        end
      elsif not dim.name.match("ndat")
        if @wavelet then
          if @wavelet == :decompose then
            [ BOAST::Dim(0, dim - 1), BOAST::Dim(0, 1) ]
          else
            [ BOAST::Dim(0, 1), BOAST::Dim(0, dim - 1) ]
          end
        else
          BOAST::Dim(0, dim - 1)
        end
      else
        BOAST::Dim(0, dim - 1)
      end
    }
    if @transpose !=0  then
      dimy = dimy.rotate(@transpose)
    end
    dimy.flatten!
    return [dimx, dimy]
  end

  def compute_inner_loop_boundaries
    if @bc.grow then
      @line_start = -@filter.upfil
      @line_end = @dim_n - @filter.lowfil - 1
    else
      @line_start = 0
      @line_end = @dim_n - 1
    end
    @border_low = -@filter.lowfil
    @border_high = @dim_n - @filter.upfil
  end

  def params(dim, index=0)
    if @wavelet then
      dim[index] /= 2
    end
    vars=[]
    varsin=[]
    varsout=[]
    nd = { "n" => dim[index], "ndat" => 1, "ndat1" => 1, "ndat2" => 1 }
    if @dims.length == 2 then
      dim.each_index { |indx|
        nd["ndat"] *= dim[indx] if indx != index
      }
    else
      dim.each_index { |indx|
        nd["ndat1"] *= dim[indx] if indx < index
        nd["ndat2"] *= dim[indx] if indx > index
      }
    end
    n_push = lambda { |varsi, varso|
      if @bc.grow then
        varsi.push(nd["n"])
        if @wavelet then
          varso.push(nd["n"] + @filter.low.length - 1)
        else
          varso.push(nd["n"] + @filter.length - 1)
        end
      elsif @bc.shrink
        if @wavelet then
          varsi.push(nd["n"] + @filter.low.length - 1)
        else
          varsi.push(nd["n"] + @filter.length - 1)
        end
        varso.push(nd["n"])
      else
        varsi.push(nd["n"])
        varso.push(nd["n"])
      end
    }
    @dims.each { |dim|
      vars.push(nd[dim.name])
      if dim.name == "n" then
        n_push.call(varsin, varsout)
      else
        varsin.push(nd[dim.name])
        varsout.push(nd[dim.name])
      end
    }
    varsout.rotate!(@transpose)
    #input and output arrays
    if @ld then
      n_push.call(vars, vars)
    end
    if @wavelet then
      vars.push(NArray.float(*varsin,2).random)
      vars.push(NArray.float(*varsout,2))
    else
      vars.push(NArray.float(*varsin).random)
      vars.push(NArray.float(*varsin).random) if @kinetic and @kinetic != :inplace
      vars.push(NArray.float(*varsout))
      vars.push(NArray.float(*varsout)) if @kinetic and @transpose != 0
    end
    #accessory scalars
    nscal=0
    nscal+=1 if @a
    nscal+=1 if @a_x
    nscal+=1 if @a_y
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
      kernel = BOAST::CKernel::new
      print_header
      p = self.procedure(optim)
      BOAST::print p
      kernel.procedure = p
      #kernel.print #if @bc.free
      kernel.build(:openmp => true)
      t_mean = 0
      dimensions = opt_space.dimensions
      par = nil
      if dimensions.length < @dims.length then
        dimensions += [dimensions[0]]*(@dims.length-dimensions.length)
      end
      dimensions.length.times { |indx|
        stats_a = []
        par = self.params(dimensions.dup, indx)
        #puts par.inspect
        opt_space.repeat.times {
          stats_a.push kernel.run(*par)
        }
        stats_a.sort_by! { |a| a[:duration] }
        stats = stats_a.first
        #puts *par[0...@dims.length]
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
    if @wavelet then
      return 2 * n * ( 2 * @filter.low.length ) * ndat
    else
      return n * ( 2 * @filter.length ) * ndat
    end
  end

  def get_tt(tt_arr, unroll)
    #try to modify tt scalars into arrays of size unroll
    if tt_arr then
      if @wavelet then
        if @wavelet == :decompose then
          tt = [ BOAST::Real("lt", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true),
                 BOAST::Real("ht", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true) ]
        else
          tt = [ BOAST::Real("et", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true),
                 BOAST::Real("ot", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true) ]
        end
      else
        tt = BOAST::Real("tt", :dim => [ BOAST::Dim(0,unroll-1)], :allocate => true)
      end
    else
      if @wavelet then
        if @wavelet == :decompose then
          tt = [ (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("lt#{index}") },
                 (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("ht#{index}") } ]
        else
          tt = [ (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("et#{index}") },
                 (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("ot#{index}") } ]
        end
      else
        tt = (1..(unroll > 0 ? unroll : 1)).collect{ |index| BOAST::Real("tt#{index}") }
      end
    end
    return tt
  end

  def get_mods(mod_arr)
    if mod_arr then
      # the mod_arr behaves as a shrink operation
      #mods=BOAST::Int("mod_arr", :allocate => true, :dim => [@dim_ngs])
      mods=BOAST::Int("mod_arr", :allocate => true, 
                      :dim => [BOAST::Dim(@filter.lowfil - @filter.upfil,
                                          @filter.upfil - @filter.lowfil - 1)])
    else
      mods=nil
    end
    return mods
  end

  def get_constants
    return [@filter.lowfil, @filter.upfil]
  end

  def decl_filters
    if @wavelet then
      if @wavelet == :decompose then
        BOAST::decl @filter.low_even.fil
        BOAST::decl @filter.low_odd.fil
        BOAST::decl @filter.high_even.fil
        BOAST::decl @filter.high_odd.fil
      else
        BOAST::decl @filter.low_reverse_even.fil
        BOAST::decl @filter.low_reverse_odd.fil
        BOAST::decl @filter.high_reverse_even.fil
        BOAST::decl @filter.high_reverse_odd.fil
      end
    else
      BOAST::decl @filter.fil
    end
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
    tt = get_tt(tt_arr, unroll)
    iters =  (1..@dims.length).collect{ |index| BOAST::Int("i#{index}")}
    mods = get_mods(mod_arr)
    constants = get_constants

    return BOAST::Procedure(function_name, vars, constants ){
      decl_filters
      BOAST::decl *iters
      BOAST::decl l
      BOAST::decl *([tt].flatten)
      if mod_arr then
        BOAST::decl mods 
        #BOAST::print BOAST::For(l, @filter.lowfil, @dim_n -1 + @filter.upfil) {
        BOAST::print BOAST::For(l, mods.dimension.first.val1, mods.dimension.first.val2) {
          BOAST::print mods[l] === BOAST::modulo(l, @dim_n)
        }
      end
      BOAST::print @dotp === 0.0 if @options[:dotp]
      if BOAST::get_lang == BOAST::FORTRAN then
        BOAST::get_output.print("!$omp parallel default(shared)&\n")
        BOAST::get_output.print("!$omp reduction(+:#{dotp})&\n") if @options[:dotp]
        BOAST::get_output.print("!$omp private(#{iters.join(",")},#{l})&\n")
        BOAST::get_output.print("!$omp private(#{([tt].flatten).join(",")})\n")
      elsif BOAST::get_lang == BOAST::C then
        BOAST::get_output.print("#pragma omp parallel default(shared) #{@options[:dotp] ? "reduction(+:#{dotp})" : ""} private(#{iters.join(",")},#{l},#{([tt].flatten).join(",")})\n{\n")
      end

      convolution1d(iters, l, tt, mods, unrolled_dim, unroll, unroll_inner)

      BOAST::get_output.print("!$omp end parallel\n") if BOAST::get_lang == BOAST::FORTRAN
      BOAST::get_output.print("}\n")  if BOAST::get_lang == BOAST::C
    }
  end

  #here follows the internal operations for the convolution 1d
  def convolution1d(iters, l, t, mods, unro, unrolling_length, unroll_inner)
    convgen= lambda { |t,tlen,reliq|
      ises0 = startendpoints(@dims[@dim_indexes[0]], unro == @dim_indexes[0], unrolling_length, reliq)
      BOAST::get_output.print("!$omp do\n") if BOAST::get_lang == BOAST::FORTRAN
      BOAST::get_output.print("#pragma omp for\n") if BOAST::get_lang == BOAST::C
      BOAST::For(iters[@dim_indexes[0]], *ises0 ) {
        if @dim_indexes.length == 3 then
          ises1 = startendpoints(@dims[@dim_indexes[1]], unro == @dim_indexes[1], unrolling_length, reliq)
          BOAST::For(iters[@dim_indexes[1]], *ises1) {
            conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
          }.print
        else
          conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
        end
      }.print
      BOAST::get_output.print("!$omp end do\n") if BOAST::get_lang == BOAST::FORTRAN
    }
    #first without the reliq
    convgen.call(t,unrolling_length,false)
    #then with the reliq but only if the unrolling patterns need it
    convgen.call(t,1,true) if (unrolling_length > 1)
  end

  #returns the starting and ending points of the convolutions according to unrolling and unrolled dimension
  def startendpoints(dim,unroll,unrolling_length,in_reliq)
    istart= (in_reliq and unroll) ? (dim/unrolling_length)*unrolling_length : 0
    iend  = (unroll and not in_reliq) ? dim-unrolling_length : dim-1
    istep = (unroll and not in_reliq) ? unrolling_length : 1 
    return [istart,iend,istep]
  end

  def conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
    # the shrink operation contains the central part only
    iter = iters[@dim_indexes[-1]]
    if @bc.shrink then
      BOAST::For(iter, @line_start, @line_end) {
        for_conv(:center, iters, l, t, tlen, unro, mods, unroll_inner)
      }.print
    else
      BOAST::For(iter, @line_start, @border_low - 1) {
        for_conv(:begin, iters, l, t, tlen, unro, mods, unroll_inner)
      }.print
      BOAST::For(iter, @border_low, @border_high - 1) {
        for_conv(:center, iters, l, t, tlen, unro, mods, unroll_inner)
      }.print
      BOAST::For(iter, @border_high, @line_end) {
        for_conv(:end, iters, l, t, tlen, unro, mods, unroll_inner)
      }.print
    end
  end

  def get_loop_start_end( side, iters )
    processed_dim = @dim_indexes[-1]
    if ( @bc.free and side == :begin) then
      loop_start = BOAST::max(-iters[processed_dim], @filter.lowfil)
      loop_end   = @filter.upfil
    elsif ( @bc.free and side == :end) then
      loop_start = @filter.lowfil
      loop_end   = BOAST::min(@filter.upfil, @dims[processed_dim] - 1 - iters[processed_dim])
    else
      loop_start=@filter.lowfil
      loop_end=@filter.upfil
    end
    return [loop_start, loop_end]
  end

  def for_conv(side, iters, l, t, tlen, unro, mods, unroll_inner)
    processed_dim = @dim_indexes[-1]
    #t.each_index { |ind|
    (0...tlen).each{ |ind|
      #WARNING: the eks conditional here can be relaxed
      if @wavelet then
        BOAST::print t[0][ind] === 0
        BOAST::print t[1][ind] === 0
      else
        BOAST::print t[ind] === 0.0
      end
    }
    loop_start, loop_end = get_loop_start_end( side, iters)
    f = BOAST::For( l, loop_start, loop_end) {
      (0...tlen).each{ |ind|
       #t.each_index { |ind|
        if @bc.free or (side == :center) then
          i_in = input_index(unro, iters, ind,processed_dim,l)
        elsif mods then
          i_in = input_index(unro, iters, ind,processed_dim,l,nil,mods,side) 
        else
          i_in = input_index(unro, iters, ind,processed_dim,l,@dims[processed_dim])
        end
        if @wavelet then
          i_in[0].flatten!
          i_in[1].flatten!
          if @wavelet == :decompose then
            BOAST::print t[0][ind] === t[0][ind] + @x[*(i_in[0])]*@filter.low_even.fil[l]
            BOAST::print t[1][ind] === t[1][ind] + @x[*(i_in[0])]*@filter.high_even.fil[l]
            BOAST::print t[0][ind] === t[0][ind] + @x[*(i_in[1])]*@filter.low_odd.fil[l]
            BOAST::print t[1][ind] === t[1][ind] + @x[*(i_in[1])]*@filter.high_odd.fil[l]
          else
            BOAST::print t[0][ind] === t[0][ind] + @x[*(i_in[0])]*@filter.low_reverse_odd.fil[l]
            BOAST::print t[1][ind] === t[1][ind] + @x[*(i_in[0])]*@filter.low_reverse_even.fil[l]
            BOAST::print t[0][ind] === t[0][ind] + @x[*(i_in[1])]*@filter.high_reverse_odd.fil[l]
            BOAST::print t[1][ind] === t[1][ind] + @x[*(i_in[1])]*@filter.high_reverse_even.fil[l]
          end
        else
          BOAST::print t[ind] === t[ind] + @x[*i_in]*@filter.fil[l]
        end
      }
    }
    if unroll_inner then
      f.unroll
    else
      f.print
    end

    #t.each_index { |ind|
    (0...tlen).each{ |ind|
      i_out = output_index(unro, iters, ind)
      i_in = input_index(unro, iters, ind)
      if @wavelet then
        i_out[0].rotate!(@transpose)
        i_out[1].rotate!(@transpose)
        i_out[0].flatten!
        i_out[1].flatten!
        BOAST::print t[0][ind] === t[0][ind] * @a if @a
        BOAST::print t[1][ind] === t[1][ind] * @a if @a
        if @accumulate then #and not @init then
          BOAST::print t[0][ind] === @y[*i_out[0]] + t[0][ind]
          BOAST::print t[1][ind] === @y[*i_out[1]] + t[1][ind]
        elsif @a_y then #and not @init then
          BOAST::print t[0][ind] === @a_y * @y[*i_out[0]] + t[0][ind]
          BOAST::print t[1][ind] === @a_y * @y[*i_out[1]] + t[1][ind]
        end
        BOAST::print @y[*i_out[0]] === t[0][ind]
        BOAST::print @y[*i_out[1]] === t[1][ind]
      else
        i_out.rotate!(@transpose)
        BOAST::print t[ind] === t[ind] * @a if @a
        if @bc.grow and (@dotp or @a_x) and side != :center then
          if side == :begin then
            BOAST::print BOAST::If(iters[processed_dim] >= 0) {
              BOAST::print @dotp === @dotp + t[ind] * @x[*i_in] if @dotp
              BOAST::print t[ind] === t[ind] + @a_x * @x[*i_in] if @a_x
            }
          elsif side == :end then
            BOAST::print BOAST::If(iters[processed_dim] < @dims[processed_dim]) {
              BOAST::print @dotp === @dotp + t[ind] * @x[*i_in] if @dotp
              BOAST::print t[ind] === t[ind] + @a_x * @x[*i_in] if @a_x
            }
          end
        else
          BOAST::print @dotp === @dotp + t[ind] * @x[*i_in] if @dotp
          BOAST::print t[ind] === t[ind] + @a_x * @x[*i_in] if @a_x
        end

        #to be controlled in the case of non-orthorhombic cells for kinetic operations
        BOAST::print t[ind] === t[ind] + @x2[*i_in] if @x2
        if @accumulate or (@kinetic == :inplace and not @options[:zero_out])  then
          BOAST::print t[ind] === t[ind] + @y[*i_out]
        elsif @a_y then
          BOAST::print t[ind] === t[ind] +  @a_y * @y[*i_out]
        end

        BOAST::print @y[*i_out] === t[ind]
        BOAST::print @y2[*i_out] === @x[*i_in] if @kinetic and @transpose
      end
    }
  end


  #returns the indices of the output array according to the starting point in the input and of the
  ## processed dimension as well as the position in the convolution
  def input_index(unrolling_dim, iters,unroll_index,processed_dim=nil,lconv_index=nil,
                        ndim_processed=nil,wrapping_array=nil,side=nil)
    if ndim_processed then
      i_in = output_index_k_mod(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,ndim_processed)
    elsif wrapping_array then
      i_in = output_index_k_mod_arr(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    elsif processed_dim then
      i_in = output_index_k(unrolling_dim, iters,unroll_index,processed_dim,lconv_index)
    else
      i_in = output_index_unroll(unrolling_dim, iters,unroll_index)
    end
    if @wavelet then
      tmp = [[], []]
      (0...iters.length).each { |indx|
        if indx == @dim_indexes[-1] then
          if @wavelet == :decompose then
            tmp[0][indx] = [0, i_in[indx]]
            tmp[1][indx] = [1, i_in[indx]]
          else
            tmp[0][indx] = [i_in[indx], 0]
            tmp[1][indx] = [i_in[indx], 1]
          end
        else
          tmp[0][indx] = i_in[indx]
          tmp[1][indx] = i_in[indx]
        end
      }
      return tmp
    else
      return i_in
    end
  end


  #returns the indices of the output array according to the starting point in the input and of the
  ## processed dimension as well as the position in the convolution
  def output_index(unrolling_dim, iters,unroll_index,processed_dim=nil,lconv_index=nil,
                        ndim_processed=nil,wrapping_array=nil,side=nil)
    if ndim_processed then
      i_out = output_index_k_mod(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,ndim_processed)
    elsif wrapping_array then
      i_out = output_index_k_mod_arr(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    elsif processed_dim then
      i_out = output_index_k(unrolling_dim, iters,unroll_index,processed_dim,lconv_index)
    else
      i_out = output_index_unroll(unrolling_dim, iters,unroll_index)
    end
    if @wavelet then
      tmp = [[], []]
      (0...iters.length).each { |indx|
        if indx == @dim_indexes[-1] then
          if @wavelet == :decompose then
            tmp[0][indx] = [i_out[indx], 0]
            tmp[1][indx] = [i_out[indx], 1]
          else
            tmp[0][indx] = [0, i_out[indx]]
            tmp[1][indx] = [1, i_out[indx]]
          end
        else
          tmp[0][indx] = i_out[indx]
          tmp[1][indx] = i_out[indx]
        end
      }
      return tmp
    else
      return i_out
    end
  end

  #returns the indices of the output according to which of the directions is unrolled
  def output_index_unroll(unrolling_dim, iters,unroll_index)
    i_out=(0...iters.length).collect { |indx| unrolling_dim == indx ? iters[unrolling_dim] + (unroll_index) : iters[indx]}
    return i_out
  end
  # index of the convolution in the internal region, k=i+l, otherwise the index in the unrolling dimension
  def output_index_k(unrolling_dim, iters,unroll_index,processed_dim,lconv_index)
    i_out=output_index_unroll(unrolling_dim, iters,unroll_index)
    return (0...iters.length).collect { |indx| processed_dim == indx ? lconv_index +iters[processed_dim] : i_out[indx]}
  end
  # index in the external region wrapped around (periodic BC), thanks to the presence of the wrapping_array
  # if the side is :begin, the recipe is the usual, otherwise (:end) the recipe is subtracted
  # the the value of the processed dim. In this way the size of mod_arr is only dependent by the size of the 
  # filter which makes life easier for the compiler
  def output_index_k_mod_arr(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    i_out=output_index_unroll(unrolling_dim, iters,unroll_index)
    if (side == :end) then
      return (0...iters.length).collect { |indx| 
        processed_dim == indx ? wrapping_array[lconv_index +iters[processed_dim] - @dims[processed_dim]] : i_out[indx]}
    else
      return (0...iters.length).collect { |indx| processed_dim == indx ? wrapping_array[lconv_index +iters[processed_dim]] : i_out[indx]}
    end
  end
  # index in the external region wrapped around (periodic BC), where the wrapping is given by the integer division
  def output_index_k_mod(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,ndim_processed)
    i_out=output_index_unroll(unrolling_dim, iters,unroll_index)
    return (0...iters.length).collect { |indx| processed_dim == indx ?  lconv_index + iters[processed_dim] - ((iters[processed_dim]+lconv_index +  ndim_processed * 2 )/ndim_processed - 2) *ndim_processed  : i_out[indx]}
  end
end

class GenericConvolutionOperator1d
  attr_accessor :procs
  attr_reader :needed_subops
  def initialize(filter,options={})
    @filter = filter
    @ld = options[:ld]
    @m = options[:m]
    @wavelet = options[:wavelet]
    @kinetic = options[:kinetic]

    @vars = []
    @vars.push @ndim  = BOAST::Int( "d",    :dir => :in )
    @vars.push @idim  = BOAST::Int( "idim", :dir => :in )
    @vars.push @dims  = BOAST::Int( "n",    :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    @vars.push @bc    = BOAST::Int( "bc",   :dir => :in )
    if @ld then
      @vars.push @nx  = BOAST::Int( "nx", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      @vars.push @ny = BOAST::Int( "ny", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    end
    @vars.push @m     = BOAST::Int( "m", :dir => :in ) if @m
    @vars.push @x     = BOAST::Real("x",  :dir => :in,  :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @y     = BOAST::Real("y",  :dir => :out, :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @a = BOAST::Real("a",:dir => :in,:dim => [ BOAST::Dim(0, @ndim - 1)]) if options[:a]
    @vars.push @a_x = BOAST::Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = BOAST::Real("a_y",:dir => :in) if options[:a_y]

    @transpose = 0
    @options=options
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a_x)
    opts.delete(:a_y)
    dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
  end
end

class GenericConvolutionOperator
  attr_accessor :procs
  attr_reader :needed_subops
  def initialize(filter,options={})
    @filter = filter
    @ld = options[:ld]
    @m = options[:m]
    @wavelet = options[:wavelet]
    @kinetic = options[:kinetic]

    @vars = []
    @vars.push @ndim  = BOAST::Int( "d",  :dir => :in )
    @vars.push @dims  = BOAST::Int( "n",  :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    @vars.push @bc    = BOAST::Int( "bc", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    if @ld then
      @vars.push @nx  = BOAST::Int( "nx", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      @vars.push @ny = BOAST::Int( "ny", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    end
    @vars.push @m     = BOAST::Int( "m", :dir => :in ) if @m
    @vars.push @x     = BOAST::Real("x",  :dir => :in,  :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @y     = BOAST::Real("y",  :dir => :out, :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @w1 = BOAST::Real("w1", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
    @vars.push @w2 = BOAST::Real("w2", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
    @vars.push @a = BOAST::Real("a",:dir => :in,:dim => [ BOAST::Dim(0, @ndim - 1)]) if options[:a]
    @vars.push @a_x = BOAST::Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = BOAST::Real("a_y",:dir => :in) if options[:a_y]
    @vars.push @eks = BOAST::Real("eks",:dir => :out,:dim =>[ BOAST::Dim(0, @ndim - 1)]) if options[:eks]

    @transpose = 0
    @transpose = options[:transpose] if options[:transpose]
    @options=options
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a_x)
    opts.delete(:a_y)
    dim_indexes_a = []
    if @transpose == 0 then
      dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
    elsif @transpose == 1
      dim_indexes_a = [ [1, 0] ]
    elsif @transpose == -1
      dim_indexes_a = [ [0, 1] ]
    end
    opt_base = [{}]
    opt_base.push( { :a => @options[:a] } ) if @options[:a]
    opt_base.push( { :a_x => @options[:a_x] } ) if @options[:a_x]
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y]
    opts_bases = []
    1...opt_base.length { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        opt_bases.push(ch)
      }
    }
    BC::CONDITIONS.each { |bc|
      dim_indexes_a.each{ |dim_indexes|
        opts_bases.each { |opt|
          op = opt.dup
          op.update(opts)
          p = ConvolutionOperator1d::new(@filter, BC::new(bc), @transpose, dim_indexes, opts)
          @needed_subops[p.base_name] = p
          @procs[p.base_name] = p.procedure
        }
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
        ndat *= dims_actual[i] if i != indx
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
        ndat *= dims_actual[i] if i < indx
        ndat2 *= dims_actual[i] if i > indx
      }
      ni = n[indx]
      ndat_ni_ndat2 = [ndat, ni, ndat2]
      indexes = [2, 0, 1]
      return [ndat_ni_ndat2, indexes]
    }
    (0...n.length).each { |indx|
      if bc[indx] == BC::SHRINK then
        if @wavelet then
          dims_actual[indx] = n[indx] * 2 + @filter.length - 2
        else
          dims_actual[indx] = n[indx] + @filter.length - 1
        end
      else
        if @wavelet then
          dims_actual[indx] = n[indx] * 2
        else
          dims_actual[indx] = n[indx]
        end
      end
    }
    change_dims = lambda { |indx|
      if bc[indx] == BC::GROW then
        if @wavelet then
          dims_actual[indx] = n[indx] * 2 + @filter.length  - 2
        else
          dims_actual[indx] = n[indx] + @filter.length - 1
        end
      else
        if @wavelet then
          dims_actual[indx] = n[indx] * 2
        else
          dims_actual[indx] = n[indx]
        end
      end
    }
    if n.length == 1 then
      d = [ n[0], 1 ]
      d_indexes = [ 1, 0 ]
      dims.reverse! if @transpose == -1
      dim_indexes.reverse! if @transpose == -1
      cost = ConvolutionOperator1d::new(@filter, BC::new(bc[0]), @transpose, d_indexes, @options).cost( *d )
    else
      cost = 0
      dims, dim_indexes = compute_ni_ndat.call(0)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[0]), @transpose, dim_indexes, @options).cost( *dims )
      change_dims.call(0)
      dims_left = n.length - 1
      while dims_left > 1 do
        if @transpose == 0 then
          dims, dim_indexes = compute_ndat_ni_ndat2.call(n.length-dims_left)
        else
          dims, dim_indexes = compute_ni_ndat.call(n.length-dims_left)
        end
        cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), @transpose, dim_indexes, @options).cost( *dims )
        change_dims.call(n.length-dims_left)
        dims_left -= 1
      end
      dims, dim_indexes = compute_ni_ndat.call(n.length-1)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), @transpose, dim_indexes, @options).cost( *dims )
    end
    return cost * m
  end

  def procedure
    function_name = ""
    if @wavelet then
      if @wavelet == :decompose then
        function_name += "dwt_"
      else
        function_name += "idwt_"
      end
    end
    function_name += @filter.name
    function_name += "_ld" if @ld
    function_name += "_m" if @m
    p = BOAST::Procedure(function_name,@vars) {
      dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, @ndim - 1) ] ) if BOAST::get_lang == BOAST::FORTRAN
      dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, 16) ] ) if BOAST::get_lang == BOAST::C
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
          if @wavelet then
            BOAST::print dims_actual[i] === @nx[i] * 2
          else
            BOAST::print dims_actual[i] === @nx[i]
          end
        else
          BOAST::print BOAST::If(@bc[i] == BC::SHRINK, lambda {
            if @wavelet then
              BOAST::print dims_actual[i] === @dims[i] * 2 + @filter.length - 2
            else
              BOAST::print dims_actual[i] === @dims[i] + @filter.length - 1
            end
          }, lambda {
            if @wavelet then
              BOAST::print dims_actual[i] === @dims[i] * 2
            else
              BOAST::print dims_actual[i] === @dims[i]
            end
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

      opts = @options.dup
      opts.delete(:a_x)
      opts.delete(:a_y)

      print_call = lambda { |indx, init, last, datas, multi_conv|
        vars = dims
        if multi_conv then
          if dim_indexes.length == 2 then
            BOAST::print ndat_tot_in === dims[dim_indexes[0]]
          else
            BOAST::print ndat_tot_in === dims[dim_indexes[0]] * dims[dim_indexes[1]]
          end
          BOAST::print ndat_tot_out === ndat_tot_in
          BOAST::print ndat_tot_in === ndat_tot_in * dims_actual[indx]
          if @ld then
            if @wavelet then
              BOAST::print ndat_tot_out === ndat_tot_out * @ny[indx] * 2
            else
              BOAST::print ndat_tot_out === ndat_tot_out * @ny[indx]
            end
          end
          f = BOAST::For(j, 0, @m-1)
        end
        if @ld then
          vars.push @nx[indx]
          vars.push @ny[indx]
        end
        indx = @ndim - 1 - indx if @transpose == -1
        vars2 = []
        opt = {}
        vars2.push( @a[indx] ) if @options[:a]
        if init and @options[:a_x] then
          vars2.push( @a_x )
          opt[:a_x] = @options[:a_x]
        end
        if last and @options[:a_y] then
          vars2.push( @a_y )
          opt[:a_y] = @options[:a_y]
        end
        vars2.push( @eks[indx] ) if @options[:eks]
        dats = []
        opt.update( opts )
        BOAST::print BOAST::Case( @bc[indx], BC::PERIODIC, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::PERIODIC), @transpose, dim_indexes, opt).base_name

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
          BOAST::print dims_actual[indx] === @ny[indx]  if @ld

        }, BC::GROW, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::GROW),     @transpose, dim_indexes, opt).base_name

          if multi_conv then
            if not @ld then
              if @wavelet then
                BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 2 )
              else
                BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 1 )
              end
            end
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.print
          else
            dats = datas.dup
          end
          BOAST::print @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          if @ld then
            BOAST::print dims_actual[indx] === @ny[indx]  if @ld
          else
            if @wavelet then
              BOAST::print dims_actual[indx] === dims_actual[indx] + @filter.length - 2
            else
              BOAST::print dims_actual[indx] === dims_actual[indx] + @filter.length - 1
            end
          end
        }, BC::SHRINK, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::SHRINK),    @transpose, dim_indexes, opt).base_name

          if multi_conv then
            if not @ld then
              if @wavelet then
                BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 2 )
              else
                BOAST::print ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 1 )
              end
            end
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.print
          else
            dats = datas.dup
          end
          BOAST::print @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          if @ld then
            BOAST::print dims_actual[indx] === @ny[indx]  if @ld
          else
            if @wavelet then
              BOAST::print dims_actual[indx] === dims_actual[indx] - @filter.length + 2
            else
              BOAST::print dims_actual[indx] === dims_actual[indx] - @filter.length + 1
            end
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
        print_call.call( 0, true, true, datas, false )
      }, lambda {
        dims, dim_indexes = compute_ni_ndat.call( 0 )
        datas = [ @x, @w1 ]
        datas = [ @x, @y ] if not @options[:work]
        print_call.call( 0, true, false, datas, @m )
        BOAST::print dims_left === dims_left - 1
        BOAST::print i === 1
        BOAST::print BOAST::While( dims_left > 2 ) {
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w1, @w2 ]
          else
            datas = [ @w1, @w2 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @m )
          BOAST::print i === i + 1
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w2, @w1 ]
          else
            datas = [ @w2, @w1 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @m )
          BOAST::print i === i + 1
          BOAST::print dims_left === dims_left - 2
        }
        BOAST::print BOAST::If( dims_left == 2, lambda {
          if @transpose == 0 then
            dims, dim_indexes = compute_ndat_ni_ndat2.call( i )
          else
            dims, dim_indexes = compute_ni_ndat.call( i )
          end
          if @kinetic then
            datas = [ @x, @w1, @w2 ]
          else
            datas = [ @w1, @w2 ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, false, datas, @m )
          BOAST::print i === i + 1
          dims, dim_indexes = compute_ni_ndat.call( i )
          if @transpose == 0 then
            dims.reverse!
            dim_indexes.reverse!
          end
          if @kinetic then
            datas = [ @x, @w2, @y ]
          else
            datas = [ @w2, @y ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, true, datas, @m )
        }, lambda {
          dims, dim_indexes = compute_ni_ndat.call( i )
          if @transpose == 0 then
            dims.reverse! 
            dim_indexes.reverse!
          end
          if @kinetic then
            datas = [ @x, @w1, @y ]
          else
            datas = [ @w1, @y ]
          end
          datas = [ @x, @y ] if not @options[:work]
          print_call.call( i, false, true, datas, @m )
        })
      })
    }
    return [ p, @procs ]
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

def print_header
  if BOAST::get_lang == BOAST::C then
    BOAST::get_output.print "inline #{BOAST::Int::new.decl} modulo( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return (a+b)%b;}\n"
    BOAST::get_output.print "inline #{BOAST::Int::new.decl} min( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return a < b ? a : b;}\n"
    BOAST::get_output.print "inline #{BOAST::Int::new.decl} max( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return a > b ? a : b;}\n"
  end
end


