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
        if e.is_a?(String) then
          if e[0] == '-' then
            filt2.push(e[1..-1])
          else
            filt2.push("-"+e)
          end
        else
          filt2.push(-e)
        end
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

class PoissonFilter
  attr_reader :low
  attr_reader :name
  attr_reader :center
  attr_reader :fil_array, :filters_array
  attr_reader :length
  attr_reader :filters
  attr_reader :fil
  # extremes of the filter, calculated via its central point (integers)
  attr_reader :lowfil_val, :upfil_val
  # extremes of the filter, calculated via its central point (BOAST object)
  attr_reader :lowfil, :upfil
  # name of the filter
  attr_reader :name
  def initialize(name, filts,nord)
    @fil_array = filts.dup
    @center = nord/2
    #@center -= @center%2
    tmp_array = []
    @filters=[]
    index=0
    @fil_array.each_with_index { |e,i|
      @filters[i]=ConvolutionFilter::new(name+"_#{i}", e, @center)
        e.each{|val|
            tmp_array[index]=val
            index= index+1
        }    
    }
    arr = BOAST::ConstArray::new(tmp_array,BOAST::Real)
    @filters_array = BOAST::Real("#{name}_fil",:constant => arr,:dim => [ BOAST::Dim(0,(2*@center+1)*(2*@center+1)-1) ])
    @fil = @filters[@center].fil

    @length = @fil_array.length
    @name = name
    @lowfil_val = -center
    @upfil_val = @filters[@center].length - center - 1
    @lowfil = BOAST::Int("lowfil",:constant => @lowfil_val)
    @upfil = BOAST::Int("upfil",:constant => @upfil_val)
  end

def filters
  return @filters
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
#                          -2 : Non periodic BC, for nabla operators
class BoundaryConditions
  # conditions names
  PERIODIC = 0
  GROW = 1
  SHRINK = -1
  NPERIODIC = -2
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
  # determine if the convolution is non periodic
  attr_reader :nper
  # original id of the boundary conditions
  attr_reader :id

  def initialize(ibc)
    @id     = ibc
    @free   = (ibc != PERIODIC && ibc != NPERIODIC)
    @nper   = (ibc == NPERIODIC)
    @grow   = (ibc == GROW)
    @shrink = (ibc == SHRINK)
    @discard = (ibc == FREE)
    if not @free then
      if not @nper then
        @name = 'p'
      else
        @name = 'np'
      end
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
  def initialize(filter, bc, dim_indexes, options={})
    @filter = filter
    @bc = bc
    @transpose = options[:transpose]
    @dim_indexes = dim_indexes
    @ld = options[:ld]
    @kinetic = options[:kinetic]
    @wavelet = options[:wavelet]
    @poisson = options[:poisson]
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
    if @kinetic and @kinetic != :inplace and not options[:zero_out] then
      if @bc.grow then
        @vars.push @x2 = BOAST::Real("x2",:dir => :in, :dim => dimy, :restrict => true)
      else
        @vars.push @x2 = BOAST::Real("x2",:dir => :in, :dim => dimx, :restrict => true)
      end
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
    @vars.push @dot_in = BOAST::Real("dot_in",:dir => :out) if options[:dot_in]
    @cost = BOAST::Int("cost", :dir => :out)
    @options = options
    @base_name = ""
    @base_name += "s_" if BOAST::default_real_size == 4
    @base_name += "d_" if BOAST::default_real_size == 8
    if @wavelet then
      if @wavelet == :decompose then
        @base_name += "dwt_"
      else
        @base_name += "idwt_"
      end
    end
    @base_name += @filter.name + "_" + @bc.name + "_#{@dim_indexes.join('')}"
    @base_name += "_a" if @a
    @base_name += "_ain" if @a_x
    @base_name += "_ay" if @a_y
    @base_name += "_acc" if @accumulate
    @base_name += "_dotin" if @dot_in
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
    case BOAST::default_real_size
    when 4
      type = NArray::SFLOAT
    when 8
      type = NArray::FLOAT
    else
      raise "Unsupported precision!"
    end

    if @wavelet then
      vars.push(NArray::new(type, *varsin,2).random)
      vars.push(NArray::new(type, *varsout,2))
    else
      vars.push(NArray::new(type, *varsin).random)
      if @kinetic and @kinetic != :inplace then
        if @bc.grow then
          vars.push(NArray::new(type, *varsout).random)
        else
          vars.push(NArray::new(type, *varsin).random)
        end
      end
      vars.push(NArray::new(type, *varsout))
      vars.push(NArray::new(type, *varsout)) if @kinetic and @transpose != 0
    end
    #accessory scalars
    nscal=0
    nscal+=1 if @a
    nscal+=1 if @a_x
    nscal+=1 if @a_y
    nscal.times{vars.push(0.5)}
    vars.push(NArray::new(type, 1).random) if @dot_in
    return vars
  end

  def optimize(opt_space)
    opt_space=GenericOptimization::new if not opt_space
    t_best=Float::INFINITY
    p_best = nil
    already_tested = {}
    opt_space.each{ |optim|
      next if optim[:unrolled_dim_index] == 1 and @dims.length < 3
      #next if optim[:mod_arr] and @bc.free
      #puts optim
      kernel = BOAST::CKernel::new
      print_header
      p = self.procedure(optim)
      BOAST::pr p
      kernel.procedure = p
      next if already_tested[p.name]
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
          puts "#{indx} - [#{par[0...@dims.length].join(", ")}] - #{kernel.procedure.name}: #{stats[:duration]*1.0e3} ms #{self.cost(*par[0...@dims.length]) / (stats[:duration]*1.0e9)} GFlops"
          puts optim
        end
        t_mean += stats[:duration]
      }
      t_mean /= dimensions.length
      puts "#{kernel.procedure.name}: #{t_mean*1.0e3} ms #{self.cost(*par[0...@dims.length]) / (t_mean*1.0e9)} GFlops"
      already_tested[p.name] = true
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
          tt = [ (0..(unroll > 0 ? unroll - 1 : 0)).collect{ |index| BOAST::Real("lt#{index}") },
                 (0..(unroll > 0 ? unroll - 1 : 0)).collect{ |index| BOAST::Real("ht#{index}") } ]
        else
          tt = [ (0..(unroll > 0 ? unroll - 1 : 0)).collect{ |index| BOAST::Real("et#{index}") },
                 (0..(unroll > 0 ? unroll - 1 : 0)).collect{ |index| BOAST::Real("ot#{index}") } ]
        end
      else
        tt = (0..(unroll > 0 ? unroll - 1 : 0)).collect{ |index| BOAST::Real("tt#{index}") }
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
    elsif @poisson and @bc.nper  then
      BOAST::decl @filter.filters_array
      BOAST::decl @filter.fil
    else
      BOAST::decl @filter.fil
    end
  end

  def procedure(options={})
    #(unroll, unrolled_dim, use_mod, tt_arr)
    #default values
    @no_temp = false
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
    util = options[:util]

    function_name = @base_name + 
      "_u#{unroll}_#{unrolled_dim}_#{mod_arr}_#{tt_arr}_#{unroll_inner}"
    function_name += "_" + util.to_s if util

    if util == :cost then
      return BOAST::Procedure(function_name, @dims + [@cost] ){
        BOAST::decl ndat_t = BOAST::Int("ndat_t")
        BOAST::pr ndat_t === 1
        @dim_indexes[0...-1].each { |indx|
          BOAST::pr ndat_t === ndat_t * @dims[indx]
        }
        if @wavelet then
          BOAST::pr @cost === @dims[@dim_indexes.last] * 2 * 2 * @filter.low.length * ndat_t
        else
          BOAST::pr @cost === @dims[@dim_indexes.last] * 2 * @filter.length * ndat_t
        end
      }
    end

    l = BOAST::Int("l")
    tt = get_tt(tt_arr, unroll)
    iters =  (1..@dims.length).collect{ |index| BOAST::Int("i#{index}")}
    mods = get_mods(mod_arr)
    constants = get_constants

    return BOAST::Procedure(function_name, vars, constants ){
      decl_filters
      BOAST::decl *iters
      BOAST::decl l
      BOAST::decl *([tt].flatten) if not @no_temp
      if mod_arr then
        BOAST::decl mods 
        #BOAST::pr BOAST::For(l, @filter.lowfil, @dim_n -1 + @filter.upfil) {
        BOAST::pr BOAST::For(l, mods.dimension.first.val1, mods.dimension.first.val2) {
          BOAST::pr mods[l] === BOAST::modulo(l, @dim_n)
        }
      end
      BOAST::pr @dot_in === 0.0 if @options[:dot_in]
      BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: (@options[:dot_in] ? {"+" => dot_in} : nil ), private: iters + ( (not @no_temp) ? [tt] : [])) { 
        convolution1d(iters, l, tt, mods, unrolled_dim, unroll, unroll_inner)
      }
    }
  end

  #here follows the internal operations for the convolution 1d
  def convolution1d(iters, l, t, mods, unro, unrolling_length, unroll_inner)
    convgen= lambda { |t,tlen,reliq|
      ises0 = startendpoints(@dims[@dim_indexes[0]], unro == @dim_indexes[0], unrolling_length, reliq)
      BOAST::For(iters[@dim_indexes[0]], ises0[0], ises0[1], step: ises0[2], openmp: true ) {
        if @dim_indexes.length == 3 then
          ises1 = startendpoints(@dims[@dim_indexes[1]], unro == @dim_indexes[1], unrolling_length, reliq)
          BOAST::For(iters[@dim_indexes[1]], ises1[0], ises1[1], step: ises1[2]) {
            conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
          }.pr
        else
          conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
        end
      }.pr
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
      }.pr
    else
      BOAST::For(iter, @line_start, @border_low - 1) {
        for_conv(:begin, iters, l, t, tlen, unro, mods, unroll_inner)
      }.pr
      BOAST::For(iter, @border_low, @border_high - 1) {
        for_conv(:center, iters, l, t, tlen, unro, mods, unroll_inner)
      }.pr
      BOAST::For(iter, @border_high, @line_end) {
        for_conv(:end, iters, l, t, tlen, unro, mods, unroll_inner)
      }.pr
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

  def init_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    (0...tlen).each{ |ind|
      if @no_temp then
        i_out = output_index(unro, iters, ind)
        if @wavelet then
          i_out[0].rotate!(@transpose)
          i_out[1].rotate!(@transpose)
          i_out[0].flatten!
          i_out[1].flatten!
          if @a_y then
            BOAST::pr @y[*i_out[0]] === @a_y * @y[*i_out[0]]
            BOAST::pr @y[*i_out[1]] === @a_y * @y[*i_out[1]]
          elsif not @accumulate then
            BOAST::pr @y[*i_out[0]] === 0.0
            BOAST::pr @y[*i_out[1]] === 0.0
          end
          BOAST::pr @y[*i_out[0]] === @y[*i_out[0]] / @a if @a and ( @accumulate or @a_y )
          BOAST::pr @y[*i_out[1]] === @y[*i_out[1]] / @a if @a and ( @accumulate or @a_y )
        else
          i_out.rotate!(@transpose)
          if @a_y then
            BOAST::pr @y[*i_out] === @a_y * @y[*i_out]
          elsif not ( @accumulate or (@kinetic == :inplace and not @options[:zero_out]) ) then
            BOAST::pr @y[*i_out] === 0.0
          end
          BOAST::pr @y[*i_out] === @y[*i_out] / @a if @a and ( @accumulate or (@kinetic == :inplace and not @options[:zero_out]) )
        end
      else
        #WARNING: the eks conditional here can be relaxed
        if @wavelet then
          BOAST::pr t[0][ind] === 0.0
          BOAST::pr t[1][ind] === 0.0
        else
          BOAST::pr t[ind] === 0.0
        end
      end
    }
  end

  def compute_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    processed_dim = @dim_indexes[-1]
    (0...tlen).each{ |ind|
      if @no_temp then
        i_out = output_index(unro, iters, ind)
        if @wavelet then
          i_out[0].rotate!(@transpose)
          i_out[1].rotate!(@transpose)
          i_out[0].flatten!
          i_out[1].flatten!
          out_even = @y[*i_out[0]]
          out_odd = @y[*i_out[1]]
        else
          i_out.rotate!(@transpose)
          out = @y[*i_out]
        end
      else
        if @wavelet then
          out_even = t[0][ind]
          out_odd = t[1][ind]
        else
          out = t[ind]
        end
      end
      if @bc.free or (side == :center) then
        i_in = input_index(unro, iters, ind, processed_dim, l, nil, nil, side)
      elsif mods then
        i_in = input_index(unro, iters, ind, processed_dim, l, nil, mods, side) 
      else
        i_in = input_index(unro, iters, ind, processed_dim, l, @dims[processed_dim],nil, side)
      end
      if @wavelet then
        i_in[0].flatten!
        i_in[1].flatten!
        if @wavelet == :decompose then
          BOAST::pr out_even === out_even + @x[*(i_in[0])]*@filter.low_even.fil[l]
          BOAST::pr out_odd === out_odd + @x[*(i_in[0])]*@filter.high_even.fil[l]
          BOAST::pr out_even === out_even + @x[*(i_in[1])]*@filter.low_odd.fil[l]
          BOAST::pr out_odd === out_odd + @x[*(i_in[1])]*@filter.high_odd.fil[l]
        else
          BOAST::pr out_even === out_even + @x[*(i_in[0])]*@filter.low_reverse_odd.fil[l]
          BOAST::pr out_odd === out_odd + @x[*(i_in[0])]*@filter.low_reverse_even.fil[l]
          BOAST::pr out_even === out_even + @x[*(i_in[1])]*@filter.high_reverse_odd.fil[l]
          BOAST::pr out_odd === out_odd + @x[*(i_in[1])]*@filter.high_reverse_even.fil[l]
        end
      elsif @poisson and @bc.nper and (side != :center) then
        if(side == :begin) then
#          BOAST::pr testindex = iters[processed_dim] - @filter.upfil
          BOAST::pr out === out + @x[*i_in]*@filter.filters_array[((iters[processed_dim] - (@filter.upfil))*(2*@filter.center+1)) + (@filter.center)*(2*@filter.center+1) + (l) + (@filter.center) ]
        else
#          BOAST::pr testindex = iters[processed_dim] + @filter.lowfil - @dims[processed_dim]+1
          BOAST::pr out === out + @x[*i_in]*@filter.filters_array[((iters[processed_dim] + (@filter.center) - @dims[processed_dim] +1)*(2*@filter.center+1)) + (@filter.center)*(2*@filter.center+1) + (l) + @filter.center]
        end

      else
        BOAST::pr out === out + @x[*i_in]*@filter.fil[l]
      end
    }
  end

  
  def post_process_and_store_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    processed_dim = @dim_indexes[-1]
    (0...tlen).each{ |ind|
      i_out = output_index(unro, iters, ind)
      i_in = input_index(unro, iters, ind)
      if @wavelet then
        i_out[0].rotate!(@transpose)
        i_out[1].rotate!(@transpose)
        i_out[0].flatten!
        i_out[1].flatten!
      else
        i_out.rotate!(@transpose)
      end
      if @no_temp then
        if @wavelet then
          out_even = @y[*i_out[0]]
          out_odd = @y[*i_out[1]]
        else
          out = @y[*i_out]
        end
      else
        if @wavelet then
          out_even = t[0][ind]
          out_odd = t[1][ind]
        else
          out = t[ind]
        end
      end
      if @wavelet then
        BOAST::pr out_even === out_even * @a if @a
        BOAST::pr out_odd === out_odd * @a if @a
        if not @no_temp then
          if @accumulate then #and not @init then
            BOAST::pr out_even === @y[*i_out[0]] + out_even
            BOAST::pr out_odd === @y[*i_out[1]] + out_odd
          elsif @a_y then #and not @init then
            BOAST::pr out_even === @a_y * @y[*i_out[0]] + out_even
            BOAST::pr out_odd === @a_y * @y[*i_out[1]] + out_odd
          end
          BOAST::pr @y[*i_out[0]] === out_even
          BOAST::pr @y[*i_out[1]] === out_odd
        end
      else
        BOAST::pr out === out * @a if @a
        if @bc.grow and (@dot_in or @a_x) and side != :center then
          if side == :begin then
            BOAST::pr BOAST::If(iters[processed_dim] >= 0) {
              BOAST::pr @dot_in === @dot_in + out * @x[*i_in] if @dot_in
              BOAST::pr out === out + @a_x * @x[*i_in] if @a_x
            }
          elsif side == :end then
            BOAST::pr BOAST::If(iters[processed_dim] < @dims[processed_dim]) {
              BOAST::pr @dot_in === @dot_in + out * @x[*i_in] if @dot_in
              BOAST::pr out === out + @a_x * @x[*i_in] if @a_x
            }
          end
        else
          BOAST::pr @dot_in === @dot_in + out * @x[*i_in] if @dot_in
          BOAST::pr out === out + @a_x * @x[*i_in] if @a_x
        end

        #to be controlled in the case of non-orthorhombic cells for kinetic operations
        BOAST::pr out === out + @x2[*i_in] if @x2
        if not @no_temp then
          if @accumulate or (@kinetic == :inplace and not @options[:zero_out])  then
            BOAST::pr out === out + @y[*i_out]
          elsif @a_y then
            BOAST::pr out === out +  @a_y * @y[*i_out]
          end
          BOAST::pr @y[*i_out] === out
        end
        BOAST::pr @y2[*i_out] === @x[*i_in] if @kinetic and @transpose != 0
      end
    }
  end

  def for_conv(side, iters, l, t, tlen, unro, mods, unroll_inner)

    init_values(side, iters, l, t, tlen, unro, mods, unroll_inner)

    loop_start, loop_end = get_loop_start_end( side, iters)

    f = BOAST::For( l, loop_start, loop_end) {
      compute_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    }
    if unroll_inner then
      f.unroll
    else
      f.pr
    end

    post_process_and_store_values(side, iters, l, t, tlen, unro, mods, unroll_inner)

  end


  #returns the indices of the output array according to the starting point in the input and of the
  ## processed dimension as well as the position in the convolution
  def input_index(unrolling_dim, iters,unroll_index,processed_dim=nil,lconv_index=nil,
                        ndim_processed=nil,wrapping_array=nil,side=nil)
    if @poisson and @bc.nper and (side != :center) then
      i_in = output_index_k_nper(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    elsif ndim_processed then
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
    if @poisson and @bc.nper and (side != :center) then
      i_out = output_index_k_nper(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    elsif ndim_processed then
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

  # index in the external region wrapped around (non periodic BC), 
  # all the indexes close to the bound are the same, (filter will differ)
  def output_index_k_nper(unrolling_dim, iters,unroll_index,processed_dim,lconv_index,wrapping_array,side)
    i_out=output_index_unroll(unrolling_dim, iters,unroll_index)
    if (side == :end) then
      return (0...iters.length).collect { |indx| 
        processed_dim == indx ? (lconv_index) + (@filter.lowfil) + @dims[processed_dim] -1 : i_out[indx]}
    else
      return (0...iters.length).collect { |indx| processed_dim == indx ? lconv_index + (@filter.upfil) : i_out[indx]}
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
    @narr = options[:narr]
    @wavelet = options[:wavelet]
    @kinetic = options[:kinetic]
    @poisson = options[:poisson]

    @vars = []
    @vars.push @ndim  = BOAST::Int( "d",    :dir => :in )
    @vars.push @idim  = BOAST::Int( "idim", :dir => :in )
    @vars.push @dims  = BOAST::Int( "n",    :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    @vars.push @bc    = BOAST::Int( "bc",   :dir => :in )
    if @ld then
      @vars.push @nx  = BOAST::Int( "nx", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
      @vars.push @ny = BOAST::Int( "ny", :dir => :in, :dim => [ BOAST::Dim(0, @ndim - 1) ] )
    end
    @vars.push @narr     = BOAST::Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = BOAST::Real("x",  :dir => :in,  :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @y     = BOAST::Real("y",  :dir => :out, :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @a = BOAST::Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = BOAST::Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = BOAST::Real("a_y",:dir => :in) if options[:a_y]
    @vars.push @dot_in = BOAST::Real("dot_in",:dir => :out) if options[:dot_in]
    @cost = BOAST::Int( "cost", :dir => :out )

    @transpose = 0
    @options=options.dup
    @options[:transpose] = 0
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a)
    opts.delete(:a_x)
    opts.delete(:a_y)
    dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
    opt_base = []
    opt_base.push( { :a   => @options[:a]   } ) if @options[:a]
    opt_base.push( { :a_x => @options[:a_x] } ) if @options[:a_x]
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y]
    opts_bases = []
    (0..opt_base.length).each { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        opts_bases.push(ch)
      }
    }
    if @poisson then
    conditions = BC::CONDITIONS << BC::NPERIODIC 
    else
    conditions = BC::CONDITIONS
    end
    conditions.each { |bc|
      dim_indexes_a.each { |dim_indexes|
        opts_bases.each { |opt|
          op = opt.dup
          op.update(opts)
          p = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, op)
          @needed_subops[p.base_name] = p
          @procs[p.base_name] = p.procedure
          @procs[p.base_name+"_cost"] = p.procedure( :util => :cost )
        }
      }
    }
    p = self.procedure(:cost).first
    @procs[p.name] = p
  end

  def optimize(opt_space=nil)
    @needed_subops.each { |name, subop|
      @procs[name] = subop.optimize(opt_space)
    }
  end

  def cost( idim, n, bc, m = 1 )
    ndim = n.length
    if ndim == 1 then
      dims = [ n[0], 1 ]
      dim_indexes = [1, 0]
    else
      ni = n[idim]
      if idim == 0 then
        ndat_left = nil
      else
        ndat_left = 1
        n[0...idim].each { |val| ndat_left *= val }
      end
      if idim == ndim - 1 then
        ndat_right = nil
      else
        ndat_right = 1
        n[(idim+1)..-1].each { |val| ndat_right *= val }
      end
      if idim == 0 then
        dim_indexes = [1, 0]
      elsif idim == ndim - 1 then
        dim_indexes = [0, 1]
      else
        dim_indexes = [2, 0, 1]
      end
      dims = []
      dims.push( ndat_left ) if ndat_left
      dims.push( ni )
      dims.push( ndat_right ) if ndat_right
    end
    cost = ConvolutionOperator1d::new(@filter, BC::new(bc[idim]), dim_indexes, @options).cost( *dims )
    return cost * m
  end

  def procedure(util = nil)
    function_name = ""
    function_name += "d_" if BOAST::default_real_size == 8
    function_name += "s_" if BOAST::default_real_size == 4
    if @wavelet then
      if @wavelet == :decompose then
        function_name += "s0s1_1d_"
      else
        function_name += "s1s0_1d_"
      end
    else
      function_name += "s0s0_1d_"
    end
    function_name += @filter.name
    function_name += "dotin" if @dot_in
    function_name += "_" + util.to_s if util

    vv = @vars
    vv += [ @cost ] if util == :cost
    
    p = BOAST::Procedure( function_name, vv ) {
      ni = BOAST::Int "ni"
      ndat_left = BOAST::Int "ndat_left"
      ndat_right = BOAST::Int "ndat_right"
      nti = BOAST::Int "nti" if @narr
      nto = BOAST::Int "nto" if @narr
      i = BOAST::Int "i"
      j = BOAST::Int "j"
      tmp_cost = BOAST::Int "c"
      BOAST::decl i, j, ni, ndat_left, ndat_right, tmp_cost
      BOAST::decl nti, nto if @narr
      if @narr and @ld then
        BOAST::pr nti === @nx[@idim]
        BOAST::pr nto === @ny[@idim]
      elsif @narr then
        BOAST::pr BOAST::If(@bc[i] == BC::SHRINK, lambda {
          if @wavelet then
            BOAST::pr nti === @dims[@idim] + @filter.low.length - 1
          else
            BOAST::pr nti === @dims[@idim] + @filter.length - 1
          end
          BOAST::pr nto === @dims[@idim]
        }, @bc[i] == BC::GROW, lambda {
          BOAST::pr nti === @dims[@idim]
          if @wavelet then
            BOAST::pr nto === @dims[@idim] + @filter.low.length - 1
          else
            BOAST::pr nto === @dims[@idim] + @filter.length - 1
          end
        })
      end
      if @narr and @wavelet then
        BOAST::pr nti === nti * 2
        BOAST::pr nto === nto * 2
      end
      dims = []
      dim_indexes = []
      dats = []
      if @narr then
        f = BOAST::For(j, 0, @narr-1)
        dats[0] = (@x[nti*j+1]).address
        dats[1] = (@y[nto*j+1]).address
      else
        dats[0] = @x
        dats[1] = @y
      end

      selected_bc = nil

      print_call_generic = lambda { |bc, a, a_x, a_y|
        opts = @options.dup
        opts.delete(:a) if not a
        opts.delete(:a_x) if not a_x
        opts.delete(:a_y) if not a_y
        vars = []
        vars.push( @a ) if a
        vars.push( @a_x ) if a_x
        vars.push( @a_y ) if a_y
        vars.push( @dot_in ) if @dot_in
        vars.push( tmp_cost.address ) if util == :cost
        lds = []
        lds.push( @nx[@idim] ) if @ld
        lds.push( @ny[@idim] ) if @ld
        procname = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, opts).base_name
        procname += "_" + util.to_s if util
        args = dims + lds + dats + vars
        if util == :cost then
          BOAST::pr @cost === 0
          args = dims + [tmp_cost]
        end
        f.open if @narr
          BOAST::pr @procs[procname].call( *args )
          BOAST::pr @cost === @cost + tmp_cost if util == :cost
        f.close if @narr
      }

      print_call_param_a_y = lambda { |bc, a, a_x|
        if @a_y then
          BOAST::pr BOAST::If(@a_y == 0.0, lambda {
            print_call_generic.call(bc, a, a_x, false)
          }, lambda {
            print_call_generic.call(bc, a, a_x, true)
          })
        else
          print_call_generic.call(bc, a, a_x, false)
        end
      }

      print_call_param_a_x = lambda { |bc, a|
        if @a_x then
          BOAST::pr BOAST::If(@a_x == 0.0, lambda {
            print_call_param_a_y.call(bc, a, false)
          }, lambda {
            print_call_param_a_y.call(bc, a, true )
          })
        else
          print_call_param_a_y.call(bc, a, false)
        end
      }

      print_call_param_a = lambda { |bc|
        if @a then
          BOAST::pr BOAST::If(@a == 1.0, lambda {
            print_call_param_a_x.call(bc, false)
          }, lambda {
            print_call_param_a_x.call(bc, true )
          })
        else
          print_call_param_a_x.call(bc, false)
        end
      }

      print_call = lambda {
        BOAST::pr BOAST::Case( @bc, BC::PERIODIC, lambda {
          print_call_param_a.call( BC::PERIODIC )
        }, BC::GROW, lambda {
          print_call_param_a.call( BC::GROW )
        }, BC::SHRINK, lambda {
          print_call_param_a.call( BC::SHRINK )
        }, BC::NPERIODIC, lambda {
          print_call_param_a.call( BC::NPERIODIC ) if @poisson 
        })
      }
      BOAST::pr BOAST::If( @idim == 0, lambda {
        BOAST::pr ndat_right === 1
        BOAST::pr BOAST::For( i, 1, @ndim - 1 ) {
          if @ld and util != :cost then
            BOAST::pr ndat_right === ndat_right * @nx[i]
          else
            BOAST::pr ndat_right === ndat_right * @dims[i]
          end
          BOAST::pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          BOAST::pr nti === nti * ndat_right
          BOAST::pr nto === nto * ndat_right
        end
        dim_indexes = [1,0]
        dims = [@dims[@idim], ndat_right]
        print_call.call
      }, @idim == @ndim - 1, lambda {
        BOAST::pr ndat_left === 1
        BOAST::pr BOAST::For( i, 0, @ndim - 2 ) {
          if @ld and util != :cost then
            BOAST::pr ndat_left === ndat_left * @nx[i]
          else
            BOAST::pr ndat_left === ndat_left * @dims[i]
          end
          BOAST::pr ndat_left === ndat_left * 2 if @wavelet
        }
        if @narr then
          BOAST::pr nti === nti * ndat_left
          BOAST::pr nto === nto * ndat_left
        end
        dim_indexes = [0,1]
        dims = [ndat_left, @dims[@idim]]
        print_call.call
      }, lambda {
        BOAST::pr ndat_left === 1
        BOAST::pr ndat_right === 1
        BOAST::pr BOAST::For( i, 0, @idim - 1 ) {
          if @ld and util != :cost then
            BOAST::pr ndat_left === ndat_left * @nx[i]
          else
            BOAST::pr ndat_left === ndat_left * @dims[i]
          end
          BOAST::pr ndat_left === ndat_left * 2 if @wavelet
        }
        BOAST::pr BOAST::For( i, @idim + 1, @ndim - 1 ) {
          if @ld and util != :cost then
            BOAST::pr ndat_right === ndat_right * @nx[i]
          else
            BOAST::pr ndat_right === ndat_right * @dims[i]
          end
          BOAST::pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          BOAST::pr nti === nti * ndat_left * ndat_right
          BOAST::pr nto === nto * ndat_left * ndat_right
        end
        dim_indexes = [2,0,1]
        dims = [ndat_left, @dims[@idim], ndat_right]
        print_call.call
      })
    }
    return [ p, @procs ]
  end

end

class GenericConvolutionOperator
  attr_accessor :procs
  attr_reader :needed_subops
  def initialize(filter,options={})
    @filter = filter
    @ld = options[:ld]
    @narr = options[:narr]
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
    @vars.push @narr     = BOAST::Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = BOAST::Real("x",  :dir => :in,  :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @y     = BOAST::Real("y",  :dir => :out, :restrict => true, :dim => [ BOAST::Dim() ] )
    @vars.push @w1 = BOAST::Real("w1", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
    @vars.push @w2 = BOAST::Real("w2", :dir => :inout, :restrict => true, :dim => [ BOAST::Dim() ] ) if options[:work]
    @vars.push @a = BOAST::Real("a",:dir => :in,:dim => [ BOAST::Dim(0, @ndim - 1)]) if options[:a]
    @vars.push @a_x = BOAST::Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = BOAST::Real("a_y",:dir => :in) if options[:a_y]
    @vars.push @dot_in = BOAST::Real("dot_in",:dir => :out,:dim =>[ BOAST::Dim(0, @ndim - 1)]) if options[:dot_in]

    @transpose = 0
    @transpose = options[:transpose] if options[:transpose]
    @options=options
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a)
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
    opt_base = []
    opt_base.push( { :a => @options[:a] } ) if @options[:a]
    opt_base.push( { :a_x => @options[:a_x] } ) if @options[:a_x]
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y]
    opts_bases = []
    (0..opt_base.length).each { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        opts_bases.push(ch)
      }
    }
    BC::CONDITIONS.each { |bc|
      dim_indexes_a.each{ |dim_indexes|
        opts_bases.each { |opt|
          op = opt.dup
          op.update(opts)
          p = ConvolutionOperator1d::new(@filter, BC::new(bc), dim_indexes, op)
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
      cost = ConvolutionOperator1d::new(@filter, BC::new(bc[0]), d_indexes, @options).cost( *d )
    else
      cost = 0
      dims, dim_indexes = compute_ni_ndat.call(0)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[0]), dim_indexes, @options).cost( *dims )
      change_dims.call(0)
      dims_left = n.length - 1
      while dims_left > 1 do
        if @transpose == 0 then
          dims, dim_indexes = compute_ndat_ni_ndat2.call(n.length-dims_left)
        else
          dims, dim_indexes = compute_ni_ndat.call(n.length-dims_left)
        end
        cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), dim_indexes, @options).cost( *dims )
        change_dims.call(n.length-dims_left)
        dims_left -= 1
      end
      dims, dim_indexes = compute_ni_ndat.call(n.length-1)
      cost += ConvolutionOperator1d::new(@filter, BC::new(bc[n.length-dims_left]), dim_indexes, @options).cost( *dims )
    end
    return cost * m
  end

  def procedure
    function_name = ""
    function_name += "d_" if BOAST::default_real_size == 8
    function_name += "s_" if BOAST::default_real_size == 4
    if @wavelet then
      if @wavelet == :decompose then
        function_name += "dwt_"
      else
        function_name += "idwt_"
      end
    end
    function_name += @filter.name
    function_name += "_ld" if @ld
    function_name += "_narr" if @narr
    p = BOAST::Procedure(function_name,@vars) {
      dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, @ndim - 1) ] ) if BOAST::get_lang == BOAST::FORTRAN
      dims_actual = BOAST::Int( "dims_actual", :allocate => true, :dim => [ BOAST::Dim(0, 16) ] ) if BOAST::get_lang == BOAST::C
      dims_left   = BOAST::Int "dims_left"
      ni = BOAST::Int "ni"
      ndat = BOAST::Int "ndat"
      ndat2 = BOAST::Int "ndat2"
      ndat_tot_in = BOAST::Int "nti" if @narr
      ndat_tot_out = BOAST::Int "nto" if @narr
      i = BOAST::Int "i"
      j = BOAST::Int "j"
      BOAST::decl i, j, dims_actual, dims_left, ni, ndat, ndat2
      BOAST::decl ndat_tot_in, ndat_tot_out if @narr
      dims = []
      dim_indexes = []
      BOAST::pr dims_left === @ndim
      BOAST::pr BOAST::For( i, 0, @ndim - 1 ) {
        if @ld then
          if @wavelet then
            BOAST::pr dims_actual[i] === @nx[i] * 2
          else
            BOAST::pr dims_actual[i] === @nx[i]
          end
        else
          BOAST::pr BOAST::If(@bc[i] == BC::SHRINK, lambda {
            if @wavelet then
              BOAST::pr dims_actual[i] === @dims[i] * 2 + @filter.length - 2
            else
              BOAST::pr dims_actual[i] === @dims[i] + @filter.length - 1
            end
          }, lambda {
            if @wavelet then
              BOAST::pr dims_actual[i] === @dims[i] * 2
            else
              BOAST::pr dims_actual[i] === @dims[i]
            end
          })
        end
      }
      compute_ni_ndat = lambda { |indx|
        indx = @ndim - 1 - indx if @transpose == -1
        BOAST::pr ni === @dims[indx]
        BOAST::pr ndat === 1
        BOAST::pr BOAST::For(j, 0, @ndim - 1) {
          BOAST::pr BOAST::If( j != indx ) {
            BOAST::pr ndat === ndat * dims_actual[j]
          }
        }
        d = [ ni, ndat ]
        d_indexes = [ 1, 0]
        d.reverse! if @transpose == -1
        d_indexes.reverse! if @transpose == -1
        return [ d, d_indexes ]
      }
      compute_ndat_ni_ndat2 = lambda { |indx|
        BOAST::pr ni === @dims[indx]
        BOAST::pr ndat === 1
        BOAST::pr ndat2 === 1
        BOAST::pr BOAST::For(j, 0, @ndim - 1) {
          BOAST::pr BOAST::If( j < indx, lambda {
            BOAST::pr ndat === ndat * dims_actual[j]
          }, j > indx, lambda {
            BOAST::pr ndat2 === ndat2 * dims_actual[j]
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
            BOAST::pr ndat_tot_in === dims[dim_indexes[0]]
          else
            BOAST::pr ndat_tot_in === dims[dim_indexes[0]] * dims[dim_indexes[1]]
          end
          BOAST::pr ndat_tot_out === ndat_tot_in
          BOAST::pr ndat_tot_in === ndat_tot_in * dims_actual[indx]
          if @ld then
            if @wavelet then
              BOAST::pr ndat_tot_out === ndat_tot_out * @ny[indx] * 2
            else
              BOAST::pr ndat_tot_out === ndat_tot_out * @ny[indx]
            end
          end
          f = BOAST::For(j, 0, @narr-1)
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
        vars2.push( @dot_in[indx] ) if @options[:dot_in]
        dats = []
        opt.update( opts )
        BOAST::pr BOAST::Case( @bc[indx], BC::PERIODIC, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::PERIODIC), dim_indexes, opt).base_name

          if multi_conv then
            BOAST::pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.pr
          else
            dats = datas.dup
          end
          BOAST::pr @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          BOAST::pr dims_actual[indx] === @ny[indx]  if @ld

        },BC::NPERIODIC, lambda {
          if@poisson then
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::NPERIODIC), dim_indexes, opt).base_name

          if multi_conv then
            BOAST::pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.pr
          else
            dats = datas.dup
          end
          BOAST::pr @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          BOAST::pr dims_actual[indx] === @ny[indx]  if @ld
          end
        }, BC::GROW, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::GROW), dim_indexes, opt).base_name

          if multi_conv then
            if not @ld then
              if @wavelet then
                BOAST::pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 2 )
              else
                BOAST::pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 1 )
              end
            end
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.pr
          else
            dats = datas.dup
          end
          BOAST::pr @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          if @ld then
            BOAST::pr dims_actual[indx] === @ny[indx]  if @ld
          else
            if @wavelet then
              BOAST::pr dims_actual[indx] === dims_actual[indx] + @filter.length - 2
            else
              BOAST::pr dims_actual[indx] === dims_actual[indx] + @filter.length - 1
            end
          end
        }, BC::SHRINK, lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::SHRINK), dim_indexes, opt).base_name

          if multi_conv then
            if not @ld then
              if @wavelet then
                BOAST::pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 2 )
              else
                BOAST::pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 1 )
              end
            end
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            f.pr
          else
            dats = datas.dup
          end
          BOAST::pr @procs[procname].call( *vars, *dats, *vars2 )
          f.close if multi_conv
          if @ld then
            BOAST::pr dims_actual[indx] === @ny[indx]  if @ld
          else
            if @wavelet then
              BOAST::pr dims_actual[indx] === dims_actual[indx] - @filter.length + 2
            else
              BOAST::pr dims_actual[indx] === dims_actual[indx] - @filter.length + 1
            end
          end
        })
      }

      BOAST::pr BOAST::If( @ndim == 1 , lambda {
        conv_number = 1
        conv_number = @narr if @narr
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
        print_call.call( 0, true, false, datas, @narr )
        BOAST::pr dims_left === dims_left - 1
        BOAST::pr i === 1
        BOAST::pr BOAST::While( dims_left > 2 ) {
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
          print_call.call( i, false, false, datas, @narr )
          BOAST::pr i === i + 1
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
          print_call.call( i, false, false, datas, @narr )
          BOAST::pr i === i + 1
          BOAST::pr dims_left === dims_left - 2
        }
        BOAST::pr BOAST::If( dims_left == 2, lambda {
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
          print_call.call( i, false, false, datas, @narr )
          BOAST::pr i === i + 1
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
          print_call.call( i, false, true, datas, @narr )
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
          print_call.call( i, false, true, datas, @narr )
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

def print_header(macro = false)
  if BOAST::get_lang == BOAST::C then
    if macro then
      BOAST::get_output.print "#define modulo(a, b) ((a+b)%(b))\n"
      BOAST::get_output.print "#define min( a, b) ((a) < (b) ? a : b)\n"
      BOAST::get_output.print "#define max( a, b) ((a) > (b) ? a : b)\n"
    else
      BOAST::get_output.print "static inline #{BOAST::Int::new.decl} modulo( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return (a+b)%b;}\n"
      BOAST::get_output.print "static inline #{BOAST::Int::new.decl} min( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return a < b ? a : b;}\n"
      BOAST::get_output.print "static inline #{BOAST::Int::new.decl} max( #{BOAST::Int::new.decl} a, #{BOAST::Int::new.decl} b) { return a > b ? a : b;}\n"
    end
  end
end


