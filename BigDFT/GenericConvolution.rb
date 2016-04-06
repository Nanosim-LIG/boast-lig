require 'BOAST'
require 'narray'
include BOAST

register_funccall("modulo")
register_funccall("min")
register_funccall("max")

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
    arr = ConstArray::new(@fil_array,Real)
    @fil = Real("#{name}_fil",:constant => arr,:dim => [ Dim((-center),(@fil_array.length - center -1)) ])
    @lowfil_val = -center
    @upfil_val = @fil_array.length - center - 1
    @lowfil = Int("lowfil",:constant => @lowfil_val)
    @upfil = Int("upfil",:constant => @upfil_val)
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
    arr = ConstArray::new(tmp_array,Real)
    @filters_array = Real("#{name}_fil",:constant => arr,:dim => [ Dim(0,(2*@center+1)*(2*@center+1)-1) ])
    @fil = @filters[@center].fil

    @length = @fil_array.length
    @name = name
    @lowfil_val = -center
    @upfil_val = @filters[@center].length - center - 1
    @lowfil = Int("lowfil",:constant => @lowfil_val)
    @upfil = Int("upfil",:constant => @upfil_val)
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
  attr_reader :openmp
  
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
    vector_length = 1
    vector_length = options[:vector_length] if options[:vector_length]
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
    @openmp = true
    @openmp = false if options[:openmp] == false

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
    space[:vector_length_range] = [vector_length].flatten
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

    @vars.push @x = Real("x",:dir => :in, :dim => dimx, :restrict => true)
    if @kinetic and @kinetic != :inplace and not options[:zero_out_work] then
      if @bc.grow then
        @vars.push @x2 = Real("x2",:dir => :in, :dim => dimy, :restrict => true)
      else
        @vars.push @x2 = Real("x2",:dir => :in, :dim => dimx, :restrict => true)
      end
    end
    @vars.push @y = Real("y",:dir => options[:a_y] ? :inout : :out, :dim => dimy, :restrict => true)
    if @kinetic and @transpose != 0 then
      @vars.push @y2 =  Real("y2", :dir => :out, :dim => dimy, :restrict => true)
    end
    @vars.push @a = Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x] #and init
    if options[:a_y] then
      if options[:a_y] == 1.0 then
        @accumulate = true
      else
        @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y]
      end
    end
    @vars.push @dot_in = Real("dot_in",:dir => :out) if options[:dot_in]
    @dot_in_tmp = nil
    @cost = Int("cost", :dir => :out)
    @options = options
    @base_name = ""
    @base_name += "s_" if default_real_size == 4
    @base_name += "d_" if default_real_size == 8
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
    @base_name += "_x2" if @x2
  end

  def compute_dims
    @dim_n = Int("n",:dir =>:in)
    if @ld then
      @nx = Int("nx",:dir =>:in)
      @ny = Int("ny",:dir =>:in)
    else
      @nx = @dim_n
      @ny = @dim_n
    end
    @dims = [@dim_n]
    @dims_in = [@nx]
    @dims_out = [@ny]
    if (dim_indexes.length == 3) then
      ndat1 = Int("ndat1",:dir =>:in)
      ndat2 = Int("ndat2",:dir =>:in)
      @dims     = [ndat1] + @dims     + [ndat2]
      @dims_in  = [ndat1] + @dims_in  + [ndat2]
      @dims_out = [ndat1] + @dims_out + [ndat2]
    elsif dim_indexes.last == 0
      ndat = Int("ndat",:dir =>:in)
      @dims     = @dims     + [ndat]
      @dims_in  = @dims_in  + [ndat]
      @dims_out = @dims_out + [ndat]
    else
      ndat = Int("ndat",:dir =>:in)
      @dims     = [ndat] + @dims
      @dims_in  = [ndat] + @dims_in
      @dims_out = [ndat] + @dims_out
    end

    if @wavelet then
      if @wavelet == :decompose then
        @dim_ngs = [ Dim(0, 1), Dim( @filter.lowfil, @dim_n + @filter.upfil - 1 ) ]
        @dim_nsg = [ Dim( -@filter.upfil, @dim_n - @filter.lowfil - 1 ), Dim(0, 1) ] 
      else
        @dim_ngs = [ Dim( @filter.lowfil, @dim_n + @filter.upfil - 1 ), Dim(0, 1) ]
        @dim_nsg = [ Dim(0, 1), Dim( -@filter.upfil, @dim_n - @filter.lowfil - 1 ) ]
      end
    else
      #growed dimension, to be used either for extremes or for mod_arr
      @dim_ngs = Dim( @filter.lowfil, @dim_n + @filter.upfil  - 1)
      #dimensions corresponding to the output of a grow operation
      @dim_nsg = Dim(-@filter.upfil,  @dim_n - @filter.lowfil - 1)
    end

  end

  def compute_dimx_dimy
    dimx = @dims_in.collect{ |dim|
      if not dim.name.match("ndat") and @bc.shrink and not @ld then
        @dim_ngs
      elsif not dim.name.match("ndat") and @bc.shrink
        if @wavelet then
          if @wavelet == :decompose then
            [ Dim(0, 1), Dim( @filter.lowfil, dim + @filter.lowfil - 1 ) ]
          else
            [ Dim( @filter.lowfil, dim + @filter.lowfil - 1 ), Dim(0, 1) ]
          end
        else
          Dim(@filter.lowfil, dim + @filter.lowfil - 1)
        end
      elsif not dim.name.match("ndat")
        if @wavelet then
          if @wavelet == :decompose then
            [ Dim(0, 1), Dim(0, dim - 1) ]
          else
            [ Dim(0, dim - 1), Dim(0, 1) ]
          end
        else
          Dim(0, dim - 1)
        end
      else
        Dim(0, dim - 1)
      end
    }
    dimx.flatten!
    dimy = @dims_out.collect{ |dim|
      if not dim.name.match("ndat") and @bc.grow and not @ld then
        @dim_nsg
      elsif not dim.name.match("ndat") and @bc.grow then
        if @wavelet then
          if @wavelet == :decompose then
            [ Dim( -@filter.upfil, dim - @filter.upfil - 1 ), Dim(0, 1) ]
          else
            [ Dim(0, 1), Dim( -@filter.upfil, dim - @filter.upfil - 1 ) ]
          end
        else
          Dim( -@filter.upfil, dim - @filter.upfil - 1)
        end
      elsif not dim.name.match("ndat")
        if @wavelet then
          if @wavelet == :decompose then
            [ Dim(0, dim - 1), Dim(0, 1) ]
          else
            [ Dim(0, 1), Dim(0, dim - 1) ]
          end
        else
          Dim(0, dim - 1)
        end
      else
        Dim(0, dim - 1)
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
    case default_real_size
    when 4
      type = NArray::SFLOAT
    when 8
      type = NArray::FLOAT
    else
      raise "Unsupported precision!"
    end

    align = 64
    if @wavelet then
      vars.push(ANArray::new(type, align, *varsin,2).random!)
      vars.push(ANArray::new(type, align, *varsout,2))
    else
      vars.push(ANArray::new(type, align, *varsin).random!)
      if @x2 then
        if @bc.grow then
          vars.push(ANArray::new(type, align, *varsout).random!)
        else
          vars.push(ANArray::new(type, align, *varsin).random!)
        end
      end
      vars.push(ANArray::new(type, align, *varsout))
      vars.push(ANArray::new(type, align, *varsout)) if @kinetic and @transpose != 0
    end
    #accessory scalars
    nscal=0
    nscal+=1 if @a
    nscal+=1 if @a_x
    nscal+=1 if @a_y
    nscal.times{vars.push(0.5)}
    vars.push(NArray::new(type, 1).random!) if @dot_in
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
      kernel = CKernel::new
      print_header
      p = self.procedure(optim)
      pr p
      kernel.procedure = p
      next if already_tested[p.name]
      #kernel.print #if @bc.free
      kernel.build(:openmp => opt_space.openmp)
      dimensions = opt_space.dimensions
      par = nil
      if dimensions.length < @dims.length then
        dimensions += [dimensions[0]]*(@dims.length-dimensions.length)
      end
      stats_a = []
      par = self.params(dimensions.dup,@dim_indexes.last)
      #puts par.inspect
      opt_space.repeat.times {
        stats_a.push kernel.run(*par)
      }
      stats_a.sort_by! { |a| a[:duration] }
      stats = stats_a.first
      #puts *par[0...@dims.length]
      if get_verbose then
        puts "#{optim} - [#{par[0...@dims.length].join(", ")}] - #{kernel.procedure.name}: #{stats[:duration]*1.0e3} ms #{self.cost(*par[0...@dims.length]) / (stats[:duration]*1.0e9)} GFlops"
        puts optim
      end
      t_min = stats[:duration]
      puts "#{kernel.procedure.name}: #{t_min*1.0e3} ms #{self.cost(*par[0...@dims.length]) / (t_min*1.0e9)} GFlops"
      already_tested[p.name] = true
      if t_best > t_min then
        t_best = t_min
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

  def get_tt(tt_arr, unroll, vec_len)
    #try to modify tt scalars into arrays of size unroll
    if tt_arr then
      if @wavelet then
        if @wavelet == :decompose then
          tt = [ Real("lt", :dim => [ Dim(0,unroll/vec_len-1)], :allocate => true, :vector_length => vec_len),
                 Real("ht", :dim => [ Dim(0,unroll/vec_len-1)], :allocate => true, :vector_length => vec_len) ]
        else
          tt = [ Real("et", :dim => [ Dim(0,unroll/vec_len-1)], :allocate => true, :vector_length => vec_len),
                 Real("ot", :dim => [ Dim(0,unroll/vec_len-1)], :allocate => true, :vector_length => vec_len) ]
        end
      else
        tt = Real("tt", :dim => [ Dim(0,unroll/vec_len-1)], :allocate => true, :vector_length => vec_len)
      end
    else
      if @wavelet then
        if @wavelet == :decompose then
          tt = [ (0..(unroll > 0 ? unroll/vec_len - 1 : 0)).collect { |index| Real("lt#{index}", :vector_length => vec_len) },
                 (0..(unroll > 0 ? unroll/vec_len - 1 : 0)).collect { |index| Real("ht#{index}", :vector_length => vec_len) } ]
        else
          tt = [ (0..(unroll > 0 ? unroll/vec_len - 1 : 0)).collect { |index| Real("et#{index}", :vector_length => vec_len) },
                 (0..(unroll > 0 ? unroll/vec_len - 1 : 0)).collect { |index| Real("ot#{index}", :vector_length => vec_len) } ]
        end
      else
        tt = (0..(unroll > 0 ? unroll/vec_len - 1 : 0)).collect { |index| Real("tt#{index}", :vector_length => vec_len) }
      end
    end
    return tt
  end

  def get_mods(mod_arr)
    if mod_arr then
      # the mod_arr behaves as a shrink operation
      #mods=Int("mod_arr", :allocate => true, :dim => [@dim_ngs])
      mods=Int("mod_arr", :allocate => true, 
                      :dim => [Dim(@filter.lowfil - @filter.upfil,
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
        decl @filter.low_even.fil
        decl @filter.low_odd.fil
        decl @filter.high_even.fil
        decl @filter.high_odd.fil
      else
        decl @filter.low_reverse_even.fil
        decl @filter.low_reverse_odd.fil
        decl @filter.high_reverse_even.fil
        decl @filter.high_reverse_odd.fil
      end
    elsif @poisson and @bc.nper  then
      decl @filter.filters_array
      decl @filter.fil
    else
      decl @filter.fil
    end
  end

  def procedure(options={})
    #(unroll, unrolled_dim, use_mod, tt_arr)
    #default values
    register_funccall("modulo")
    unroll = 1
    vec_len = 1
    mod_arr = true
    tt_arr = false
    unroll_inner = true
    unroll = options[:unroll] if options[:unroll]
    mod_arr = options[:mod_arr] if not options[:mod_arr].nil?
    tt_arr = options[:tt_arr] if not options[:tt_arr].nil?
    unroll_inner = options[:unroll_inner] if not options[:unroll_inner].nil?

    unrolled_dim=@dim_indexes[0]
    unrolled_dim=@dim_indexes[options[:unrolled_dim_index]] if @dim_indexes.length > 2 and options[:unrolled_dim_index]

    vec_len = options[:vector_length] if get_lang == C and not @dot_in and unrolled_dim == 0 and options[:vector_length] and options[:vector_length] <= @filter.length and unroll % options[:vector_length] == 0

    mod_arr = false if @bc.free
    util = options[:util]

    function_name = @base_name + 
      "_u#{unroll}_v#{vec_len}_#{unrolled_dim}_#{mod_arr}_#{tt_arr}_#{unroll_inner}"
    function_name += "_" + util.to_s if util

    if util == :cost then
      return Procedure(function_name, @dims + [@cost] ){
        decl ndat_t = Int("ndat_t")
        pr ndat_t === 1
        @dim_indexes[0...-1].each { |indx|
          pr ndat_t === ndat_t * @dims[indx]
        }
        if @wavelet then
          pr @cost === @dims[@dim_indexes.last] * 2 * 2 * @filter.low.length * ndat_t
        else
          pr @cost === @dims[@dim_indexes.last] * 2 * @filter.length * ndat_t
        end
      }
    end

    l = Int("l")
    tt = get_tt(tt_arr, unroll, vec_len)
    iters =  (1..@dims.length).collect{ |index| Int("i#{index}")}
    mods = get_mods(mod_arr)
    constants = get_constants

    return Procedure(function_name, vars, constants ){
      decl_filters
      decl *iters
      decl l
      decl *([tt].flatten)
      decl @filter_val = tt[0].copy("filter_val") if not @wavelet
      if mod_arr then
        decl mods 
        #pr For(l, @filter.lowfil, @dim_n -1 + @filter.upfil) {
        pr For(l, mods.dimension.first.val1, mods.dimension.first.val2) {
          pr mods[l] === modulo(l, @dim_n)
        }
      end
      vec_len = [tt].flatten[0].type.vector_length
      if @options[:dot_in] and vec_len > 1 then
        decl @dot_in_tmp = @dot_in.copy("dot_in_tmp", :vector_length => vec_len, :dir => nil, :direction => nil)
      else
        @dot_in_tmp = @dot_in
      end
      pr @dot_in.set(0.0) if @options[:dot_in]
      pr OpenMP::Parallel(default: :shared, reduction: (@options[:dot_in] ? {"+" => dot_in} : nil ), private: iters + [l] + [tt] + ( @filter_val ? [@filter_val] : [] )) { 
        convolution1d(iters, l, tt, mods, unrolled_dim, unroll, unroll_inner)
      }
    }
  end

  #here follows the internal operations for the convolution 1d
  def convolution1d(iters, l, t, mods, unro, unrolling_length, unroll_inner)
    vec_len = [t].flatten[0].type.vector_length
    convgen= lambda { |t,tlen,reliq|
      ises0 = startendpoints(@dims[@dim_indexes[0]], unro == @dim_indexes[0], unrolling_length, reliq, vec_len)
      pr For(iters[@dim_indexes[0]], ises0[0], ises0[1], step: ises0[2], openmp: true ) {
        if @dim_indexes.length == 3 then
          ises1 = startendpoints(@dims[@dim_indexes[1]], unro == @dim_indexes[1], unrolling_length, reliq, vec_len)
          pr For(iters[@dim_indexes[1]], ises1[0], ises1[1], step: ises1[2]) {
            conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
          }
        else
          conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
        end
        if @options[:dot_in] and vec_len > 1 then
          Reduce(@dot_in_tmp, @dot_in)
        end
      }
    }
    #first without the reliq
    convgen.call(t,unrolling_length,false)
    #then with the reliq but only if the unrolling patterns need it
    convgen.call(t,vec_len,true) if (unrolling_length > vec_len)
  end

  #returns the starting and ending points of the convolutions according to unrolling and unrolled dimension
  def startendpoints(dim,unroll,unrolling_length,in_reliq,vec_len)
    istart= (in_reliq and unroll) ? (dim/unrolling_length)*unrolling_length : 0
    iend  = (unroll and not in_reliq) ? dim-unrolling_length : dim-1
    istep = (unroll and not in_reliq) ? unrolling_length : ( unroll ? vec_len : 1 )
    return [istart,iend,istep]
  end

  def conv_lines(iters, l, t, tlen, unro, mods, unroll_inner)
    # the shrink operation contains the central part only
    iter = iters[@dim_indexes[-1]]
    if @bc.shrink then
      pr For(iter, @line_start, @line_end) {
        for_conv(:center, iters, l, t, tlen, unro, mods, unroll_inner)
      }
    else
      pr For(iter, @line_start, @border_low - 1) {
        for_conv(:begin, iters, l, t, tlen, unro, mods, unroll_inner)
      }
      pr For(iter, @border_low, @border_high - 1) {
        for_conv(:center, iters, l, t, tlen, unro, mods, unroll_inner)
      }
      pr For(iter, @border_high, @line_end) {
        for_conv(:end, iters, l, t, tlen, unro, mods, unroll_inner)
      }
    end
  end

  def get_loop_start_end( side, iters )
    register_funccall("min")
    register_funccall("max")
    processed_dim = @dim_indexes[-1]
    if ( @bc.free and side == :begin) then
      loop_start = max(-iters[processed_dim], @filter.lowfil)
      loop_end   = @filter.upfil
    elsif ( @bc.free and side == :end) then
      loop_start = @filter.lowfil
      loop_end   = min(@filter.upfil, @dims[processed_dim] - 1 - iters[processed_dim])
    else
      loop_start=@filter.lowfil
      loop_end=@filter.upfil
    end
    return [loop_start, loop_end]
  end

  def init_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    vec_len = [t].flatten[0].type.vector_length
    (0...tlen).step(vec_len).each{ |ind|
      #WARNING: the eks conditional here can be relaxed
      tt_ind = ind/vec_len
      if @wavelet then
        pr t[0][tt_ind].set( 0.0 )
        pr t[1][tt_ind].set( 0.0 )
      else
        pr t[tt_ind].set( 0.0 )
      end
    }
  end

  def compute_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    processed_dim = @dim_indexes[-1]
    vec_len = [t].flatten[0].type.vector_length
    (0...tlen).step(vec_len).each{ |ind|
      tt_ind = ind/vec_len
      if @wavelet then
        out_even = t[0][tt_ind]
        out_odd  = t[1][tt_ind]
      else
        out = t[tt_ind]
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
          pr out_even === FMA(Load(@x[*(i_in[0])], out_even), Set(@filter.low_even.fil[l] , out_even), out_even)
          pr out_odd  === FMA(Load(@x[*(i_in[0])], out_odd ), Set(@filter.high_even.fil[l], out_odd ), out_odd )
          pr out_even === FMA(Load(@x[*(i_in[1])], out_even), Set(@filter.low_odd.fil[l]  , out_even), out_even)
          pr out_odd  === FMA(Load(@x[*(i_in[1])], out_odd ), Set(@filter.high_odd.fil[l] , out_odd ), out_odd )
        else
          pr out_even === FMA(Load(@x[*(i_in[0])], out_even), Set(@filter.low_reverse_odd.fil[l]  , out_even), out_even)
          pr out_odd  === FMA(Load(@x[*(i_in[0])], out_odd ), Set(@filter.low_reverse_even.fil[l] , out_odd ), out_odd )
          pr out_even === FMA(Load(@x[*(i_in[1])], out_even), Set(@filter.high_reverse_odd.fil[l] , out_even), out_even)
          pr out_odd  === FMA(Load(@x[*(i_in[1])], out_odd ), Set(@filter.high_reverse_even.fil[l], out_odd ), out_odd )
        end
      else
        pr out === FMA(Load(@x[*i_in].set_align(out.type.total_size), out), @filter_val, out) # + @x[*i_in]*@filter.fil[l]
      end
    }
  end

  
  def post_process_and_store_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    processed_dim = @dim_indexes[-1]
    vec_len = [t].flatten[0].type.vector_length
    (0...tlen).step(vec_len).each{ |ind|
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
      tt_ind = ind/vec_len
      if @wavelet then
        out_even = t[0][tt_ind]
        out_odd  = t[1][tt_ind]
      else
        out = t[tt_ind]
      end
      if @wavelet then
        pr out_even === out_even * Set(@a, out_even) if @a
        pr out_odd  === out_odd  * Set(@a, out_odd ) if @a
        if @accumulate then #and not @init then
          pr out_even === Load(@y[*i_out[0]], out_even) + out_even
          pr out_odd  === Load(@y[*i_out[1]], out_odd)  + out_odd
        elsif @a_y then #and not @init then
          pr out_even === FMA(Load(@y[*i_out[0]], out_even), Set(@a_y, out_even), out_even)
          pr out_odd  === FMA(Load(@y[*i_out[1]], out_odd ), Set(@a_y, out_odd ), out_odd )
        end
        pr @y[*i_out[0]] === out_even
        pr @y[*i_out[1]] === out_odd
      else
        pr out === out * Set(@a, out) if @a
        finish_block = lambda {
          pr @dot_in_tmp === FMA(Load(@x[*i_in].set_align(out.type.total_size), out), out, @dot_in_tmp) if @dot_in_tmp #reduction !!!!!!!!!!!!!!!!!!!!!!!
          pr out === FMA(Load(@x[*i_in], out), Set(@a_x, out), out) if @a_x
        }
        if @bc.grow and (@dot_in or @a_x) and side != :center then
          if side == :begin then
            pr If(iters[processed_dim] >= 0, &finish_block)
          elsif side == :end then
            pr If(iters[processed_dim] < @dims[processed_dim], &finish_block)
          end
        else
          finish_block.call
        end

        #to be controlled in the case of non-orthorhombic cells for kinetic operations
        pr out === out + Load(@x2[*i_in], out) if @x2
        if @accumulate or (@kinetic == :inplace and not @options[:zero_out])  then
          pr out === out + Load(@y[*i_out].set_align(out.type.total_size), out)
        elsif @a_y then
          pr out === FMA(Set(@a_y, out), Load(@y[*i_out].set_align(out.type.total_size), out), out)
        end
        pr @y[*i_out].set_align(out.type.total_size) === out
        pr @y2[*i_out] === Load(@x[*i_in], out) if @kinetic and @transpose != 0
      end
    }
  end

  def for_conv(side, iters, l, t, tlen, unro, mods, unroll_inner)

    init_values(side, iters, l, t, tlen, unro, mods, unroll_inner)

    loop_start, loop_end = get_loop_start_end( side, iters)

    f = For( l, loop_start, loop_end) {
      if not (@poisson and @bc.nper and (side != :center))
        pr @filter_val === Set(@filter.fil[l], t[0]) unless @wavelet 
      else 
        processed_dim = @dim_indexes[-1]
        if side == :begin then
          pr @filter_val === Set(@filter.filters_array[((iters[processed_dim] - (@filter.upfil))*(2*@filter.center+1)) + (@filter.center)*(2*@filter.center+1) + (l) + (@filter.center) ], t[0])
        elsif side == :end then
          pr @filter_val === Set(@filter.filters_array[((iters[processed_dim] + (@filter.center) - @dims[processed_dim] +1)*(2*@filter.center+1)) + (@filter.center)*(2*@filter.center+1) + (l) + @filter.center], t[0])
        end
      end
      compute_values(side, iters, l, t, tlen, unro, mods, unroll_inner)
    }
    if unroll_inner then
      pr f.unroll
    else
      pr f
    end

    post_process_and_store_values(side, iters, l, t, tlen, unro, mods, unroll_inner)

  end


  #returns the indices of the input array according to the starting point in the input and of the
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
    @vars.push @ndim  = Int( "d",    :dir => :in )
    @vars.push @idim  = Int( "idim", :dir => :in )
    @vars.push @dims  = Int( "n",    :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    @vars.push @bc    = Int( "bc",   :dir => :in )
    if @ld then
      @vars.push @nx  = Int( "nx", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
      @vars.push @ny = Int( "ny", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    end
    @vars.push @narr     = Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = Real("x",  :dir => :in,  :restrict => true, :dim => [ Dim() ] )
    @vars.push @y     = Real("y",  :dir => options[:a_y] ? :inout : :out, :restrict => true, :dim => [ Dim() ] )
    @vars.push @a = Real("a",:dir => :in) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y] and options[:a_y] != 1.0
    @vars.push @dot_in = Real("dot_in",:dir => :out) if options[:dot_in]
    @cost = Int( "cost", :dir => :out )

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
    opt_base.push( { :a_y => @options[:a_y] } ) if @options[:a_y] and @options[:a_y] != 1.0
    opts_bases = []
    (0..opt_base.length).each { |indx|
      opt_base.combination(indx) { |c|
        ch = {}
        c.each { |item|
          ch.update(item)
        }
        ch.update( { :a_y => @options[:a_y] } ) if @options[:a_y] and @options[:a_y] == 1.0
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
    function_name += "d_" if default_real_size == 8
    function_name += "s_" if default_real_size == 4
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
    
    p = Procedure( function_name, vv ) {
      ndat_left = Int "ndat_left"
      ndat_right = Int "ndat_right"
      nti = Int "nti" if @narr
      nto = Int "nto" if @narr
      i = Int "i"
      j = Int "j"
      tmp_cost = Int "c"
      decl i, ndat_left, ndat_right
      decl tmp_cost if util == :cost
      decl nti, nto, j if @narr
      if @narr and @ld then
        pr nti === @nx[@idim]
        pr nto === @ny[@idim]
      elsif @narr then
        pr If(@bc[i] == BC::SHRINK => lambda {
          if @wavelet then
            pr nti === @dims[@idim] + @filter.low.length - 1
          else
            pr nti === @dims[@idim] + @filter.length - 1
          end
          pr nto === @dims[@idim]
        }, @bc[i] == BC::GROW => lambda {
          pr nti === @dims[@idim]
          if @wavelet then
            pr nto === @dims[@idim] + @filter.low.length - 1
          else
            pr nto === @dims[@idim] + @filter.length - 1
          end
        })
      end
      if @narr and @wavelet then
        pr nti === nti * 2
        pr nto === nto * 2
      end
      dims = []
      dim_indexes = []
      dats = []
      if @narr then
        f = For(j, 0, @narr-1)
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
        opts.update( { :a_y => 1.0 } ) if @options[:a_y] and @options[:a_y] == 1.0
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
          pr @cost === 0
          args = dims + [tmp_cost.address]
        end
        opn f if @narr
          pr @procs[procname].call( *args )
          pr @cost === @cost + tmp_cost if util == :cost
        close f if @narr
      }

      print_call_param_a_y = lambda { |bc, a, a_x|
        if @a_y then
          pr If(@a_y == 0.0 => lambda {
            print_call_generic.call(bc, a, a_x, false)
          }, else: lambda {
            print_call_generic.call(bc, a, a_x, true)
          })
        else
          print_call_generic.call(bc, a, a_x, false)
        end
      }

      print_call_param_a_x = lambda { |bc, a|
        if @a_x then
          pr If(@a_x == 0.0 => lambda {
            print_call_param_a_y.call(bc, a, false)
          }, else: lambda {
            print_call_param_a_y.call(bc, a, true )
          })
        else
          print_call_param_a_y.call(bc, a, false)
        end
      }

      print_call_param_a = lambda { |bc|
        if @a then
          pr If(@a == 1.0 => lambda {
            print_call_param_a_x.call(bc, false)
          }, else: lambda {
            print_call_param_a_x.call(bc, true )
          })
        else
          print_call_param_a_x.call(bc, false)
        end
      }

      print_call = lambda {
        case_args = {
          BC::PERIODIC => lambda {
            print_call_param_a.call( BC::PERIODIC )
          },
          BC::GROW => lambda {
            print_call_param_a.call( BC::GROW )
          },
          BC::SHRINK => lambda {
            print_call_param_a.call( BC::SHRINK )
          }
        }
        case_args[BC::NPERIODIC] = lambda { print_call_param_a.call( BC::NPERIODIC ) } if @poisson
        pr Case( @bc, case_args)
      }
      pr If( @idim == 0 => lambda {
        pr ndat_right === 1
        pr For( i, 1, @ndim - 1 ) {
          if @ld and util != :cost then
            pr ndat_right === ndat_right * @nx[i]
          else
            pr ndat_right === ndat_right * @dims[i]
          end
          pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_right
          pr nto === nto * ndat_right
        end
        dim_indexes = [1,0]
        dims = [@dims[@idim], ndat_right]
        print_call.call
      }, @idim == @ndim - 1 => lambda {
        pr ndat_left === 1
        pr For( i, 0, @ndim - 2 ) {
          if @ld and util != :cost then
            pr ndat_left === ndat_left * @nx[i]
          else
            pr ndat_left === ndat_left * @dims[i]
          end
          pr ndat_left === ndat_left * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_left
          pr nto === nto * ndat_left
        end
        dim_indexes = [0,1]
        dims = [ndat_left, @dims[@idim]]
        print_call.call
      }, else: lambda {
        pr ndat_left === 1
        pr ndat_right === 1
        pr For( i, 0, @idim - 1 ) {
          if @ld and util != :cost then
            pr ndat_left === ndat_left * @nx[i]
          else
            pr ndat_left === ndat_left * @dims[i]
          end
          pr ndat_left === ndat_left * 2 if @wavelet
        }
        pr For( i, @idim + 1, @ndim - 1 ) {
          if @ld and util != :cost then
            pr ndat_right === ndat_right * @nx[i]
          else
            pr ndat_right === ndat_right * @dims[i]
          end
          pr ndat_right === ndat_right * 2 if @wavelet
        }
        if @narr then
          pr nti === nti * ndat_left * ndat_right
          pr nto === nto * ndat_left * ndat_right
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
    @vars.push @ndim  = Int( "d",  :dir => :in )
    @vars.push @dims  = Int( "n",  :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    @vars.push @bc    = Int( "bc", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    if @ld then
      @vars.push @nx  = Int( "nx", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
      @vars.push @ny = Int( "ny", :dir => :in, :dim => [ Dim(0, @ndim - 1) ] )
    end
    @vars.push @narr     = Int( "narr", :dir => :in ) if @narr
    @vars.push @x     = Real("x",  :dir => :in,  :restrict => true, :dim => [ Dim() ] )
    @vars.push @y     = Real("y",  :dir => options[:a_y] ? :inout : :out, :restrict => true, :dim => [ Dim() ] )
    @vars.push @w1 = Real("w1", :dir => :inout, :restrict => true, :dim => [ Dim() ] ) if options[:work]
    @vars.push @w2 = Real("w2", :dir => :inout, :restrict => true, :dim => [ Dim() ] ) if options[:work]
    @vars.push @a = Real("a",:dir => :in,:dim => [ Dim(0, @ndim - 1)]) if options[:a]
    @vars.push @a_x = Real("a_x",:dir => :in) if options[:a_x]
    @vars.push @a_y = Real("a_y",:dir => :in) if options[:a_y]
    @vars.push @dot_in = Real("dot_in",:dir => :out,:dim =>[ Dim(0, @ndim - 1)]) if options[:dot_in]

    @transpose = 0
    @transpose = options[:transpose] if options[:transpose]
    @options=options
    @procs = {}
    @needed_subops = {}
    opts = @options.dup
    opts.delete(:a)
    opts.delete(:a_x)
    opts.delete(:a_y)
    opts.delete(:zero_out_work)
    dim_indexes_a = []
    if @transpose == 0 then
      dim_indexes_a = [ [1, 0], [0, 1], [2, 0, 1] ]
    elsif @transpose == 1
      dim_indexes_a = [ [1, 0] ]
    elsif @transpose == -1
      dim_indexes_a = [ [0, 1] ]
    end
    opt_base = []
    init_options = {}
    init_options[:a_x] = @options[:a_x] if @options[:a_x]
    init_options[:zero_out_work] = @options[:zero_out_work] if @options[:zero_out_work]
    opt_base.push( init_options ) if init_options.length > 0
    opt_base.push( { :a => @options[:a] } ) if @options[:a]
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
    function_name += "d_" if default_real_size == 8
    function_name += "s_" if default_real_size == 4
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
    p = Procedure(function_name,@vars) {
      dims_actual = Int( "dims_actual", :allocate => true, :dim => [ Dim(0, @ndim - 1) ] ) if get_lang == FORTRAN
      dims_actual = Int( "dims_actual", :allocate => true, :dim => [ Dim(0, 16) ] ) if get_lang == C
      dims_left   = Int "dims_left"
      ni = Int "ni"
      ndat = Int "ndat"
      ndat2 = Int "ndat2"
      ndat_tot_in = Int "nti" if @narr
      ndat_tot_out = Int "nto" if @narr
      i = Int "i"
      j = Int "j"
      decl i, j, dims_actual, dims_left, ni, ndat, ndat2
      decl ndat_tot_in, ndat_tot_out if @narr
      dims = []
      dim_indexes = []
      pr dims_left === @ndim
      pr For( i, 0, @ndim - 1 ) {
        if @ld then
          if @wavelet then
            pr dims_actual[i] === @nx[i] * 2
          else
            pr dims_actual[i] === @nx[i]
          end
        else
          pr If(@bc[i] == BC::SHRINK => lambda {
            if @wavelet then
              pr dims_actual[i] === @dims[i] * 2 + @filter.length - 2
            else
              pr dims_actual[i] === @dims[i] + @filter.length - 1
            end
          }, else: lambda {
            if @wavelet then
              pr dims_actual[i] === @dims[i] * 2
            else
              pr dims_actual[i] === @dims[i]
            end
          })
        end
      }
      compute_ni_ndat = lambda { |indx|
        indx = @ndim - 1 - indx if @transpose == -1
        pr ni === @dims[indx]
        pr ndat === 1
        pr For(j, 0, @ndim - 1) {
          pr If( j != indx ) {
            pr ndat === ndat * dims_actual[j]
          }
        }
        d = [ ni, ndat ]
        d_indexes = [ 1, 0]
        d.reverse! if @transpose == -1
        d_indexes.reverse! if @transpose == -1
        return [ d, d_indexes ]
      }
      compute_ndat_ni_ndat2 = lambda { |indx|
        pr ni === @dims[indx]
        pr ndat === 1
        pr ndat2 === 1
        pr For(j, 0, @ndim - 1) {
          pr If( j < indx => lambda {
            pr ndat === ndat * dims_actual[j]
          }, j > indx => lambda {
            pr ndat2 === ndat2 * dims_actual[j]
          })
        }
        return [ [ndat, ni, ndat2], [2, 0, 1] ]
      }

      opts = @options.dup
      opts.delete(:a_x)
      opts.delete(:a_y)
      opts.delete(:zero_out_work)

      print_call = lambda { |indx, init, last, datas, multi_conv|
        vars = dims
        if multi_conv then
          if dim_indexes.length == 2 then
            pr ndat_tot_in === dims[dim_indexes[0]]
          else
            pr ndat_tot_in === dims[dim_indexes[0]] * dims[dim_indexes[1]]
          end
          pr ndat_tot_out === ndat_tot_in
          pr ndat_tot_in === ndat_tot_in * dims_actual[indx]
          if @ld then
            if @wavelet then
              pr ndat_tot_out === ndat_tot_out * @ny[indx] * 2
            else
              pr ndat_tot_out === ndat_tot_out * @ny[indx]
            end
          end
          f = For(j, 0, @narr-1)
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
        if init and @options[:zero_out_work]
          opt[:zero_out_work] = @options[:zero_out_work]
        end
        if last and @options[:a_y] then
          vars2.push( @a_y )
          opt[:a_y] = @options[:a_y]
        end
        vars2.push( @dot_in[indx] ) if @options[:dot_in]
        dats = []
        opt.update( opts )
        case_args = {
          BC::PERIODIC => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::PERIODIC), dim_indexes, opt).base_name

            if multi_conv then
              pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            pr dims_actual[indx] === @ny[indx]  if @ld
          },
          BC::GROW => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::GROW), dim_indexes, opt).base_name

            if multi_conv then
              if not @ld then
                if @wavelet then
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 2 )
                else
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] + @filter.length - 1 )
                end
              end
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            if @ld then
              pr dims_actual[indx] === @ny[indx]  if @ld
            else
              if @wavelet then
                pr dims_actual[indx] === dims_actual[indx] + @filter.length - 2
              else
                pr dims_actual[indx] === dims_actual[indx] + @filter.length - 1
              end
            end
          },
          BC::SHRINK => lambda {
            procname = ConvolutionOperator1d::new(@filter, BC::new(BC::SHRINK), dim_indexes, opt).base_name

            if multi_conv then
              if not @ld then
                if @wavelet then
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 2 )
                else
                  pr ndat_tot_out === ndat_tot_out * ( dims_actual[indx] - @filter.length + 1 )
                end
              end
              dats[0] = (datas[0][ndat_tot_in*j+1]).address
              dats[1] = (datas[1][ndat_tot_out*j+1]).address
              pr f
            else
              dats = datas.dup
            end
            pr @procs[procname].call( *vars, *dats, *vars2 )
            close f if multi_conv
            if @ld then
              pr dims_actual[indx] === @ny[indx]  if @ld
            else
              if @wavelet then
                pr dims_actual[indx] === dims_actual[indx] - @filter.length + 2
              else
                pr dims_actual[indx] === dims_actual[indx] - @filter.length + 1
              end
            end
          }
        }
        case_args[BC::NPERIODIC] = lambda {
          procname = ConvolutionOperator1d::new(@filter, BC::new(BC::NPERIODIC), dim_indexes, opt).base_name

          if multi_conv then
            pr ndat_tot_out === ndat_tot_out * dims_actual[indx] if not @ld
            dats[0] = (datas[0][ndat_tot_in*j+1]).address
            dats[1] = (datas[1][ndat_tot_out*j+1]).address
            pr f
          else
            dats = datas.dup
          end
          pr @procs[procname].call( *vars, *dats, *vars2 )
          close f if multi_conv
          pr dims_actual[indx] === @ny[indx]  if @ld
        } if @poisson
        pr Case( @bc[indx], case_args)
      }

      pr If( @ndim == 1 => lambda {
        conv_number = 1
        conv_number = @narr if @narr
        dims = [ @dims[0], conv_number ]
        dim_indexes = [ 1, 0 ]
        dims.reverse! if @transpose == -1
        dim_indexes.reverse! if @transpose == -1
        datas = [ @x, @y ]
        print_call.call( 0, true, true, datas, false )
      }, else: lambda {
        dims, dim_indexes = compute_ni_ndat.call( 0 )
        datas = [ @x, @w1 ]
        datas = [ @x, @y ] if not @options[:work]
        print_call.call( 0, true, false, datas, @narr )
        pr dims_left === dims_left - 1
        pr i === 1
        pr While( dims_left > 2 ) {
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
          pr i === i + 1
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
          pr i === i + 1
          pr dims_left === dims_left - 2
        }
        pr If( dims_left == 2 => lambda {
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
          pr i === i + 1
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
        }, else: lambda {
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
  if get_lang == C then
    get_output.puts "#include <immintrin.h>" if get_architecture == X86
    get_output.puts "#include <arm_neon.h>" if get_architecture == ARM
    if macro then
      get_output.print "#define modulo(a, b) ((a+b)%(b))\n"
      get_output.print "#define min( a, b) ((a) < (b) ? a : b)\n"
      get_output.print "#define max( a, b) ((a) > (b) ? a : b)\n"
    else
      get_output.print "static inline #{Int::new.decl} modulo( #{Int::new.decl} a, #{Int::new.decl} b) { return (a+b)%b;}\n"
      get_output.print "static inline #{Int::new.decl} min( #{Int::new.decl} a, #{Int::new.decl} b) { return a < b ? a : b;}\n"
      get_output.print "static inline #{Int::new.decl} max( #{Int::new.decl} a, #{Int::new.decl} b) { return a > b ? a : b;}\n"
    end
  end
end


