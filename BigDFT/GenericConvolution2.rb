require 'BOAST'
require 'narray'
include BOAST

register_funccall("modulo")
register_funccall("min")
register_funccall("max")

class Filter
  # List of the floating point values of the convolution filter
  attr_reader :fil_array
  # central point of the filter
  attr_reader :center
  attr_reader :length
  # name of the filter
  attr_reader :name
  attr_reader :base_name

  def elem_number
    n = (@unroll + @vec_len - 1)/@vec_len
    return n <= 0 ? 1 : n
  end

  def init_optims(tt_arr, unroll, vec_len)
    @tt_arr = tt_arr
    @unroll = unroll
    @vec_len = vec_len
    @tt_n = elem_number
    @tt = init_tt
    @filter_val = init_filter_val
  end

  def get_tt
    return @tt
  end

  def get_filter_val
    return @filter_val
  end

  def decl_filter_val
    decl *([@filter_val].flatten)
  end

end

class ConvolutionFilter < Filter
  # BOAST variables
  # Filter array (to be used on BOAST functions)
  attr_reader :fil
  # extremes of the filter, calculated via its central point (integers)
  attr_reader :lowfil_val, :upfil_val
  # extremes of the filter, calculated via its central point (BOAST object)
  attr_reader :lowfil, :upfil

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
    @base_name = "s0s0"
    @length = @fil_array.length
  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return Real("tt", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len)
    else
      return (0..(@tt_n-1)).collect { |index| Real("tt#{index}", :vector_length => @vec_len) }
    end
  end

  def set_zero_tt( tt_ind )
    pr @tt[tt_ind].set( 0.0 )
  end

  def init_filter_val
    return @tt[0].copy("filter_val")
  end

  def set_filter_val( index, options = {} )
    pr @filter_val === Set(@fil[index], @tt[0])
  end

  def compute_values( tt_ind, indexes, data_in)
    out = @tt[tt_ind]
    pr out === FMA(Load(data_in[*indexes].set_align(out.type.total_size), out), @filter_val, out) # + @x[*i_in]*@filter.fil[l]
  end

  def post_process_and_store_values( tt_ind, indexes, data_out, transpose, options = {} )
    i_out = indexes.rotate(transpose)
    i_in = options[:indexes_in]
    out = @tt[tt_ind]
    pr out === out * Set(options[:a], out) if options[:a]

    finish_block = lambda {
      pr options[:dot_in_tmp] === FMA(Load(options[:data_in][*i_in].set_align(out.type.total_size), out), out, options[:dot_in_tmp]) if options[:dot_in_tmp] #reduction !
      pr out === FMA(Load(options[:data_in][*i_in], out), Set(options[:a_x], out), out) if options[:a_x]
    }
    if options[:bc].grow and (options[:dot_in_tmp] or options[:a_x]) and options[:side] != :center then
      if options[:side] == :begin then
        pr If(options[:position] >= 0, &finish_block)
      elsif options[:side] == :end then
        pr If(options[:position] < options[:dim], &finish_block)
      end
    else
      finish_block.call
    end

    #to be controlled in the case of non-orthorhombic cells for kinetic operations
    pr out === out + Load(options[:data_in2][*i_in], out) if options[:data_in2]

    if options[:accumulate] or (options[:kinetic] == :inplace and not options[:zero_out])  then
      pr out === out + Load(data_out[*i_out].set_align(out.type.total_size), out)
    elsif options[:a_y] then
      pr out === FMA(Set(options[:a_y], out), Load(data_out[*i_out].set_align(out.type.total_size), out), out)
    end
    pr data_out[*i_out].set_align(out.type.total_size) === out
    pr options[:data_out2][*i_out] === Load(options[:data_in][*i_in], out) if options[:data_out2]
  end

  def get_input_dim( dim )
    return Dim(0, dim - 1)
  end

  alias get_output_dim get_input_dim

  def get_input_dim_shrink( dim )
    return Dim( lowfil, dim + upfil  - 1)
  end

  def get_output_dim_grow( dim )
    return Dim(-upfil,  dim - lowfil - 1)
  end

  def get_input_dim_ld_shrink( dim )
    return Dim(lowfil, dim + lowfil - 1)
  end

  def get_output_dim_ld_grow( dim )
    return Dim( -upfil, dim - upfil - 1)
  end

  def decl_filters( options = {} )
    decl fil
  end

  def cost
    return @length * 2
  end

end #class ConvolutionFilter

class PoissonFilter < ConvolutionFilter
  attr_reader :filters_array
  attr_reader :filters

  def initialize(name, filts, nord)
    @center = nord/2
    tmp_array = []
    @filters=[]
    index=0
    filts.each_with_index { |e,i|
      @filters[i]=ConvolutionFilter::new(name+"_#{i}", e, @center)
        e.each{|val|
            tmp_array[index]=val
            index= index+1
        }    
    }
    arr = ConstArray::new(tmp_array,Real)
    @filters_array = Real("#{name}_fil",:constant => arr,:dim => [ Dim(0,(2*@center+1)*(2*@center+1)-1) ])

    super(filters[@center].name, @filters[@center].fil_array, @center)
    @name = name

    @fil_array = filts.dup

  end

  def decl_filters( options = {} )
    decl filters_array if options[:bc].nper
    decl @filters[@center].fil
  end

  def set_filter_val( index, options = {} )
    if options[:bc].nper then
      if options[:side] == :center then
        super
      elsif options[:side] == :begin then
        pr @filter_val === Set(@filters_array[((options[:position] - (upfil)                       )*(2*@center+1)) + (@center)*(2*@center+1) + index + @center], @tt[0])
      elsif options[:side] == :end then
        pr @filter_val === Set(@filters_array[((options[:position] + (@center) - options[:dim] + 1 )*(2*@center+1)) + (@center)*(2*@center+1) + index + @center], @tt[0])
      end
    else
      super
    end
  end

end

class WaveletFilter < Filter
  attr_reader :low, :high
  attr_reader :low_even, :low_odd
  attr_reader :high_even, :high_odd
  attr_reader :convolution
  def initialize(name, filt1, opts = {})
    @convolution = opts[:convolution_filter]
    @fil_array = filt1.dup
    if @convolution then
      require 'bigdecimal'
      lp = @fil_array.collect { |e| BigDecimal.new(e) }
      hp = lp.each_with_index.collect { |e,i| i % 2 == 0 ? e : -e }.reverse
      cv = @convolution.fil_array.collect { |e| BigDecimal.new(e) }

      lp_length = lp.length
      hp_length = hp.length
      cv_length = @convolution.length

      lp_center = lp_length / 2
      hp_center = hp_length / 2
      cv_center = @convolution.center

      cv_bounds = [ -cv_center, cv_length - cv_center - 1 ]
      lp_bounds = [ -lp_center, lp_length - lp_center - 1 ]
      hp_bounds = [ -hp_center, hp_length - hp_center - 1 ]

      lp_cv_bounds = [ cv_bounds[0] - lp_bounds[1], cv_bounds[1] - lp_bounds[0] ]
      hp_cv_bounds = [ cv_bounds[0] - hp_bounds[1], cv_bounds[1] - hp_bounds[0] ]

      lp_cv_bounds[0] -= 1 if lp_cv_bounds[0] % 2 != 0
      hp_cv_bounds[0] -= 1 if hp_cv_bounds[0] % 2 != 0
      lp_cv_bounds[1] += 1 if lp_cv_bounds[1] % 2 == 0
      hp_cv_bounds[1] += 1 if hp_cv_bounds[1] % 2 == 0

      lp_cv_center = - lp_cv_bounds[0]
      hp_cv_center = - hp_cv_bounds[0]

      lp_cv = (lp_cv_bounds[0]..lp_cv_bounds[1]).collect { |i|
        sum = BigDecimal.new(0)
        (lp_bounds[0]..lp_bounds[1]).each { |j|
          sum += lp[j+lp_center] * cv[i-j-1+cv_center] if (0...cv_length).include?(i-j-1+cv_center)
        }
        sum
      }

      hp_cv = (hp_cv_bounds[0]..hp_cv_bounds[1]).collect { |i|
        sum = BigDecimal.new(0)
        (hp_bounds[0]..hp_bounds[1]).each { |j|
          sum += hp[j+hp_center] * cv[i-j-1+cv_center] if (0...cv_length).include?(i-j-1+cv_center)
        }
        sum
      }
      @center = lp_cv_center
      @low = ConvolutionFilter::new(name+"_l", lp_cv.collect { |e| e.truncate(30).to_s }, lp_cv_center)
      @high = ConvolutionFilter::new(name+"_h", hp_cv.collect { |e| e.truncate(30).to_s }, hp_cv_center)
    else
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

    end

    @length = @low.fil_array.length
    @name = name
  end

  def lowfil(wavelet=nil)
    return @low_even.lowfil
  end

  def upfil(wavelet=nil)
    return @low_even.upfil
  end

  def init_filter_val
    return (0..3).collect { |index| @tt[0][0].copy("filter_val#{index}") }
  end

  def compute_values( tt_ind, indexes, data)
    out_even = @tt[0][tt_ind]
    out_odd  = @tt[1][tt_ind]
    i_in0 = indexes[0].flatten
    i_in1 = indexes[1].flatten
    pr out_even === FMA(Load(data[*i_in0], out_even), @filter_val[0], out_even)
    pr out_odd  === FMA(Load(data[*i_in0], out_odd ), @filter_val[1], out_odd )
    pr out_even === FMA(Load(data[*i_in1], out_even), @filter_val[2], out_even)
    pr out_odd  === FMA(Load(data[*i_in1], out_odd ), @filter_val[3], out_odd )
  end

  def post_process_and_store_values( tt_ind, indexes, data, transpose, options = {} )
    out_even = @tt[0][tt_ind]
    out_odd  = @tt[1][tt_ind]
    i_out0 = indexes[0].rotate(transpose).flatten
    i_out1 = indexes[1].rotate(transpose).flatten
    if options[:a] then
      pr out_even === out_even * Set( options[:a], out_even)
      pr out_odd  === out_odd  * Set( options[:a], out_odd )
    end

    if options[:accumulate] then
      pr out_even === Load(data[*i_out0], out_even) + out_even
      pr out_odd  === Load(data[*i_out1], out_odd ) + out_odd
    elsif options[:a_y] then
      pr out_even === FMA(Load(data[*i_out0], out_even), Set(options[:a_y], out_even), out_even)
      pr out_odd  === FMA(Load(data[*i_out1], out_odd ), Set(options[:a_y], out_odd ), out_odd )
    end

    pr data[*i_out0] === out_even
    pr data[*i_out1] === out_odd
  end

  def set_zero_tt( tt_ind )
    pr @tt[0][tt_ind].set( 0.0 )
    pr @tt[1][tt_ind].set( 0.0 )
  end

  def decl_filters( options = {} )
    decl low_even.fil
    decl low_odd.fil
    decl high_even.fil
    decl high_odd.fil
  end

  def cost
    return @length * 2 * 2
  end

end

class WaveletFilterDecompose < WaveletFilter

  def initialize(name, filt1, opts = {})
    super

    @base_name = "s0s1"

    center_half = @center / 2

    filt_1 = @low.fil_array.values_at(*(0..(@low.fil_array.length-1)).step(2).collect)
    @low_even = ConvolutionFilter::new(name+"_le", filt_1, center_half)

    filt_2 = @low.fil_array.values_at(*(1..(@low.fil_array.length-1)).step(2).collect)
    @low_odd = ConvolutionFilter::new(name+"_lo", filt_2, center_half)

    filt_3 = @high.fil_array.values_at(*(0..(@high.fil_array.length-1)).step(2).collect)
    @high_even = ConvolutionFilter::new(name+"_he", filt_3, center_half)

    filt_4 = @high.fil_array.values_at(*(1..(@high.fil_array.length-1)).step(2).collect)
    @high_odd = ConvolutionFilter::new(name+"_ho", filt_4, center_half)

  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return [ Real("lt", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len),
               Real("ht", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len) ]
    else
      return [ (0..(@tt_n-1)).collect { |index| Real("lt#{index}", :vector_length => @vec_len) },
               (0..(@tt_n-1)).collect { |index| Real("ht#{index}", :vector_length => @vec_len) } ]
    end
  end

  def set_filter_val( index, options = {} )
    pr @filter_val[0] === Set(@low_even.fil[index], @tt[0][0])
    pr @filter_val[1] === Set(@high_even.fil[index], @tt[0][0])
    pr @filter_val[2] === Set(@low_odd.fil[index], @tt[0][0])
    pr @filter_val[3] === Set(@high_odd.fil[index], @tt[0][0])
  end

  def get_input_dim( dim )
    return [ Dim(0, 1), Dim( 0, dim - 1 ) ]
  end

  def get_input_dim_shrink( dim )
    return [ Dim(0, 1), Dim( lowfil, dim + upfil - 1 ) ]
  end

  def get_input_dim_ld_shrink( dim )
    return [ Dim(0, 1), Dim( lowfil, dim + lowfil - 1 ) ]
  end

  def get_output_dim( dim )
    return [ Dim( 0, dim - 1 ), Dim(0, 1) ]
  end

  def get_output_dim_grow( dim )
    return [ Dim( -upfil, dim - lowfil - 1 ), Dim(0, 1) ]
  end

  def get_output_dim_ld_grow( dim )
    return [ Dim( -upfil, dim - upfil - 1 ), Dim(0, 1) ]
  end

end

class WaveletFilterRecompose < WaveletFilter

  def initialize(name, filt1, opts = {})
    super

    @base_name = "s1s0"

    center_half = (@low.fil_array.length - @center - 1)/2

    filt_1 = @low.fil_array.reverse.values_at(*(0..(@low.fil_array.length-1)).step(2).collect)
    @low_even = ConvolutionFilter::new(name+"_lre", filt_1, center_half)

    filt_2 = @low.fil_array.reverse.values_at(*(1..(@low.fil_array.length-1)).step(2).collect)
    @low_odd = ConvolutionFilter::new(name+"_lro", filt_2, center_half)

    filt_3 = @high.fil_array.reverse.values_at(*(0..(@high.fil_array.length-1)).step(2).collect)
    @high_even = ConvolutionFilter::new(name+"_hre", filt_3, center_half)

    filt_4 = @high.fil_array.reverse.values_at(*(1..(@high.fil_array.length-1)).step(2).collect)
    @high_odd = ConvolutionFilter::new(name+"_hro", filt_4, center_half)

  end

  def init_tt
    #try to modify tt scalars into arrays of size unroll
    if @tt_arr then
      return [ Real("et", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len),
               Real("ot", :dim => [ Dim(0,@tt_n-1)], :allocate => true, :vector_length => @vec_len) ]
    else
      return [ (0..(@tt_n-1)).collect { |index| Real("et#{index}", :vector_length => @vec_len) },
               (0..(@tt_n-1)).collect { |index| Real("ot#{index}", :vector_length => @vec_len) } ]
    end
  end

  def set_filter_val( index, options = {} )
    pr @filter_val[0] === Set(@low_odd.fil[index], @tt[0][0])
    pr @filter_val[1] === Set(@low_even.fil[index], @tt[0][0])
    pr @filter_val[2] === Set(@high_odd.fil[index], @tt[0][0])
    pr @filter_val[3] === Set(@high_even.fil[index], @tt[0][0])
  end

  def get_input_dim( dim )
    return [ Dim( 0, dim - 1 ), Dim(0, 1) ]
  end

  def get_input_dim_shrink( dim )
    return [ Dim( lowfil, dim + upfil - 1 ), Dim(0, 1) ]
  end

  def get_input_dim_ld_shrink( dim )
    return [ Dim( lowfil, dim + lowfil - 1 ), Dim(0, 1) ]
  end

  def get_output_dim( dim )
    return [ Dim(0, 1), Dim( 0, dim - 1 ) ]
  end

  def get_output_dim_grow( dim )
    return [ Dim(0, 1), Dim( -upfil, dim - lowfil - 1 ) ]
  end

  def get_output_dim_ld_grow( dim )
    return [ Dim(0, 1), Dim( -upfil, dim - upfil - 1 ) ]
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
    @filter = filter.dup
    @bc = bc
    @transpose = options[:transpose]
    @dim_indexes = dim_indexes
    @ld = options[:ld]
    @kinetic = options[:kinetic]
    @wavelet = options[:wavelet]
    @poisson = options[:poisson]

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
    @base_name += @filter.base_name + "_"
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

    @dim_ngs = @filter.get_input_dim_shrink( @dim_n )
    @dim_nsg = @filter.get_output_dim_grow( @dim_n )

  end

  def compute_dimx_dimy
    dimx = @dims_in.collect{ |dim|
      if not dim.name.match("ndat") and @bc.shrink and not @ld then
        @filter.get_input_dim_shrink( dim )
      elsif not dim.name.match("ndat") and @bc.shrink
        @filter.get_input_dim_ld_shrink( dim )
      elsif not dim.name.match("ndat")
        @filter.get_input_dim( dim )
      else
        Dim(0, dim - 1)
      end
    }
    dimx.flatten!
    dimy = @dims_out.collect{ |dim|
      if not dim.name.match("ndat") and @bc.grow and not @ld then
        @filter.get_output_dim_grow( dim )
      elsif not dim.name.match("ndat") and @bc.grow then
        @filter.get_output_dim_ld_grow( dim )
      elsif not dim.name.match("ndat")
        @filter.get_output_dim( dim )
      else
        Dim(0, dim - 1)
      end
    }
    dimy.rotate!(@transpose).flatten!
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
    vars.push(0.0) if @dot_in
    return vars
  end

  def optimize(opt_space)
    opt_space=GenericOptimization::new if not opt_space
    t_best=Float::INFINITY
    p_best_optim = nil
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
        p_best_optim = optim
      end
    }
    return self.procedure(p_best_optim)
  end

  def cost(*dimens)
    n = dimens[@dim_indexes.last]
    ndat = 1
    @dim_indexes[0...-1].each { |indx|
      ndat *= dimens[indx]
    }
    return n * @filter.cost * ndat
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

  def stos(bool)
   if bool then
    return 't'
   else
    return 'f'
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

    vec_len = options[:vector_length] if not @dot_in and unrolled_dim == 0 and options[:vector_length] and options[:vector_length] <= @filter.length and unroll % options[:vector_length] == 0

    mod_arr = false if @bc.free
    util = options[:util]

    function_name = @base_name   + 
      "_u#{unroll}_v#{vec_len}_#{unrolled_dim}" ###_#{mod_arr}_#{tt_arr}_#{unroll_inner}"
    function_name += '_'+stos(mod_arr)+'_'+stos(tt_arr)+'_'+stos(unroll_inner)

    function_name += "_" + util.to_s if util

    if util == :cost then
      return Procedure(function_name, @dims + [@cost] ){
        decl ndat_t = Int("ndat_t")
        pr ndat_t === 1
        @dim_indexes[0...-1].each { |indx|
          pr ndat_t === ndat_t * @dims[indx]
        }
        pr @cost === @dims[@dim_indexes.last] * @filter.cost * ndat_t
      }
    end

    l = Int("l")
    @filter.init_optims(tt_arr, unroll, vec_len)
    tt = @filter.get_tt
    iters =  (1..@dims.length).collect{ |index| Int("i#{index}")}
    mods = get_mods(mod_arr)
    constants = get_constants

    return Procedure(function_name, vars, :constants => constants ){
      @filter.decl_filters( :bc => @bc )
      decl *iters
      decl l
      decl *([tt].flatten)
      @filter.decl_filter_val
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
      pr OpenMP::Parallel(default: :shared, reduction: (@options[:dot_in] ? {"+" => @dot_in} : nil ), private: iters + [l] + [tt] + ( @filter.get_filter_val ? [@filter.get_filter_val].flatten : [] )) { 
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

  def init_values(side, iters, l, t, tlen, unro, mods)
    vec_len = [t].flatten[0].type.vector_length
    (0...tlen).step(vec_len).each_with_index{ |_,tt_ind|
      #WARNING: the eks conditional here can be relaxed
      @filter.set_zero_tt(tt_ind)
    }
  end

  def compute_values(side, iters, l, tlen, unro, mods)
    processed_dim = @dim_indexes[-1]
    vec_len = [@filter.get_tt].flatten[0].type.vector_length

    (0...tlen).step(vec_len).each_with_index{ |ind,tt_ind|

      if @bc.free or (side == :center) then
        i_in = input_index(unro, iters, ind, processed_dim, l, nil, nil, side)
      elsif mods then
        i_in = input_index(unro, iters, ind, processed_dim, l, nil, mods, side) 
      else
        i_in = input_index(unro, iters, ind, processed_dim, l, @dims[processed_dim],nil, side)
      end

      @filter.compute_values(tt_ind, i_in, @x)
    }
  end

  
  def post_process_and_store_values(side, iters, l, t, tlen, unro, mods)
    processed_dim = @dim_indexes[-1]
    vec_len = [t].flatten[0].type.vector_length
    (0...tlen).step(vec_len).each_with_index{ |ind,tt_ind|
      i_out = output_index(unro, iters, ind)
      i_in = input_index(unro, iters, ind)
      @filter.post_process_and_store_values( tt_ind,
                                             i_out,
                                             @y,
                                             @transpose,
                                             :bc => @bc,
                                             :side => side,
                                             :position => iters[processed_dim],
                                             :dim => @dims[processed_dim],
                                             :accumulate => @accumulate,
                                             :a => @a,
                                             :a_y => @a_y,
                                             :a_x => @a_x,
                                             :dot_in_tmp => @dot_in_tmp,
                                             :kinetic => @kinetic,
                                             :zero_out => @options[:zero_out],
                                             :data_in => @x,
                                             :data_in2 => @x2,
                                             :indexes_in => i_in,
                                             :data_out2 => @y2 )
    }
  end

  def for_conv(side, iters, l, t, tlen, unro, mods, unroll_inner)

    init_values(side, iters, l, t, tlen, unro, mods)

    loop_start, loop_end = get_loop_start_end( side, iters)

    pr For( l, loop_start, loop_end, :unroll => unroll_inner ) {
      @filter.set_filter_val(l, :side => side, :position => iters[@dim_indexes[-1]], :dim => @dims[@dim_indexes[-1]], :bc => @bc)
      compute_values(side, iters, l, tlen, unro, mods)
    }

    post_process_and_store_values(side, iters, l, t, tlen, unro, mods)

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
    function_name += @filter.base_name
    function_name += "_1d_"
    function_name += @filter.name
    function_name += "_dotin" if @dot_in
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


