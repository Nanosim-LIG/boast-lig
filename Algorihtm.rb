module ConvolutionGenerator
  FORTRAN = 1
  C = 2
  $output = STDOUT
  $lang = FORTRAN
  $replace_constants = true
  $default_int_size = 4
  $default_real_size = 8
  $indent_level = 0
  $indent_increment = 2

  def ConvolutionGenerator::set_lang(lang)
    $lang = lang
  end

  def set_output(output)
    $output = output
  end

  class Expression
    attr_reader :operator
    attr_reader :operand1
    attr_reader :operand2
    def initialize(operator, operand1, operand2)
      @operator = operator
      @operand1 = operand1
      @operand2 = operand2
    end
    def to_s
      self.to_str
    end

    def ===(x)
      return Expression::new("=",self,x)
    end
 
    def +(x)
      return Expression::new("+",self,x)
    end
 
    def *(x)
      return Expression::new("*",self,x)
    end
 
    def -(x)
      return Expression::new("-",self,x)
    end

    def to_str
      s = ""
      s += @operand1.to_s if @operand1
      s += @operator.to_s
      s += @operand2.to_s
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang==C
      $output.puts s if final
      return s
    end
  end

  class Mult < Expression
    def initialize(operand1, operand2)
      super("*", operand1, operand2)
    end
  end

  class Add < Expression
    def initialize(operand1, operand2)
      super("+", operand1, operand2)
    end
  end

  class Affect < Expression
    def initialize(operand1, operand2)
      super("=", operand1, operand2)
    end
  end

  class Index < Expression
    attr_reader :source
    attr_reader :indexes
    def initialize(source, indexes)
      @source = source
      @indexes = indexes
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C
    end
    def to_str_fortran
      s = ""
      s += "#{@source}(#{@indexes.join(", ")})"
      return s
    end
    def to_str_c
      s = ""
      s += "#{@source}["
      dim = @source.dimension.first
      if dim.val2 then
        start = dim.val1
      else
        start = 1
      end
      s += "#{@indexes.first}-#{start}"
      i=1
      ss = ""
      @source.dimension[0..-2].each{ |d|
        if d.val2 then
          ss += "*(#{d.val2}-#{d.val1}+1)"
        else
          ss += "*#{d.val1}"
        end
        dim = @source.dimension[i]
        if dim.val2 then
          start = dim.val1
        else
          start = 1
        end
        s += "+(#{@indexes[i]}-#{start})"+ss
        i+=1
      }
      s += "]"
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C
      $output.puts s if final
      return s
    end
  end

  class Variable
    attr_reader :name
    attr_reader :direction
    attr_accessor :constant
    attr_reader :type
    attr_reader :dimension
    def initialize(name,type,hash={})
      @name = name
      @direction = hash[:direction]
      @constant = hash[:constant]
      @dimension = hash[:dimension]
      @type = type::new(hash)
    end
  
    def to_s
      self.to_str
    end    

    def to_str
      return @constant.to_s if @constant and $replace_constants and not @dimension
      return @name.to_str
    end

    def ===(x)
      return Expression::new("=",self,x)
    end
 
    def +(x)
      return Expression::new("+",self,x)
    end
 
    def *(x)
      return Expression::new("*",self,x)
    end
 
    def -(x)
      return Expression::new("-",self,x)
    end
 
    def -@
      return Expression::new("-",nil,self)
    end

    def [](*args)
      return Index::new(self,args)
    end
 
    def indent
       return " "*$indent_level
    end

    def finalize
       s = ""
       s += ";" if $lang == C
       s+="\n"
       return s
    end

    def decl(final=true)
      return self.decl_fortran(final) if $lang == FORTRAN
      return self.decl_c(final) if $lang == C
    end

    def decl_c(final=true)
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += @type.decl
      if(@dimension and not @constant) then
        s += " *"
      end
      s += " #{@name}"
      if(@dimension and @constant) then
        s += "[]"
      end
      s += " = #{@constant}" if @constant
      s += self.finalize if final
      $output.print s if final
      return s
    end

    def decl_fortran(final=true)
      s = ""
      s += self.indent if final
      s += @type.decl
      s += ", intent(#{@direction})" if @direction 
      s += ", parameter" if @constant
      if(@dimension) then
        s += ", dimension("
        s += @dimension[0]
        @dimension[1..-1].each { |d|
           s += ", "
           s += d
        }
        s += ")"
      end
      s += " :: #{@name}"
      s += " = #{@constant}" if @constant
      s += self.finalize if final
      $output.print s if final
      return s
    end

  end

  class Real
    attr_reader :size
    def initialize(hash={})
      if hash["size"] then
        @size = hash["size"]
      else
        @size = $default_real_size
      end
    end
    def decl
      return "real(kind=#{@size})" if $lang == FORTRAN
      if $lang == C then
        return "float" if @size == 4
        return "double" if @size == 8
      end
    end
  end

  class Procedure
    attr_reader :name
    attr_reader :parameters
    attr_reader :constants
    def initialize(name, parameters=[], constants=[], &block)
      @name = name
      @parameters = parameters
      @constants = constants
      @block = block
    end
    def decl(final=true)
      return self.decl_fortran(final) if $lang==FORTRAN
      return self.decl_c(final) if $lang==C
    end
    def close(final=true)
      return self.close_fortran(final) if $lang==FORTRAN
      return self.close_c(final) if $lang==C
    end
    def close_c(final=true)
      $indent_level -= $indent_increment
      s = ""
      s += "}"
      $output.puts s if final
      return s
    end
    def close_fortran(final=true)
      $indent_level -= $indent_increment
      s = ""
      s += "END SUBROUTINE #{name}"
      $output.puts s if final
      return s
    end

    def print(final=true)
      s = self.decl
      if @block then
        @block.call
        s += self.close
      end
      return s
    end

    def decl_c(final=true)
      s = ""
      s += "void #{@name}("
      if parameters.first then
        s += parameters.first.decl(false)
        parameters[1..-1].each { |p|
          s += ", "+p.decl(false)
        }
      end
      s += "){\n"
      $indent_level += $indent_increment
      constants.each { |c|
        s += " "*$indent_level
        s += c.decl(false)
        s += ";\n"
      }
      $output.print s if final
      return s
    end
    def decl_fortran(final=true)
      s = ""
      s += "subroutine #{@name}("
      if parameters.first then
        s += parameters.first
        parameters[1..-1].each { |p|
          s += ", "+p
        }
      end
      s += ")\n"
      $indent_level += $indent_increment
      constants.each { |c|
        s += " "*$indent_level
        s += c.decl(false)
        s += "\n"
      }
      parameters.each { |p|
        s += " "*$indent_level
        s += p.decl(false)
        s += "\n"
      }
      $output.print s if final
      return s
    end
  end

  class Dimension
    attr_reader :val1
    attr_reader :val2
    def initialize(val1,val2=nil)
      @val1 = val1
      @val2 = val2
    end
    def to_str
      s = ""
      if val2 then
        s += val1.to_s
        s += ":"
        s += val2.to_s
      else
        s += val1.to_s
      end
    end
    def to_s
      self.to_str
    end
  end
 
  class Int
    attr_reader :size
    def initialize(hash={})
      if hash["size"] then
        @size = hash["size"]
      else
        @size = $default_int_size
      end
    end
    def decl
      return "integer(kind=#{@size})" if $lang == FORTRAN
      return "int#{8*@size}_t" if $lang == C
    end
  end

  class ConstArray < Array
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C
    end
    def to_str_fortran
      s = ""
      return s if not self.first
      s += "(/ &\n"
      s += self.first 
      self[1..-1].each { |v|
        s += ", &\n"+v
      }
      s += " /)"
    end
    def to_str_c
      s = ""
      return s if not self.first
      s += "{\n"
      s += self.first 
      self[1..-1].each { |v|
        s += ",\n"+v
      }
      s += "}"
    end
  end
 
  class FuncCall
    attr_reader :func_name
    attr_reader :args

    def initialize(func_name, *args)
      @func_name = func_name
      @args = args
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C
    end
    def to_str_fortran
      s = ""
      s += "#{func_name}(#{@args.join(", ")})"
    end
    def to_str_c
      s = ""
      s += "#{func_name}(#{@args.join(", ")})"
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      $output.puts s if final
      return s
    end
  end
 
  class For
    attr_reader :iterator
    attr_reader :begin
    attr_reader :end
    attr_reader :step
    def initialize(i, b, e, s=1, &block)
      @iterator = i
      @begin = b
      @end = e
      @step = s
      @block = block
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C
    end
    def to_str_fortran
      s = ""
      s += "do #{@iterator}=#{@begin}, #{@end}"
      s += ", #{@step}" if @step != 1
      return s
    end
    def to_str_c
      s = ""
      s += "for(#{@iterator}=#{@begin}; #{@iterator}<=#{@end}; #{@iterator}+=#{@step}){"
      return s
    end

    def unroll(final=true)
      raise "Block not given!" if not @block
      if @begin.kind_of?(Variable) then
        start = @begin.constant
      else
        start = @begin.to_i
      end
      if @end.kind_of?(Variable) then
        e = @end.constant
      else
        e = @end.to_i
      end
      if @step.kind_of?(Variable) then
        step = @step.constant
      else
        step = @step.to_i
      end
      raise "Invalid bounds (not constants)!" if not ( start and e and step )
      range = start..e
      range.step(step) { |i|
        @iterator.constant = i
        @block.call
      }
      @iterator.constant = nil
    end

    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      $indent_level += $indent_increment
      $output.puts s if final
      if @block then
        s += "\n"
        @block.call
        s += self.close
      end
      return s
    end
    def close(final=true)
      return self.close_fortran(final) if $lang == FORTRAN
      return self.close_c(final) if $lang == C
    end
    def close_c(final=true)
      s = ""
      $indent_level -= $indent_increment
      s += " "*$indent_level if final
      s += "}"
      $output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      $indent_level -= $indent_increment
      s += " "*$indent_level if final
      s += "enddo"
      $output.puts s if final
      return s
    end
  end

  def ConvolutionGenerator::MagicFilter(filt, center, unroll, invert, free=false )
    function_name = "magicfilter"
    if free then
      function_name += "_free"
    else 
      function_name += "_per"
    end
    if invert then
      function_name += "_inv"
    end
    if unroll>0 then
      function_name += "_u#{unroll}"
    end

    n = Variable::new("n",Int,{:direction => :in})
    ndat = Variable::new("ndat",Int,{:direction => :in})
    if invert then
      lowfil = Variable::new("lowfil",Int,{:constant => 1-center})
      upfil = Variable::new("upfil",Int,{:constant => filt.length-center})
    else
      lowfil = Variable::new("lowfil",Int,{:constant => center-filt.length})
      upfil = Variable::new("upfil",Int,{:constant => center-1})
    end

    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    if free then
      if invert then
        dim_in_min = lowfil
        dim_in_max = n-1+upfil
      else
        dim_out_min = -upfil
        dim_out_max = n-1-lowfil
      end
    end

    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    i = Variable::new("i",Int)
    j = Variable::new("j",Int)
    k = Variable::new("k",Int)
    l = Variable::new("l",Int)
    tt = [Variable::new("tt1",Real)]
    2.upto(unroll) { |index|
      tt.push(Variable::new("tt#{index}",Real))
    }
    if invert then filt = filt.reverse end 
    arr = ConstArray::new(filt)
    fil = Variable::new("fil",Real,{:constant => arr,:dimension => [ Dimension::new(lowfil,upfil) ]})
  
    p = Procedure::new(function_name, [n,ndat,x,y], [lowfil,upfil]) {
  
      i.decl
      j.decl
      k.decl
      l.decl
      tt.each { |t|
        t.decl
      }
  
      fil.decl
  
      if unroll>1 then 
        for1 = For::new(j,1,ndat-(unroll-1), unroll)
      else
        for1 = For::new(j,1,ndat)
      end
      for1.print
      for2 = For::new(i,dim_out_min, dim_out_max) {
        tt.each{ |t|
          (t === 0.0).print
        }
        if free and not invert
          for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) )
        else
          for3 = For::new(l, lowfil, upfil)
        end
        for3.print
        if not free then
          (k === FuncCall::new( "modulo", i+l, n)).print
        else
          (k === i+l).print
        end
        tt.each_index{ |index|
          (tt[index] === tt[index] + x[k,j+index]*fil[l]).print
        }
        for3.close
        tt.each_index{ |index|
          (y[j+index,i] === tt[index]).print
        }
      }.print
      for1.close
  
      if unroll>1 then
        for1 = For::new(j,ndat-FuncCall::new("modulo",ndat,unroll)+1,ndat) {
          for2 = For::new(i,dim_out_min,dim_out_max) {
            (tt[0] === 0.0).print
            if free and not invert then
              for3 = For::new(l, FuncCall::new("max", -i, lowfil), FuncCall::new("min", upfil, n-1-i) ) {
                (k === i+l).print
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.print
            else
              for3 = For::new(l, lowfil, upfil) {
                if not free then
                  (k === FuncCall::new( "modulo", i+l, n)).print
                else
                  (k === i+l).print
                end
                (tt[0] === tt[0] + x[k,j] * fil[l] ).print
              }.unroll
            end
            (y[j,i] === tt[0]).print
          }.print
        }.print
      end
    }.print
  end
 
def ConvolutionGenerator::AnaRotPer(filt, center, unroll, invert, free=false )

  function_name = "ana_rot_per"

  if unroll>0 then
    function_name += "_u#{unroll}"
  end

  n = Variable::new( "n", Int, {:direction => :in} )
  ndat = Variable::new( "ndat", Int, {:direction => :in} )
  x = Variable::new( "x", Real, {:direction => :in, :dimension => [ Dimension::new( 0, n*2+1 ), Dimension::new( ndat ) ]} )
  y = Variable::new( "y", Real, {:direction => :out, :dimension => [ Dimension::new( ndat ), Dimension::new( 0, n*2+1 ) ]} )

  i = Variable::new( "i", Int )
  j = Variable::new( "j", Int )
  k = Variable::new( "k", Int )
  l = Variable::new( "l", Int )

  ci = Variable::new( "ci", Real )
  di = Variable::new( "di", Real )

  ch = Variable::new( "ch", Real, {:dimension => [ Dimension::new( -7,8 )]} )
  cg = Variable::new( "cg", Real, {:dimension => [ Dimension::new( -7,8 )]} )

  chdata = [ "-0.0033824159510050025955_wp", 
       "-0.00054213233180001068935_wp", 
      "0.031695087811525991431_wp", 
       "0.0076074873249766081919_wp",
         "-0.14329423835127266284_wp", 
       "-0.061273359067811077843_wp", 
        "0.48135965125905339159_wp", 
       "0.77718575169962802862_wp",
        "0.36444189483617893676_wp",
       "-0.051945838107881800736_wp",
       "-0.027219029917103486322_wp",
       "0.049137179673730286787_wp",
        "0.0038087520138944894631_wp",
       "-0.014952258337062199118_wp",
      "-0.00030292051472413308126_wp",
       "0.0018899503327676891843_wp" ]

  charray = ConstArray::new( chdata )

  cgdata = [ "-0.0018899503327676891843_wp",
      "-0.00030292051472413308126_wp",
      "0.014952258337062199118_wp", 
      "0.0038087520138944894631_wp", 
      "-0.049137179673730286787_wp", 
      "-0.027219029917103486322_wp", 
      "0.051945838107881800736_wp", 
      "0.36444189483617893676_wp", 
      "-0.77718575169962802862_wp",
      "0.48135965125905339159_wp", 
      "0.061273359067811077843_wp",
      "-0.14329423835127266284_wp", 
      "-0.0076074873249766081919_wp",
      "0.031695087811525991431_wp", 
      "0.00054213233180001068935_wp",
      "-0.0033824159510050025955_wp", ]
  cgarray = ConstArray::new( cgdata )

  ch = Variable::new( "ch", Real, {:constant => charray, :dimension => [ Dimension::new( -7, 8 ) ]} )

  cg = Variable::new( "cg", Real, {:constant => cgarray, :dimension => [ Dimension::new( -7, 8 ) ]} )

p = Procedure::new( function_name, [n, ndat, x, y] ) {

  i.decl
  j.decl
  k.decl
  l.decl

  ci.decl
  di.decl

  ch.decl
  cg.decl

  for1 = For::new( j, 1, ndat, 1 )
  for1.print
  for2 = For::new( i, 0, n, 1 ) {
    (ci === "0.e0_wp").print
    (di === "0.e0_wp").print

    for3 = For::new( l, -7, 8, 1 ) {
      (k === FuncCall::new( "modulo", l + i * 2, n * 2 + 2 )).print
      (ci === ci + ch[l] * x[k, j]).print
      (di === di + cg[l] * x[k, j]).print
    }.print

    (y[j, i] === ci).print
    (y[j,n + 1 + i ] === di).print 
  }.print
  for1.close
}.print

  end



def ConvolutionGenerator::SynRotPer(filt, center, unroll, invert, free=false )

    function_name = "syn_rot_per"

  if unroll>0 then
    function_name += "_u#{unroll}"
  end

  n = Variable::new( "n", Int, {:direction => :in} )
  ndat = Variable::new( "ndat", Int, {:direction => :in} )
  x = Variable::new( "x", Real, {:direction => :in, :dimension => [ Dimension::new( 0, n*2+1 ), Dimension::new( ndat ) ]} )
  y = Variable::new( "y", Real, {:direction => :out, :dimension => [ Dimension::new( ndat ), Dimension::new( 0, n*2+1 ) ]} )

  i = Variable::new( "i", Int )
  j = Variable::new( "j", Int )
  k = Variable::new( "k", Int )
  l = Variable::new( "l", Int )

  so = Variable::new( "so", Real )
  se = Variable::new( "se", Real )

  ch = Variable::new( "ch", Real, {:dimension => [ Dimension::new( -7,8 )]} )
  cg = Variable::new( "cg", Real, {:dimension => [ Dimension::new( -7,8 )]} )

  chdata = [ "-0.0033824159510050025955_wp", 
       "-0.00054213233180001068935_wp", 
      "0.031695087811525991431_wp", 
       "0.0076074873249766081919_wp",
         "-0.14329423835127266284_wp", 
       "-0.061273359067811077843_wp", 
        "0.48135965125905339159_wp", 
       "0.77718575169962802862_wp",
        "0.36444189483617893676_wp",
       "-0.051945838107881800736_wp",
       "-0.027219029917103486322_wp",
       "0.049137179673730286787_wp",
        "0.0038087520138944894631_wp",
       "-0.014952258337062199118_wp",
      "-0.00030292051472413308126_wp",
       "0.0018899503327676891843_wp" ]

  charray = ConstArray::new( chdata )

  cgdata = [ "-0.0018899503327676891843_wp",
      "-0.00030292051472413308126_wp",
      "0.014952258337062199118_wp", 
      "0.0038087520138944894631_wp", 
      "-0.049137179673730286787_wp", 
      "-0.027219029917103486322_wp", 
      "0.051945838107881800736_wp", 
      "0.36444189483617893676_wp", 
      "-0.77718575169962802862_wp",
      "0.48135965125905339159_wp", 
      "0.061273359067811077843_wp",
      "-0.14329423835127266284_wp", 
      "-0.0076074873249766081919_wp",
      "0.031695087811525991431_wp", 
      "0.00054213233180001068935_wp",
      "-0.0033824159510050025955_wp", ]
  cgarray = ConstArray::new( cgdata )

  ch = Variable::new( "ch", Real, {:constant => charray, :dimension => [ Dimension::new( -8, 9 ) ]} )

  cg = Variable::new( "cg", Real, {:constant => cgarray, :dimension => [ Dimension::new( -8, 9 ) ]} )

p = Procedure::new( function_name, [n, ndat, x, y] ) {

  i.decl
  j.decl
  k.decl
  l.decl

  so.decl
  se.decl

  ch.decl
  cg.decl

  for1 = For::new( j, 1, ndat, 1 )
  for1.print
  for2 = For::new( i, 0, n, 1 ) {
    (se === "0.e0_wp").print
    (so === "0.e0_wp").print

    for3 = For::new( l, -4, 4, 1 ) {
      (k === FuncCall::new( "modulo", i - l, n + 1)).print
      (se === se + ch[l * 2] * x[k, j] + cg[l * 2] * x[n + 1 + k, j] ).print
      (so === so + cg[l * 2 + 1] * x[k, j] + cg[l * 2 + 1] * x[n + 1 + k, j] ).print
    }.print

    (y[j, i * 2] === se).print
    (y[j, i * 2 + 1] === so).print 
  }.print
  for1.close
}.print

end

def ConvolutionGenerator::ConvrotNPer(filt, center, unroll, invert, free=false )
  
    function_name = "convrot_n_per"

  if unroll>0 then
    function_name += "_u#{unroll}"
  end

  n1 = Variable::new( "n1", Int, {:direction => :in} )
  ndat = Variable::new( "ndat", Int, {:direction => :in} )
  x = Variable::new( "x", Real, {:direction => :in, :dimension => [ Dimension::new( 0, n1 ), Dimension::new( ndat ) ]} )
  y = Variable::new( "y", Real, {:direction => :out, :dimension => [ Dimension::new( ndat ), Dimension::new( 0, n1 ) ]} )

  lowfil = Variable::new( "lowfil", Int, {:constant => -8} )
  lupfil = Variable::new( "lupfil", Int, {:constant => 7} )

  i = Variable::new( "i", Int )
  j = Variable::new( "j", Int )
  k = Variable::new( "k", Int )
  l = Variable::new( "l", Int )

  tt = Variable::new( "tt", Real )

  fildata = [ "8.4334247333529341094733325815816e-7_wp",
    "-0.1290557201342060969516786758559028e-4_wp",
    "0.8762984476210559564689161894116397e-4_wp",
    "-0.30158038132690463167163703826169879e-3_wp",
    "0.174723713672993903449447812749852942e-2_wp",
    "-0.942047030201080385922711540948195075e-2_wp",
    "0.2373821463724942397566389712597274535e-1_wp",
    "0.612625895831207982195380597e-1_wp",
    "0.9940415697834003993178616713_wp",
    "-0.604895289196983516002834636e-1_wp ",
    "-0.2103025160930381434955489412839065067e-1_wp",
    "0.1337263414854794752733423467013220997e-1_wp",
    "-0.344128144493493857280881509686821861e-2_wp",
    "0.49443227688689919192282259476750972e-3_wp",
    "-0.5185986881173432922848639136911487e-4_wp",
    "2.72734492911979659657715313017228e-6_wp" ]

  filarray = ConstArray::new( fildata )

  fil = Variable::new( "fil", Real, {:size => "wp", :dimension => [ Dimension::new( "lowfil", "lupfil" ) ], :constant => filarray } )

p = Procedure::new( function_name, [n1, ndat, x, y] ) {

  lowfil.decl
  lupfil.decl

  i.decl
  j.decl
  k.decl
  l.decl
  
  tt.decl
  fil.decl

  for1 = For::new( j, 1, ndat, 1 ) {
    for2 = For::new( i, 0, n1, 1 ) {
      (tt === "0.e0_wp").print
      for3 = For::new( l, "lowfil", "lupfil", 1 ) {
        (k === FuncCall::new( "modulo", i + l, n1 + 1)).print
        (tt === tt + x[k, j] * fil[l]).print
      }.print
    (y[j, i] === tt).print
    }.print
  }.print
}.print

end

end

FILTER = [ "8.4334247333529341094733325815816e-7",
       "-0.1290557201342060969516786758559028e-4",
       "0.8762984476210559564689161894116397e-4",
       "-0.30158038132690463167163703826169879e-3",
       "0.174723713672993903449447812749852942e-2",
       "-0.942047030201080385922711540948195075e-2",
       "0.2373821463724942397566389712597274535e-1",
       "0.612625895831207982195380597e-1",
       "0.9940415697834003993178616713",
       "-0.604895289196983516002834636e-1",
       "-0.2103025160930381434955489412839065067e-1",
       "0.1337263414854794752733423467013220997e-1",
       "-0.344128144493493857280881509686821861e-2",
       "0.49443227688689919192282259476750972e-3",
       "-0.5185986881173432922848639136911487e-4",
       "2.72734492911979659657715313017228e-6" ]


#ConvolutionGenerator::MagicFilter(FILTER,8,0,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,5,true)
#ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true)
#ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true)
ConvolutionGenerator::set_lang( ConvolutionGenerator::FORTRAN )
#ConvolutionGenerator::set_lang( ConvolutionGenerator::C )
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,5,true)
#ConvolutionGenerator::MagicFilter(FILTER,8,3,false,true)
#ConvolutionGenerator::MagicFilter(FILTER,8,4,true,true)


#ConvolutionGenerator::MagicFilter(FILTER,8,0,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,5,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,0,false)

#ConvolutionGenerator::AnaRotPer(FILTER, 8, 0, false)
#ConvolutionGenerator::SynRotPer(FILTER, 8, 0, false)
ConvolutionGenerator::ConvrotNPer(FILTER, 8, 0, false)
