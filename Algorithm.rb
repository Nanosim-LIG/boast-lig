
module ConvolutionGenerator
  FORTRAN = 1
  C = 2
  OpenCL = 3
  $output = STDOUT
  $lang = FORTRAN
  $replace_constants = true
  $default_int_size = 4
  $default_int_signed = true
  $default_real_size = 8
  $indent_level = 0
  $indent_increment = 2
  $array_start = 1

  def ConvolutionGenerator::set_lang(lang)
    $lang = lang
  end

  def ConvolutionGenerator::get_lang
    return $lang
  end

  def ConvolutionGenerator::set_output(output)
    $output = output
  end

  def ConvolutionGenerator::get_output(output)
    return $output
  end

  class Expression
    attr_reader :operator
    attr_reader :operand1
    attr_reader :operand2
    attr_reader :operand3
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

    def ==(x)
      return Expression::new("==",self,x)
    end
 
    def +(x)
      return Expression::new("+",self,x)
    end
 
    def <(x)
      return Expression::new("<",self,x)
    end
 
    def >=(x)
      return Expression::new(">=",self,x)
    end
 
    def *(x)
      return Expression::new("*",self,x)
    end

    def /(x)
      return Expression::new("/",self,x)
    end
 
    def -(x)
      return Expression::new("-",self,x)
    end

    def to_str
      s = ""
      if @operand1 then
        s += "(" if (@operator == "*" or @operator == "/") 
        s += @operand1.to_s
        s += ")" if (@operator == "*" or @operator == "/") 
      end        
      s += " " + @operator.to_s + " "
      s += "(" if (@operator == "*" or @operator == "/" or @operator == "-" or @operator == "+") 
#      s += "(" if (@operator == "*" or @operator == "/" or @operator == "-")
      s += @operand2.to_s
#      s += ")" if (@operator == "*" or @operator == "/" or @operator == "-") 
      s += ")" if (@operator == "*" or @operator == "/" or @operator == "-" or @operator == "+") 
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang==C or $lang==OpenCL
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
      return self.to_str_c if $lang == C or $lang == OpenCL
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
      sub = "#{@indexes.first} - #{start}"
      i=1
      ss = ""
      @source.dimension[0..-2].each{ |d|
        if d.val2 then
          ss += " * (#{d.val2} - (#{d.val1}) + 1)"
        else
          ss += " * #{d.val1}"
        end
        dim = @source.dimension[i]
        if dim.val2 then
          start = dim.val1
        else
          start = 1
        end
        sub += " + (#{@indexes[i]} - (#{start}))"+ss
        i+=1
      }
      s += sub + "]"
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C or $lang = OpenCL
      $output.puts s if final
      return s
    end
  end

  class Variable
    attr_reader :name
    attr_accessor :direction
    attr_accessor :constant
    attr_reader :type
    attr_reader :dimension
    attr_reader :local
    def initialize(name,type,hash={})
      @name = name
      @direction = hash[:direction]
      @constant = hash[:constant]
      @dimension = hash[:dimension]
      @local = hash[:local]
      @type = type::new(hash)
    end

    def copy
      return Variable::new(@name, @type.class, {:direction => @direction, :constant => @constant, :dimension => @dimension, :local => @local})
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
 
    def ==(x)
      return Expression::new("==",self,x)
    end

    def <(x)
      return Expression::new("<",self,x)
    end
 
    def >=(x)
      return Expression::new(">=",self,x)
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
       s += ";" if $lang == C or $lang == OpenCL
       s+="\n"
       return s
    end

    def header(lang=C,final=true)
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and $lang == OpenCL
      s += "__local " if @local and $lang == OpenCL
      s += @type.decl
      if(@dimension and not @constant and not @local) then
        s += " *"
      end
      if not @dimension and lang == FORTRAN then
        s += " *"
      end
      s += " #{@name}"
      if(@dimension and @constant) then
        s += "[]"
      end
      if(@dimension and @local) then
         s +="["
         s += @dimension.reverse.join("*")
         s +="]"
      end 
      s += " = #{@constant}" if @constant
      s += self.finalize if final
      $output.print s if final
      return s
    end

    def decl(final=true)
      return self.decl_fortran(final) if $lang == FORTRAN
      return self.decl_c(final) if $lang == C or $lang == OpenCL      
    end

    def decl_c(final=true)
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and $lang == OpenCL
      s += "__local " if @local and $lang == OpenCL
      s += @type.decl
      if(@dimension and not @constant and not @local) then
        s += " *"
      end
      s += " #{@name}"
      if(@dimension and @constant) then
        s += "[]"
      end
      if(@dimension and @local) then
         s +="["
         s += @dimension.reverse.join("*")
         s +="]"
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
      if $lang == OpenCL then
        return "cl_float" if @size == 4
        return "cl_double" if @size == 8
      end
    end
  end

  class CodeBlock
     def initialize(&block)
       @block = block
     end

     def print(final=true)
      s=""
      s += " "*$indent_level if final
      $indent_level += $indent_increment
      $output.puts s if final
      if @block then
        s += "\n"
        @block.call
      end
      return s
    end 
  end

  class Procedure
    attr_reader :name
    attr_reader :parameters
    attr_reader :constants
    attr_reader :properties
    def initialize(name, parameters=[], constants=[], properties={}, &block)
      @name = name
      @parameters = parameters
      @constants = constants
      @block = block
      @properties = properties
    end
    def header(lang=C,final=true)
      s = ""
      if $lang == OpenCL then
        s += "__kernel " if $lang == OpenCL
        wgs = @properties[:reqd_work_group_size]
        if wgs then
          s += "__attribute__((reqd_work_group_size(#{wgs[0]},#{wgs[1]},#{wgs[2]}))) "
        end
      end
      trailer = ""
      trailer += "_" if lang == FORTRAN
      s += "void #{@name}#{trailer}("
      if parameters.first then
        s += parameters.first.header(lang,false)
        parameters[1..-1].each { |p|
          s += ", "
          s += p.header(lang,false)
        }
      end
      s += ")"
      s += ";\n" if final
      $output.print s if final
      return s
    end
    def decl(final=true)
      return self.decl_fortran(final) if $lang==FORTRAN
      return self.decl_c(final) if $lang==C or $lang==OpenCL
    end
    def close(final=true)
      return self.close_fortran(final) if $lang==FORTRAN
      return self.close_c(final) if $lang==C or $lang==OpenCL
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
      if $lang == OpenCL then
        s += "__kernel " if $lang == OpenCL
        wgs = @properties[:reqd_work_group_size]
        if wgs then
          s += "__attribute__((reqd_work_group_size(#{wgs[0]},#{wgs[1]},#{wgs[2]}))) "
        end
      end
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
      s += " "*$indent_level + "integer, parameter :: wp=kind(1.0d0)\n"
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
      if not val2 then
        @val1 = $array_start
        @val2 = val1 + $array_start - 1
      else
        @val1 = val1
        @val2 = val2
      end
    end
    def to_str
      s = ""
      if val2 then
        if $lang == FORTRAN then
          s += val1.to_s
          s += ":"
          s += val2.to_s
        elsif $lang == C or $lang == OpenCL then
          s += (val2 - val1 + 1).to_s
        end
      else
        s += val1.to_s
      end
    end
    def to_s
      self.to_str
    end
  end

  class Sizet
    attr_reader :signed
    def initialize(hash={})
      if hash[:signed] != nil then
        @signed = hash[:signed]
      end
    end
    def decl
      return "integer(kind=#{$default_int_signed})" if $lang == FORTRAN
      if not @signed then
        return "size_t" if $lang == C or $lang == OpenCL
      else
        return "ptrdiff_t" if $lang == C or $lang == OpenCL
      end
    end
  end
 
  class Int
    attr_reader :size
    attr_reader :signed
    def initialize(hash={})
      if hash[:size] then
        @size = hash[:size]
      else
        @size = $default_int_size
      end
      if hash[:signed] != nil then
        @signed = hash[:signed]
      else
        @signed = $default_int_signed
      end
    end
    def decl
      return "integer(kind=#{@size})" if $lang == FORTRAN
      return "int#{8*@size}_t" if $lang == C
      if $lang == OpenCL then
        char="cl_"
        char += "u" if not @signed
        return char += "char" if @size==1
        return char += "short" if @size==2
        return char += "int" if @size==4
        return char += "long" if @size==8
      end
    end
  end

  class ConstArray < Array
    def initialize(array,type = nil)
      super(array)
      @type = type::new if type
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C or $lang == OpenCL
    end
    def to_str_fortran
      s = ""
      return s if not self.first
      s += "(/ &\n"
      s += self.first
      s += "_wp" if @type and @type.size == 8
      self[1..-1].each { |v|
        s += ", &\n"+v
        s += "_wp" if @type and @type.size == 8
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

  class Ternary
    attr_reader :operand1
    attr_reader :operand2
    attr_reader :operand3
    
    def initialize(x,y,z)
      @operand1 = x
      @operand2 = y
      @operand3 = z
    end

    def ==(x)
      return Expression::new("==",self,x)
    end

    def +(x)
      return Expression::new("+",self,x)
    end
 
    def <(x)
      return Expression::new("<",self,x)
    end
 
    def >=(x)
      return Expression::new(">=",self,x)
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
 
    def to_s
      self.to_str
    end

    def to_str
      raise "Ternary operator unsupported in Fortran!" if $lang == FORTRAN
      return self.to_str_c if $lang == C or $lang == OpenCL
    end
    def to_str_c
      s = ""
      s += "(#{@operand1} ? #{@operand2} : #{@operand3})"
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C or $lang = OpenCL
      $output.puts s if final
      return s
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

    def ==(x)
      return Expression::new("==",self,x)
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
 

    def to_str
      return self.to_str_fortran if $lang == FORTRAN
      return self.to_str_c if $lang == C or $lang == OpenCL
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
      s += ";" if final and $lang == C or $lang = OpenCL
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
      return self.to_str_c if $lang == C or $lang == OpenCL
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

    def unroll(*args)
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
        @block.call(*args)
      }
      @iterator.constant = nil
    end

    def print(*args)
      final = true
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      $indent_level += $indent_increment
      $output.puts s if final
      if @block then
        s += "\n"
        @block.call(*args)
        s += self.close
      end
      return s
    end

    def close(final=true)
      return self.close_fortran(final) if $lang == FORTRAN
      return self.close_c(final) if $lang == C or $lang == OpenCL
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

end


#ConvolutionGenerator::MagicFilter(FILTER,8,0,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,5,false)
#ConvolutionGenerator::MagicFilter(FILTER,8,2,false, true)

#ConvolutionGenerator::AnaRotPer(FILTER, 8, 0, false)
