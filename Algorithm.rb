
module ConvolutionGenerator
  FORTRAN = 1
  C = 2
  CL = 3
  CUDA = 4
  $output = STDOUT
  $lang = FORTRAN
  $replace_constants = true
  $default_int_size = 4
  $default_int_signed = true
  $default_real_size = 8
  $indent_level = 0
  $indent_increment = 2
  $array_start = 1

  def ConvolutionGenerator::set_default_real_size(size)
    $default_real_size = size
  end

  def ConvolutionGenerator::get_default_real_size
    return $default_real_size
  end

  def ConvolutionGenerator::set_lang(lang)
    $lang = lang
  end

  def ConvolutionGenerator::get_lang
    return $lang
  end

  def ConvolutionGenerator::set_output(output)
    $output = output
  end

  def ConvolutionGenerator::get_output
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
 
    def dereference
      return Expression::new("*",nil,self)
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
      s += " " if @operator.to_s != "++"
      s += @operator.to_s 
      s += " "
      if @operand2 then
        s += "(" if (@operator == "*" or @operator == "/" or @operator == "-" or @operator == "+") 
        s += @operand2.to_s
        s += ")" if (@operator == "*" or @operator == "/" or @operator == "-" or @operator == "+") 
      end
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang==C or $lang==CL
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
      return self.to_str_c if $lang == C or $lang == CL
    end
    def to_str_fortran
      s = ""
      s += "#{@source}(#{@indexes.join(", ")})"
      return s
    end
    def to_str_c
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
      if $replace_constants then
        begin
#         puts sub
         indx = eval(sub)
         indx = indx.to_i
#         puts indx
         return "#{@source.constant[indx]}"
        rescue Exception => e
        end
      end
      s = "#{@source}["
      s += sub + "]"
      return s
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C or $lang == CL
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
      @hash = hash
    end

    def copy
      return Variable::new(@name, @type.class, @hash)
    end
  
    def to_s
      self.to_str
    end    

    def to_str
      if @constant and $replace_constants and not @dimension then
        s = @constant.to_s 
        s += "_wp" if $lang == FORTRAN and @type and @type.size == 8
        return s
      end
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
 
    def /(x)
      return Expression::new("/",self,x)
    end
 
    def -(x)
      return Expression::new("-",self,x)
    end
 
    def -@
      return Expression::new("-",nil,self)
    end
   
    def dereference
      return Expression::new("*",nil,self)
    end
   
    def inc
      return Expression::new("++",self,nil)
    end

    def [](*args)
      return Index::new(self,args)
    end
 
    def indent
       return " "*$indent_level
    end

    def finalize
       s = ""
       s += ";" if $lang == C or $lang == CL
       s+="\n"
       return s
    end

    def header(lang=C,final=true)
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and $lang == CL
      s += "__local " if @local and $lang == CL
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
      return self.decl_c(final) if $lang == C or $lang == CL      
    end

    def decl_c(final=true)
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and $lang == CL
      s += "__local " if @local and $lang == CL
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
      if @constant
        s += " = #{@constant}"
        s += "_wp" if not @dimension and @type and @type.size == 8
      end
      s += self.finalize if final
      $output.print s if final
      return s
    end

  end

  class CustomType
    attr_reader :size, :name, :vector_length
    def initialize(hash={})
      @name = hash[:type_name]
      @size = hash[:size]
      @vector_length = hash[:vector_length]
    end
    def decl
      return "#{@name}" if $lang == C or $lang == CL
    end
  end

  class Real
    attr_reader :size
    def initialize(hash={})
      if hash[:size] then
        @size = hash[:size]
      else
        @size = $default_real_size
      end
      if hash[:vector_length] and hash[:vector_length] > 1 then
        @vector_length = hash[:vector_length]
      else
        @vector_length = 1
      end
    end
    def decl
      return "real(kind=#{@size})" if $lang == FORTRAN
      if ($lang == C or $lang == CL or $lang == CUDA) and @vector_length == 1 then
        return "float" if @size == 4
        return "double" if @size == 8
      elsif $lang == CL and @vector_length > 1 then
        return "float#{@vector_length}" if @size == 4
        return "double#{@vector_length}" if @size == 8
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
    attr_reader :headers
    def initialize(name, parameters=[], constants=[], properties={}, &block)
      @name = name
      @parameters = parameters
      @constants = constants
      @block = block
      @properties = properties
      @headers = properties[:headers]
      @headers = [] if not @headers
    end

    def header(lang=C,final=true)
      s = ""
      headers.each { |h|
        s += "#include <#{h}>\n"
      }
      if $lang == CL then
        s += "__kernel " if $lang == CL
        wgs = @properties[:reqd_work_group_size]
        if wgs then
          s += "__attribute__((reqd_work_group_size(#{wgs[0]},#{wgs[1]},#{wgs[2]}))) "
        end
      end
      trailer = ""
      trailer += "_" if lang == FORTRAN
      trailer += "_wrapper" if lang == CUDA
      if @properties[:return] then
        s += "#{@properties[:return].type.decl} "
      elsif lang == CUDA
        s += "unsigned long long int "
      else
        s += "void "
      end
      s += "#{@name}#{trailer}("
      if parameters.first then
        s += parameters.first.header(lang,false)
        parameters[1..-1].each { |p|
          s += ", "
          s += p.header(lang,false)
        }
      end
      if lang == CUDA then
        s += ", " if parameters.first
        s += "size_t *block_number, size_t *block_size"
      end
      s += ")"
      s += ";\n" if final
      $output.print s if final
      return s
    end

    def call(*parameters)
      prefix = ""
      prefix += "call " if $lang==FORTRAN
      f = FuncCall::new(@name, *parameters)
      f.prefix = prefix
      return f
    end
    def decl(final=true)
      return self.decl_fortran(final) if $lang==FORTRAN
      return self.decl_c(final) if $lang==C or $lang==CL
    end
    def close(final=true)
      return self.close_fortran(final) if $lang==FORTRAN
      return self.close_c(final) if $lang==C or $lang==CL
    end
    def close_c(final=true)
      $indent_level -= $indent_increment
      s = ""
      s += "  return #{@properties[:return]};\n" if @properties[:return]
      s += "}"
      $output.puts s if final
      return s
    end
    def close_fortran(final=true)
      $indent_level -= $indent_increment
      s = ""
      if @properties[:return] then
        s += "  #{@name} = #{@properties[:return]}\n"
        s += "END FUNCTION #{@name}"
      else
        s += "END SUBROUTINE #{@name}"
      end
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
      s += self.header($lang,false)
      s += ";\n"
      if $lang == CL then
        s += "__kernel " if $lang == CL
        wgs = @properties[:reqd_work_group_size]
        if wgs then
          s += "__attribute__((reqd_work_group_size(#{wgs[0]},#{wgs[1]},#{wgs[2]}))) "
        end
      end
      if @properties[:return] then
        s += "#{@properties[:return].type.decl} "
      else
        s += "void "
      end
      s += "#{@name}("
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
      if @properties[:return] then
        s += "#{@properties[:return].type.decl} FUNCTION "
      else
        s += "SUBROUTINE "
      end
      s += "#{@name}("
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
        elsif $lang == C or $lang == CL then
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
        return "size_t" if $lang == C or $lang == CL
      else
        return "ptrdiff_t" if $lang == C or $lang == CL
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
      if $lang == CL then
        #char="cl_"
        char=""
        char += "u" if not @signed
        return char += "char" if @size==1
        return char += "short" if @size==2
        return char += "int" if @size==4
        return char += "long" if @size==8
      end
      if $lang == CUDA then
        char = ""
        char += "unsigned " if not @signed
        return char += "char" if @size==1
        return char += "short" if @size==2
        return char += "int" if @size==4
        return char += "long long" if @size==8
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
      return self.to_str_c if $lang == C or $lang == CL
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
      return self.to_str_c if $lang == C or $lang == CL
    end
    def to_str_c
      s = ""
      s += "(#{@operand1} ? #{@operand2} : #{@operand3})"
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C or $lang == CL
      $output.puts s if final
      return s
    end

  end
 
  class FuncCall
    attr_reader :func_name
    attr_reader :args
    attr_accessor :prefix

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
      return self.to_str_c if $lang == C or $lang == CL
    end
    def to_str_fortran
      s = ""
      s += @prefix if @prefix
      s += "#{func_name}(#{@args.join(", ")})"
    end
    def to_str_c
      s = ""
      s += @prefix if @prefix
      s += "#{func_name}(#{@args.join(", ")})"
    end
    def print(final=true)
      s=""
      s += " "*$indent_level if final
      s += self.to_str
      s += ";" if final and $lang == C or $lang == CL
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
      return self.to_str_c if $lang == C or $lang == CL
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
      begin
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
      rescue Exception => e
        return self.print(*args) if not ( start and e and step )
      end
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
      return self.close_c(final) if $lang == C or $lang == CL
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

