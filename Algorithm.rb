class Object
  alias_method :orig_method_missing, :method_missing
  
  def method_missing(m, *a, &b)
    s=nil
    klass = begin
      s = (self.is_a?(Module) ? self : self.class)
      s.const_get(m)
    rescue NameError
    end
    
    return klass.send(:parens, *a, &b)  if klass.respond_to? :parens

    return BOAST::FuncCall::new(m,*a,&b) if s == BOAST

    orig_method_missing m, *a, &b
  end
end

module BOAST

  FORTRAN = 1
  C = 2
  CL = 3
  CUDA = 4
  @@output = STDOUT
  @@lang = FORTRAN
  @@replace_constants = true
  @@default_int_size = 4
  @@default_int_signed = true
  @@default_real_size = 8
  @@indent_level = 0
  @@indent_increment = 2
  @@array_start = 1
  @@chain_code = false

  @@env = Hash.new{|h, k| h[k] = []}


  def BOAST::push_env(vars = {})
    vars.each { |key,value|
      var = nil
      begin
        var = BOAST::class_variable_get("@@"+key.to_s)
      rescue
        raise "Unknown module variable #{key}!"
      end
      @@env[key].push(var)
      BOAST::class_variable_set("@@"+key.to_s, value)
    }
  end

  def BOAST::pop_env(*vars)
    vars.each { |key|
      raise "Unknown module variable #{key}!" unless @@env.has_key?(key)
      ret = @@env[key].pop
      raise "No stored value for #{key}!" if ret.nil?
      BOAST::class_variable_set("@@"+key.to_s, ret)
    }
  end

  def BOAST::print(a)
    a.print
  end

  def BOAST::decl(*a)
    a.each { |d|
      d.decl
    }
  end

  def BOAST::close(a)
    a.close
  end

  def BOAST::set_indent_level(level)
    @@indent_level = level
  end

  def BOAST::get_indent_level
    return @@indent_level
  end

  def BOAST::get_indent_increment
    return @@indent_increment
  end

  def BOAST::increment_indent_level(increment = @@indent_increment)
    @@indent_level += increment
  end
  
  def BOAST::decrement_indent_level(increment = @@indent_increment)
    @@indent_level -= increment
  end
  
  def BOAST::set_replace_constants(replace_constants)
    @@replace_constants = replace_constants
  end

  def BOAST::get_replace_constants
    return @@replace_constants
  end

  def BOAST::set_default_int_signed(signed)
    @@default_int_signed = signed
  end

  def BOAST::get_default_int_signed
    return @@default_int_signed
  end

  def BOAST::set_default_int_size(size)
    @@default_int_size = size
  end

  def BOAST::get_default_int_size
    return @@default_int_size
  end

  def BOAST::set_default_real_size(size)
    @@default_real_size = size
  end

  def BOAST::get_default_real_size
    return @@default_real_size
  end

  def BOAST::set_lang(lang)
    @@lang = lang
  end

  def BOAST::get_lang
    return @@lang
  end

  def BOAST::set_output(output)
    @@output = output
  end

  def BOAST::get_output
    return @@output
  end

  def BOAST::set_chain_code(chain_code)
    @@chain_code = chain_code
  end

  def BOAST::get_chain_code
    return @@chain_code
  end

  def BOAST::set_array_start(array_start)
    @@array_start = array_start
  end

  def BOAST::get_array_start
    return @@array_start
  end

  class Pragma
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :name
    attr_reader :options

    def initialize(name, options)
      @name = name
      @options = options
    end

    def to_s
      self.to_str
    end

    def to_str
      s = ""
      if BOAST::get_lang == FORTRAN then
        s += "$!"
      else
        s += "#pragma"
      end
      @options.each{ |opt|
        s += " #{opt}"
      }
      return s
    end

    def print(final = true)
      s=""
      s += self.to_str
      BOAST::get_output.puts s if final
      return s
    end
  end

  def BOAST::Return(value)
    return Expression("return",nil, value)
  end

  class Expression
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

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

    def ==(x)
      return Expression::new("==",self,x)
    end
 
    def +(x)
      return Expression::new("+",self,x)
    end
 
    def >(x)
      return Expression::new(">",self,x)
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

    def struct_reference(x)
      return Expression::new(".",self,x)
    end
 
    def -(x)
      return Expression::new("-",self,x)
    end

    def -@
      return Expression::new("-",nil,self)
    end
   
    def to_str
      s = ""
      if @operand1 then
        s += "(" if (@operator == "*" or @operator == "/") 
        s += @operand1.to_s
        s += ")" if (@operator == "*" or @operator == "/") 
      end        
      s += " " unless @operator.to_s == "++" or @operator.to_s == "."
      s += @operator.to_s 
      s += " " unless @operator.to_s == "." or @operator.to_s == "&"
      if @operand2 then
        s += "(" if (@operator == "*" or @operator == "/" or @operator == "-") 
        s += @operand2.to_s
        s += ")" if (@operator == "*" or @operator == "/" or @operator == "-") 
      end
      return s
    end
    def print(final=true)
      s=""
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      s += ";" if final and [C, CL, CUDA].include?( BOAST::get_lang ) 
      BOAST::get_output.puts s if final
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
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_fortran
      s = ""
      s += "#{@source}(#{@indexes.join(", ")})"
      return s
    end
    def to_str_texture
      raise "Unsupported language #{BOAST::get_lang} for texture!" if not [CL, CUDA].include?( BOAST::get_lang )
      raise "Write is unsupported for textures!" if not ( @source.constant or @source.direction == :in )
      dim_number = 1
      if @source.dimension then
        dim_number == @source.dimension.size
      end
      raise "Unsupported number of dimension: #{dim_number}!" if dim_number > 3
      s = ""
      if BOAST::get_lang == CL then
        s += "as_#{@source.type.decl}("
        s += "read_imageui(#{@source}, #{@source.sampler}, "
        if dim_number == 1 then
          s += "int2(#{@indexes[0]},0)"
        else
          if dim_number == 2 then
            s += "int2("
          else
            s += "int3("
          end
          s += "#{@indexes.join(", ")})"
        end
        s += ")"
        if @source.type.size == 4 then
          s += ".x"
        elsif @source.type.size == 8 then
          s += ".xy"
        end
        s += ")"
      else
        s += "tex#{dim_number}Dfetch(#{@source},"
        if dim_number == 1 then
          s += "#{@indexes[0]}"
        else
          if dim_number == 2 then
            s += "int2("
          else
            s += "int3("
          end
          s += "#{@indexes.join(", ")})"
        end
        s += ")"
      end
      return s
    end
    def to_str_c
      return to_str_texture if @source.texture
      dim = @source.dimension.first
      if dim.val2 then
        start = dim.val1
      else
        start = BOAST::get_array_start
      end
      sub = "#{@indexes.first} - #{start}"
      i=1
      ss = ""
      @source.dimension[0..-2].each{ |d|
        if d.size then
          ss += " * (#{d.size})"
        elsif d.val2 then
          ss += " * (#{d.val2} - (#{d.val1}) + 1)"
        else
          raise "Unkwown dimension size!"
        end
        dim = @source.dimension[i]
        if dim.val2 then
          start = dim.val1
        else
          start = BOAST::get_array_start
        end
        sub += " + (#{@indexes[i]} - (#{start}))"+ss
        i+=1
      }
      if BOAST::get_replace_constants then
        begin
#         puts sub
         indx = eval(sub)
         indx = indx.to_i
#         puts indx
         return "#{@source.constant[indx]}"
        rescue Exception => e
        end
      end
      s = "#{@source}[" + sub + "]"
      return s
    end
    def print(final=true)
      s=""
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      s += ";" if final and [C, CL, CUDA].include?( BOAST::get_lang )
      BOAST::get_output.puts s if final
      return s
    end
  end

  class CStruct
    attr_reader :name, :members, :members_array
    def self.parens(*args,&block)
      return Variable::new(args[0], self, *args[1..-1], &block)
    end

    def initialize(hash={})
      @name = hash[:type_name]
      @members = {}
      @members_array = []
      hash[:members].each { |m|
        mc = m.copy
        @members_array.push(mc)
        @members[mc.name] = mc
      }
    end

    def decl
      return "struct #{@name}" if [C, CL, CUDA].include?( BOAST::get_lang )
      return "TYPE(#{@name})" if BOAST::get_lang == FORTRAN
    end

    def finalize
       s = ""
       s += ";" if [C, CL, CUDA].include?( BOAST::get_lang )
       s+="\n"
       return s
    end

    def indent
       return " "*BOAST::get_indent_level
    end

    def header
      return header_c if [C, CL, CUDA].include?( BOAST::get_lang )
      return header_fortran if BOAST::get_lang == FORTRAN
      raise "Unsupported language!"
    end

    def header_c(final = true)
      s = ""
      s += self.indent if final
      s += self.decl + " {\n"
      @members_array.each { |value|
         s+= self.indent if final
         s+= " "*BOAST::get_indent_increment + value.decl(false)+";\n"
      }
      s += self.indent if final
      s += "}"
      s += self.finalize if final
      BOAST::get_output.print s if final
      return s
    end
    
    def header_fortran(final = true)
      s = ""
      s += self.indent if final
      s += "TYPE :: #{@name}\n"
      members_array.each { |value|
         s+= self.indent if final
         s+= " "*BOAST::get_indent_increment + value.decl(false)+"\n"
      }
      s += self.indent if final
      s += "END TYPE #{@name}"
      s += self.finalize if final
      BOAST::get_output.print s if final
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
      return "#{@name}" if [C, CL, CUDA].include?( BOAST::get_lang )
    end
  end


  class Variable
    alias_method :orig_method_missing, :method_missing

    def method_missing(m, *a, &b)
      return self.struct_reference(type.members[m.to_s]) if type.members[m.to_s]
      return self.orig_method_missing(m, *a, &b)
    end

    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :name
    attr_accessor :direction
    attr_accessor :constant
    attr_reader :type
    attr_reader :dimension
    attr_reader :local
    attr_reader :texture
    attr_reader :sampler
    attr_reader :replace_constant

    def initialize(name,type,hash={})
      @name = name
      @direction = hash[:direction] ? hash[:direction] : hash[:dir]
      @constant = hash[:constant] ? hash[:constant]  : hash[:const]
      @dimension = hash[:dimension] ? hash[:dimension] : hash[:dim]
      @local = hash[:local]
      @texture = hash[:texture]
      if not hash[:replace_constant].nil? then
        @replace_constant = hash[:replace_constant]
      else
        @replace_constant = true
      end
      if @texture and BOAST::get_lang == CL then
        @sampler = Variable::new("sampler_#{name}", BOAST::CustomType,:type_name => "sampler_t" ,:replace_constant => false, :constant => "CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST")
      else
        @sampler = nil
      end
      @type = type::new(hash)
      @hash = hash
    end

    def copy(name=@name)
      return Variable::new(name, @type.class, @hash)
    end
  
    def to_s
      self.to_str
    end    

    def to_str
      if @replace_constant and @constant and BOAST::get_replace_constants and not @dimension then
        s = @constant.to_s 
        s += "_wp" if BOAST::get_lang == FORTRAN and @type and @type.size == 8
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

    def >(x)
      return Expression::new(">",self,x)
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

    def !
      return Expression::new("!",nil,self)
    end
 
    def -@
      return Expression::new("-",nil,self)
    end

    def address
      return Expression::new("&",nil,self)
    end
   
    def dereference
      return Expression::new("*",nil,self)
    end
   
    def struct_reference(x)
      return x.copy(self.name+"."+x.name) if [C, CL, CUDA].include?( BOAST::get_lang )
      return x.copy(self.name+"%"+x.name) if BOAST::get_lang == FORTRAN
    end
 
    def inc
      return Expression::new("++",self,nil)
    end

    def [](*args)
      return Index::new(self,args)
    end
 
    def indent
       return " "*BOAST::get_indent_level
    end

    def finalize
       s = ""
       s += ";" if [C, CL, CUDA].include?( BOAST::get_lang )
       s+="\n"
       return s
    end

    def decl_c(final=true)
      return decl_texture(final) if @texture
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and BOAST::get_lang == CL
      s += "__local " if @local and BOAST::get_lang == CL
      s += "__shared__ " if @local and BOAST::get_lang == CUDA
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
      BOAST::get_output.print s if final
      return s
    end

    def header(lang=C,final=true)
      return decl_texture(final) if @texture
      s = ""
      s += self.indent if final
      s += "const " if @constant or @direction == :in
      s += "__global " if @direction and @dimension and BOAST::get_lang == CL
      s += "__local " if @local and BOAST::get_lang == CL
      s += "__shared__ " if @local and BOAST::get_lang == CUDA
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
      BOAST::get_output.print s if final
      return s
    end

    def decl(final=true)
      return self.decl_fortran(final) if BOAST::get_lang == FORTRAN
      return self.decl_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end

    def decl_texture(final=true)
      raise "Unsupported language #{BOAST::get_lang} for texture!" if not [CL, CUDA].include?( BOAST::get_lang )
      raise "Write is unsupported for textures!" if not (@constant or @direction == :in)
      dim_number = 1
      if @dimension then
        dim_number == @dimension.size
      end
      raise "Unsupported number of dimension: #{dim_number}!" if dim_number > 3
      s = ""
      s += self.indent if final
      if BOAST::get_lang == CL then
        s += "__read_only "
        if dim_number < 3 then
          s += "image2d_t " #from OCL 1.2+ image1d_t is defined
        else
          s += "image3d_t "
        end
      else
        s += "texture<#{@type.decl}, cudaTextureType#{dim_number}D, cudaReadModeElementType> "
      end
      s += @name
      s += self.finalize if final
      BOAST::get_output.print s if final
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
        dim = @dimension[0].to_str
        if dim then
          s += dim
          @dimension[1..-1].each { |d|
             s += ", "
             s += d
          }
        else
          s += ":"
        end
        s += ")"
      end
      s += " :: #{@name}"
      if @constant
        s += " = #{@constant}"
        s += "_wp" if not @dimension and @type and @type.size == 8
      end
      s += self.finalize if final
      BOAST::get_output.print s if final
      return s
    end

  end



  class Real
    def self.parens(*args,&block)
      return Variable::new(args[0], self, *args[1..-1], &block)
    end

    attr_reader :size
    def initialize(hash={})
      if hash[:size] then
        @size = hash[:size]
      else
        @size = BOAST::get_default_real_size
      end
      if hash[:vector_length] and hash[:vector_length] > 1 then
        @vector_length = hash[:vector_length]
      else
        @vector_length = 1
      end
    end
    def decl
      return "real(kind=#{@size})" if BOAST::get_lang == FORTRAN
      if [C, CL, CUDA].include?( BOAST::get_lang ) and @vector_length == 1 then
        return "float" if @size == 4
        return "double" if @size == 8
      elsif [CL, CUDA].include?(BOAST::get_lang) and @vector_length > 1 then
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
      s += " "*BOAST::get_indent_level if final
      BOAST::increment_indent_level
      BOAST::get_output.puts s if final
      if @block then
        s += "\n"
        @block.call
      end
      return s
    end 
  end

  class Procedure
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

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
      if BOAST::get_lang == CL then
        s += "__kernel "
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
      BOAST::get_output.print s if final
      return s
    end

    def call(*parameters)
      prefix = ""
      prefix += "call " if BOAST::get_lang==FORTRAN
      f = FuncCall::new(@name, *parameters)
      f.prefix = prefix
      return f
    end
    def decl(final=true)
      return self.decl_fortran(final) if BOAST::get_lang==FORTRAN
      return self.decl_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang==FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      BOAST::decrement_indent_level
      s = ""
      s += "  return #{@properties[:return]};\n" if @properties[:return]
      s += "}"
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      BOAST::decrement_indent_level
      s = ""
      if @properties[:return] then
        s += "  #{@name} = #{@properties[:return]}\n"
        s += "END FUNCTION #{@name}"
      else
        s += "END SUBROUTINE #{@name}"
      end
      BOAST::get_output.puts s if final
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
#      s += self.header(BOAST::get_lang,false)
#      s += ";\n"
      if BOAST::get_lang == CL then
        s += "__kernel "
        wgs = @properties[:reqd_work_group_size]
        if wgs then
          s += "__attribute__((reqd_work_group_size(#{wgs[0]},#{wgs[1]},#{wgs[2]}))) "
        end
      elsif BOAST::get_lang == CUDA then
        if @properties[:local] then
          s += "__device__ "
        else
          s += "__global__ "
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
      BOAST::increment_indent_level
      constants.each { |c|
        s += " "*BOAST::get_indent_level
        s += c.decl(false)
        s += ";\n"
      }
      BOAST::get_output.print s if final
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
      BOAST::increment_indent_level
      s += " "*BOAST::get_indent_level + "integer, parameter :: wp=kind(1.0d0)\n"
      constants.each { |c|
        s += " "*BOAST::get_indent_level
        s += c.decl(false)
        s += "\n"
      }
      parameters.each { |p|
        s += " "*BOAST::get_indent_level
        s += p.decl(false)
        s += "\n"
      }
      BOAST::get_output.print s if final
      return s
    end
  end

  class Dimension
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :val1
    attr_reader :val2
    attr_reader :size
    def initialize(v1=nil,v2=nil)
      @size = nil
      @val1 = nil
      @val2 = nil
      if v2.nil? and v1 then
        @val1 = BOAST::get_array_start
        @val2 = v1 + BOAST::get_array_start - 1
        @size = v1
      else
        @val1 = v1
        @val2 = v2
      end
    end
    def to_str
      s = ""
      if @val2 then
        if BOAST::get_lang == FORTRAN then
          s += @val1.to_s
          s += ":"
          s += @val2.to_s
        elsif [C, CL, CUDA].include?( BOAST::get_lang ) then
          s += (@val2 - @val1 + 1).to_s
        end
      elsif @val1.nil? then
        return nil
      else
        s += @val1.to_s
      end
      return s
    end
    def to_s
      self.to_str
    end
  end

  class Sizet
    def self.parens(*args,&block)
      return Variable::new(args[0], self, *args[1..-1], &block)
    end

    attr_reader :signed
    def initialize(hash={})
      if hash[:signed] != nil then
        @signed = hash[:signed]
      end
    end
    def decl
      return "integer(kind=#{BOAST::get_default_int_signed})" if BOAST::get_lang == FORTRAN
      if not @signed then
        return "size_t" if [C, CL, CUDA].include?( BOAST::get_lang )
      else
        return "ptrdiff_t" if [C, CL, CUDA].include?( BOAST::get_lang )
      end
    end
  end
 
  class Int
    def self.parens(*args,&block)
      return Variable::new(args[0], self, *args[1..-1], &block)
    end

    attr_reader :size
    attr_reader :signed
    def initialize(hash={})
      if hash[:size] then
        @size = hash[:size]
      else
        @size = BOAST::get_default_int_size
      end
      if hash[:signed] != nil then
        @signed = hash[:signed]
      else
        @signed = BOAST::get_default_int_signed
      end
    end
    def decl
      return "integer(kind=#{@size})" if BOAST::get_lang == FORTRAN
      return "int#{8*@size}_t" if BOAST::get_lang == C
      if BOAST::get_lang == CL then
        #char="cl_"
        char=""
        char += "u" if not @signed
        return char += "char" if @size==1
        return char += "short" if @size==2
        return char += "int" if @size==4
        return char += "long" if @size==8
      elsif BOAST::get_lang == CUDA then
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
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_fortran
      s = ""
      return s if self.first.nil?
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
      return s if self.first.nil?
      s += "{\n"
      s += self.first 
      self[1..-1].each { |v|
        s += ",\n"+v
      }
      s += "}"
    end
  end

  class Ternary

    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

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
      raise "Ternary operator unsupported in FORTRAN!" if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_c
      s = ""
      s += "(#{@operand1} ? #{@operand2} : #{@operand3})"
    end
    def print(final=true)
      s=""
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      s += ";" if final and [C, CL, CUDA].include?( BOAST::get_lang )
      BOAST::get_output.puts s if final
      return s
    end

  end
 
  class FuncCall
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

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
 
    def /(x)
      return Expression::new("/",self,x)
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
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
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
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      s += ";" if final and [C, CL, CUDA].include?( BOAST::get_lang )
      BOAST::get_output.puts s if final
      return s
    end
  end

  class While
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :condition
    def initialize(condition, &block)
      @condition = condition
      @block = block
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_fortran
      s = ""
      s += "do while( #{@condition} )"
      return s
    end
    def to_str_c
      s = ""
      s += "while(#{@condition}){"
      return s
    end
    def print(*args)
      final = true
      s=""
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      BOAST::increment_indent_level      
      BOAST::get_output.puts s if final
      if @block then
        s += "\n"
        @block.call(*args)
        s += self.close
      end
      return s
    end
    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang == FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "}"
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "end do"
      BOAST::get_output.puts s if final
      return s
    end

  end

  class Else
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :condition
    def initialize(condition=nil, &block)
      @condition = condition
      @block = block
    end
    def to_s
      self.to_str
    end
    def to_str
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_fortran
      s = ""
      if @condition then
        s += "else if #{@condition} then"
      else
        s += "else"
      end
      return s
    end
    def to_str_c
      s = ""
      if @condition then
        s += "else if(#{@condition}){"
      else
        s += "else {"
      end
      return s
    end
    def print(*args)
      final = true
      s=""
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      BOAST::increment_indent_level      
      BOAST::get_output.puts s if final
      if @block then
        s += "\n"
        @block.call(*args)
        s += self.close
      end
      return s
    end
    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang == FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "}"
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "end if"
      BOAST::get_output.puts s if final
      return s
    end

  end

  class Case
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :expression
    attr_reader :constants_list

    def initialize(expression, *control)
      @expression = expression
      @constants_list = []
      @blocks = []
      if control.size < 1 then
        raise "No block given!"
      elsif control.size.even? then
        (0..control.size-1).step(2) { |i|
          @constants_list[i/2] = [control[i]].flatten
          @blocks[i/2] = control[i+1]
        }
      else
        (0..control.size-2).step(2) { |i|
          @constants_list[i/2] = [control[i]].flatten
          @blocks[i/2] = control[i+1]
        }
        @blocks.push(control.last)
      end
    end

    def to_s(*args)
      self.to_str(*args)
    end

    def to_str(constants, first= true)
      return self.to_str_fortran(constants, first) if BOAST::get_lang == FORTRAN
      return self.to_str_c(constants, first) if [C, CL, CUDA].include?( BOAST::get_lang )
    end

    def to_str_fortran(constants, first)
      s = ""
      if first then
        s += "select case #{@expression}\n"
        BOAST::increment_indent_level
      else
        BOAST::decrement_indent_level
      end
      s += " "*BOAST::get_indent_level
      if constants and constants.size>0 then
        s += "case #{constants.join(" : ")}"
      else
        s += "case default"
      end
      BOAST::increment_indent_level
      return s
    end

    def to_str_c(constants, first)
      s = ""
      if first then
      s += " "*BOAST::get_indent_level
        s += "switch(#{@expression}){\n"
        BOAST::increment_indent_level
      else
        s += " "*BOAST::get_indent_level + "break;\n"
        BOAST::decrement_indent_level
      end
      s += " "*BOAST::get_indent_level
      if constants and constants.size>0 then
        s += "case #{constants.join(" : case")} :"
      else
        s += "default :"
      end
      BOAST::increment_indent_level
      return s
    end

    def print(*args)
      first = true
      @blocks.each_index { |indx|
        s = self.to_str(@constants_list[indx],first)
        BOAST::get_output.puts s
        @blocks[indx].call(*args)
        first = false
      }
      self.close
      return self
    end
    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang == FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      s = ""
      s += " "*BOAST::get_indent_level if final
      s += "break;\n"
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "}"
      BOAST::decrement_indent_level      
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "end select"
      BOAST::decrement_indent_level      
      BOAST::get_output.puts s if final
      return s
    end

  end 
  class If
    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

    attr_reader :conditions
    def initialize(*conditions, &block)
      @conditions = []
      @blocks = []
      if conditions.size == 0 then
        raise "Illegal if construct!"
      elsif conditions.size == 1 then
        @conditions.push(conditions[0])
        @blocks.push(block)
      elsif conditions.size.even? then
        (0..conditions.size-1).step(2) { |i|
          @conditions[i/2] = conditions[i]
          @blocks[i/2] = conditions[i+1]
        }
      else
        (0..conditions.size-2).step(2) { |i|
          @conditions[i/2] = conditions[i]
          @blocks[i/2] = conditions[i+1]
        }
        @blocks.push(conditions.last)
      end
    end
    def to_s(*args)
      self.to_str(*args)
    end
    def to_str(condition, first= true)
      return self.to_str_fortran(condition, first) if BOAST::get_lang == FORTRAN
      return self.to_str_c(condition, first) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def to_str_fortran(condition, first)
      s = ""
      if first then
        s += "if #{condition} then"
      else
        if condition then
          s += "else if #{condition} then"
        else
          s += "else"
        end
      end
      return s
    end
    def to_str_c(condition, first)
      s = ""
      if first then
        s += "if(#{condition}){"
      else
        if condition then
          s += "} else if(#{condition}){"
        else
          s += "} else {"
        end
      end
      return s
    end
    def print(*args)
      s=""
      s += " "*BOAST::get_indent_level
      s += self.to_str(@conditions.first)
      BOAST::increment_indent_level      
      BOAST::get_output.puts s
      if @blocks.size > 0 then
        if @blocks[0] then
          @blocks[0].call(*args)
        end
        @blocks[1..-1].each_index { |indx|
          BOAST::decrement_indent_level      
          s=""
          s += " "*BOAST::get_indent_level 
          s += self.to_str(@conditions[1..-1][indx],false)
          BOAST::increment_indent_level
          BOAST::get_output.puts s
          @blocks[1..-1][indx].call(*args)
        }
        self.close
      end
      return self
    end
    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang == FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "}"
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "end if"
      BOAST::get_output.puts s if final
      return s
    end

  end
 
  class For
    attr_reader :iterator
    attr_reader :begin
    attr_reader :end
    attr_reader :step

    def self.parens(*args,&block)
      return self::new(*args,&block)
    end

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
      return self.to_str_fortran if BOAST::get_lang == FORTRAN
      return self.to_str_c if [C, CL, CUDA].include?( BOAST::get_lang )
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
        elsif @begin.kind_of?(Expression) then
          start = eval "#{@begin}"
        else
          start = @begin.to_i
        end
        if @end.kind_of?(Variable) then
          e = @end.constant
        elsif @end.kind_of?(Expression) then
          e = eval "#{@end}"
        else
          e = @end.to_i
        end
        if @step.kind_of?(Variable) then
          step = @step.constant
        elsif @step.kind_of?(Expression) then
          step = eval "#{@step}"
        else
          step = @step.to_i
        end
        raise "Invalid bounds (not constants)!" if not ( start and e and step )
      rescue Exception => ex
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
      s += " "*BOAST::get_indent_level if final
      s += self.to_str
      BOAST::increment_indent_level      
      BOAST::get_output.puts s if final
      if @block then
        s += "\n"
        @block.call(*args)
        s += self.close
      end
      return s
    end

    def close(final=true)
      return self.close_fortran(final) if BOAST::get_lang == FORTRAN
      return self.close_c(final) if [C, CL, CUDA].include?( BOAST::get_lang )
    end
    def close_c(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "}"
      BOAST::get_output.puts s if final
      return s
    end
    def close_fortran(final=true)
      s = ""
      BOAST::decrement_indent_level      
      s += " "*BOAST::get_indent_level if final
      s += "enddo"
      BOAST::get_output.puts s if final
      return s
    end
  end
  Var = Variable
  Dim = Dimension
  Call = FuncCall
end
ConvolutionGenerator = BOAST
