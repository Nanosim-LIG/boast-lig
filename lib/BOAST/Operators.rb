module BOAST

  class Operator
    def Operator.get_vector_name(type)
      case BOAST::get_architecture
      when X86
        case type
        when Int
          size = "#{type.size*8}"
          name = ""
          if type.total_size*8 > 64
            name += "e"
          end
          if type.vector_length > 1 then
            name += "p"
          else
            name = "s"
          end
          if type.signed then
            name += "i"
          else
            name += "u"
          end
          return name += size
        when Real
          case type.size
          when 4
            return "ps" if type.vector_length > 1
            return "ss"
          when 8
            return "pd" if type.vector_length > 1
            return "sd"
          end
        else
          raise "Undefined vector type!"
        end
      else
        raise "Unsupported architecture!"
      end
    end

    def Operator.convert(arg, type)
      case BOAST::get_architecture
      when X86
        s1 = arg.type.total_size*8
        s2 = type.total_size*8
        n1 = get_vector_name(arg.type)
        n2 = get_vector_name(type)
        if s1 <= 128 and s2 <= 128 then
          return "_mm_cvt#{n1}_#{n2}( #{arg} )"
        elsif [s1, s2].max <= 256 then
          return "_mm256_cvt#{n1}_#{n2}( #{arg} )"
        elsif [s1, s2].max <= 512 then
          return "_mm512_cvt#{n1}_#{n2}( #{arg} )"
        end
      else
        raise "Unsupported architecture!"
      end
    end
  end

  class BasicBinaryOperator < Operator

    def BasicBinaryOperator.to_s(arg1, arg2, return_type)
      #puts "#{arg1.class} * #{arg2.class} : #{arg1} * #{arg2}"
      if BOAST::get_lang == C and (arg1.class == Variable and arg2.class == Variable) and (arg1.type.vector_length > 1 or arg2.type.vector_length > 1) then
        raise "Vectors have different length: #{arg1} #{arg1.type.vector_length}, #{arg2} #{arg2.type.vector_length}" if arg1.type.vector_length != arg2.type.vector_length
        #puts "#{arg1.type.signed} #{arg2.type.signed} #{return_type.type.signed}"
        if arg1.type != return_type.type
          a1 = convert(arg1, return_type.type)
        else
          a1 = "#{arg1}"
        end
        if arg2.type != return_type.type
          a2 = convert(arg2, return_type.type)
        else
          a2 = "#{arg2}"
        end
	  return_name = get_vector_name(return_type.type)
        case BOAST::get_architecture
        when X86
          intr_name = "_mm"
          size = return_type.type.total_size * 8
          if size > 128 then
            intr_name += "#{size}"
          end
          intr_name += "_#{intr_name_X86}_#{return_name}"
          return "#{intr_name}( #{a1}, #{a2} )"
        else
          raise "Unsupported architecture!"
        end
      else
        return basic_usage( arg1, arg2 )
      end
    end
  end

  class Affectation < Operator
    def Affectation.to_s(arg1, arg2, return_type)
      if BOAST::get_lang == C then
        if arg1.class == Variable and arg1.type.vector_length > 1 then
          #puts "#{arg1.type.vector_length} #{arg2.type.vector_length}"
          if arg1.type == arg2.type then
            return basic_usage(arg1, arg2)
          elsif arg1.type.vector_length == arg2.type.vector_length then
            return "#{arg1} = #{convert(arg2, arg1.type)}"
          elsif arg2.type.vector_length == 1 then
            case BOAST::get_architecture
            when X86
              intr_name = "_mm"
              size = arg1.type.total_size*8
              if size > 128 then
                intr_name += "#{size}"
              end
              intr_name += "_load_"
              if arg2.type.class == Int then
                intr_name += "si#{size}"
              else
                intr_name += "#{get_vector_name(arg1.type)}"
              end
              a2 = "#{arg2}"
              if a2[0] != "*" then
                a2 = "&" + a2
              else
                a2 = a2[1..-1]
              end
              return "#{arg1} = #{intr_name}( (#{arg1.type.decl} * ) #{a2} )"
            else
              raise "Unsupported architecture!"
            end
          else
            raise "Unknown convertion between vectors of different length!"
          end
        elsif arg2.class == Variable and arg2.type.vector_length > 1 then
          case BOAST::get_architecture
          when X86
            intr_name = "_mm"
            size = arg2.type.total_size*8
            if size > 128 then
              intr_name += "#{size}"
            end
            intr_name += "_store_"
            if arg2.type.class == Int then
              intr_name += "si#{size}"
            else
              intr_name += "#{get_vector_name(arg2.type)}"
            end
            a1 = "#{arg1}"
            if a1[0] != "*" then
              a1 = "&" + a1
            else
              a1 = a1[1..-1]
            end
            return "#{intr_name}((#{arg2.type.decl} * ) #{a1}, #{arg2} )"
          else
            raise "Unsupported architecture!"
          end
        else
          return basic_usage(arg1, arg2)
        end
      else
        return basic_usage(arg1, arg2)
      end
    end

    def Affectation.basic_usage(arg1, arg2)
      return "#{arg1} = #{arg2}"
    end
  end

  class Multiplication < BasicBinaryOperator
    class << self

      def symbol
        return "*"
      end

      def intr_name_X86
        return "mul"
      end

      def basic_usage(arg1, arg2)
        return "(#{arg1}) * (#{arg2})" 
      end
  
    end
  end

  class Addition < BasicBinaryOperator
    class << self

      def symbol
        return "+"
      end

      def intr_name_X86
        return "add"
      end
  
      def basic_usage(arg1, arg2)
        return "#{arg1} + #{arg2}" 
      end
  
    end
  end

  class Substraction < BasicBinaryOperator
    class << self

      def symbol
        return "-"
      end

      def intr_name_X86
        return "sub"
      end
  
      def basic_usage(arg1, arg2)
        return "#{arg1} - (#{arg2})" 
      end
  
    end
  end

  class Division < BasicBinaryOperator
    class << self

      def symbol
        return "/"
      end

      def intr_name_X86
        return "div"
      end
  
      def basic_usage(arg1, arg2)
        return "(#{arg1}) / (#{arg2})" 
      end
  
    end
  end

  class Minus < Operator
    def Minus.to_s(arg1, arg2, return_type)
      return " -(#{arg2})"
    end
  end

  class Not < Operator
    def Not.to_s(arg1, arg2, return_type)
      return " ! #{arg2}"
    end
  end

end
