require './Algorithm.rb'
require 'stringio'
require 'rake'
require 'tempfile'
require 'rbconfig'

module ConvolutionGenerator
  class CKernel
    include Rake::DSL
    attr_accessor :code
    attr_accessor :procedure
    attr_accessor :lang
    attr_accessor :binary
    
    def initialize
      @code = StringIO::new
    end

    def print
      @code.rewind
      puts @code.read
    end
    def build
      includes = "-I#{RbConfig::CONFIG["archdir"]}"
      includes += " -I#{RbConfig::CONFIG["rubyhdrdir"]} -I#{RbConfig::CONFIG["rubyhdrdir"]}/#{RbConfig::CONFIG["arch"]}"
      ldflags = "-L#{RbConfig::CONFIG["libdir"]} #{RbConfig::CONFIG["LIBRUBYARG"]}"
      cflags = "-Wall -fPIC #{includes}"
      fcflags = "-Wall -fPIC"
      rule '.o' => '.c' do |t|
        sh "cc #{cflags} -c -o #{t.name} #{t.source}"
      end
      rule '.o' => '.f90' do |t|
        sh "gfortran #{fcflags} -c -o #{t.name} #{t.source}"
      end

      extension = ".c" if @lang == ConvolutionGenerator::C
      extension = ".f90" if @lang == ConvolutionGenerator::FORTRAN
      extension = ".cl" if @lang == ConvolutionGenerator::OpenCL
      source_file = Tempfile::new([@procedure.name,extension])
      path = source_file.path
      target = path.chomp(File::extname(path))+".o"
      @code.rewind
      source_file.puts "#include <inttypes.h>" if @lang == ConvolutionGenerator::C
      source_file.write @code.read
      source_file.close

      previous_lang = $lang
      previous_output = $output
      ConvolutionGenerator::set_lang(ConvolutionGenerator::C)
      module_file_name = File::split(path.chomp(File::extname(path)))[0] + "/Mod_" + File::split(path.chomp(File::extname(path)))[1].gsub("-","_") + ".c"
      module_name = File::split(module_file_name.chomp(File::extname(module_file_name)))[1]
      module_file = File::open(module_file_name,"w+")
      ConvolutionGenerator::set_output(module_file)
      module_file.write <<EOF
#include "ruby.h"
#include <inttypes.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif
static VALUE rb_cVArray;

enum vector_type {
  VA_NONE,
  VA_CHAR,
  VA_UCHAR,
  VA_SHORT,
  VA_USHORT,
  VA_INT,
  VA_UINT,
  VA_LONG,
  VA_ULONG,
  VA_FLOAT,
  VA_DOUBLE,
  VA_ERROR
};

typedef struct _struct_varray {
  void* ptr;
  enum vector_type type;
  unsigned int n;
  size_t length;
  size_t size;
  VALUE obj;
} struct_varray;
EOF
      @procedure.header(previous_lang)
      module_file.write <<EOF
VALUE #{module_name} = Qnil;
void Init_#{module_name}();
VALUE method_run(int argc, VALUE *argv, VALUE self);
void Init_#{module_name}() {
  #{module_name} = rb_define_module("#{module_name}");
  rb_define_method(#{module_name}, "run", method_run, -1);
}
VALUE method_run(int argc, VALUE *argv, VALUE self) {
  if( argc != #{@procedure.parameters.length} )
    rb_raise(rb_eArgError, "wrong number of arguments for #{@procedure.name} (%d for #{@procedure.parameters.length})", argc);
  VALUE rb_ptr;
EOF
      argc = @procedure.parameters.length
      argv = Variable::new("argv",Real,{:dimension => [ Dimension::new(0,argc-1) ] })
      rb_ptr = Variable::new("rb_ptr",Int)
      @procedure.parameters.each { |param| 
        param_copy = param.copy
        param_copy.constant = nil
        param_copy.direction = nil
        param_copy.decl
      }
      @procedure.parameters.each_index do |i|
        param = @procedure.parameters[i]
        if not param.dimension then
          case param.type
            when Int 
              (param === FuncCall::new("NUM2INT", argv[i])).print if param.type.size == 4
              (param === FuncCall::new("NUM2LONG", argv[i])).print if param.type.size == 8
            when Real
              (param === FuncCall::new("NUM2DBL", argv[i])).print
          end
        else
          (rb_ptr === argv[i]).print
          module_file.print <<EOF
  if (TYPE(rb_ptr) == T_STRING) {
    #{param.name} = (void *) RSTRING_PTR(rb_ptr);
  } else if (CLASS_OF(rb_ptr) == rb_cVArray) {
    struct_varray *s_vary;
    Data_Get_Struct(rb_ptr, struct_varray, s_vary);
    #{param.name} = (void *) s_vary->ptr;
  } else
    rb_raise(rb_eArgError, "wrong type of argument %d", #{i});
EOF
        end
      end
      module_file.print "  #{@procedure.name}"
      module_file.print "_" if previous_lang == ConvolutionGenerator::FORTRAN
      module_file.print "("
      if(previous_lang == ConvolutionGenerator::FORTRAN) then
        params = []
        @procedure.parameters.each { |param|
          if param.dimension then
            params.push( param.name )
          else
            params.push( "&"+param.name )
          end
        }
        module_file.print params.join(", ")
      else
        module_file.print @procedure.parameters.join(", ") 
      end
      module_file.print "  );\n"
      module_file.print "  return Qnil;"
      module_file.print  "}"
      module_file.rewind
#      puts module_file.read
      module_file.close
      ConvolutionGenerator::set_lang(previous_lang)
      ConvolutionGenerator::set_output(previous_output)
      module_target = module_file_name.chomp(File::extname(module_file_name))+".o"
      module_final = module_file_name.chomp(File::extname(module_file_name))+".so"
      file module_final => [module_target, target] do
        sh "cc -shared -o #{module_final} #{module_target} #{target} -Wl,-Bsymbolic-functions -Wl,-z,relro -rdynamic -Wl,-export-dynamic #{ldflags}" 
      end
      Rake::Task[module_final].invoke
      require(module_final)
      eval "self.extend(#{module_name})"
#      self.run(32,32,"0"*32*32*8, "0"*32*32*8)
      f = File::open(target,"rb")
      @binary = StringIO::new
      @binary.write( f.read )
      f.close
      File.unlink(target)
      File.unlink(module_target)
      File.unlink(module_file_name)
      File.unlink(module_final)
    end
  end

end
