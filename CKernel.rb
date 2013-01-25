require './Algorithm.rb'
require 'stringio'
require 'rubygems'
require 'rake'
require 'tempfile'
require 'rbconfig'
require 'systemu'

module ConvolutionGenerator
  $verbose = false
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
    def setup_compiler(options = {})
      Rake::Task::clear
      verbose = options[:verbose]
      verbose = $verbose if not verbose
      Rake::verbose(verbose)
      Rake::FileUtilsExt.verbose_flag=verbose
      f_compiler = options[:FC]
      f_compiler = "gfortran" if not f_compiler
      c_compiler = options[:CC]
      c_compiler = "cc" if not c_compiler
      cxx_compiler = options[:CXX]
      cxx_compiler = "g++" if not cxx_compiler
      f_flags = options[:FCFLAGS]
      f_flags = "-O2 -Wall" if not f_flags
      f_flags += " -fPIC"
      f_flags += " -fno-second-underscore" if f_compiler == 'g95'
      ld_flags = options[:LDFLAGS]
      ld_flags = "" if not ld_flags


      includes = "-I#{RbConfig::CONFIG["archdir"]}"
      includes += " -I#{RbConfig::CONFIG["rubyhdrdir"]} -I#{RbConfig::CONFIG["rubyhdrdir"]}/#{RbConfig::CONFIG["arch"]}"
      ld_flags += " -L#{RbConfig::CONFIG["libdir"]} #{RbConfig::CONFIG["LIBRUBYARG"]} -lrt"
      narray_path = nil
      begin
        spec = Gem::Specification::find_by_name('narray')
        narray_path = spec.full_gem_path
      rescue Gem::LoadError => e
      rescue NoMethodError => e
        spec = Gem::available?('narray')
        if spec then
          require 'narray' 
          narray_path = Gem.loaded_specs['narray'].full_gem_path
        end
      end
      includes += " -I#{narray_path}" if narray_path
      cflags = "-O2 -Wall -fPIC #{includes}"
      cxxflags = String::new(cflags)
      cflags += " -DHAVE_NARRAY_H" if narray_path
      cflags += options[:CFLAGS] if options[:CFLAGS]
      fcflags = f_flags

      runner = lambda { |t, call_string|
        if verbose then
          sh call_string
        else
          status, stdout, stderr = systemu call_string
          if not status.success? then
            puts stderr
            fail "#{t.source}: compilation failed"
          end
          status.success?
        end
      }

      rule '.o' => '.c' do |t|
        c_call_string = "#{c_compiler} #{cflags} -c -o #{t.name} #{t.source}"
        runner.call(t, c_call_string)
      end

      rule '.o' => '.f90' do |t|
        f_call_string = "#{f_compiler} #{fcflags} -c -o #{t.name} #{t.source}"
        runner.call(t, f_call_string)
      end

      rule '.o' => '.cpp' do |t|
        cxx_call_string = "#{cxx_compiler} #{cxxflags} -c -o #{t.name} #{t.source}"
        runner.call(t, cxx_call_string)
      end

      return ld_flags
    end

    def build_opencl(options)
      require 'opencl'
      platform = nil
      platforms = OpenCL::Platform::get_platforms
      if options[:platform_vendor] then
        platforms.each{ |p|
          platform = p if p.vendor.match(options[:platform_vendor])
        }
      else
        platform = platforms.first
      end
      device = nil
      type = options[:device_type] ? options[:device_type] : OpenCL::Device::TYPE_ALL
      devices = OpenCL::Device.get_devices(platform, type)
      if options[:device_name] then
        devices.each{ |d|
          device = d if d.name.match(options[:device_name])
        }
      else
        device = devices.first
      end
      context = OpenCL::Context::new(nil,[device])
      program = OpenCL::Program::create_with_source(context, [@code.string])
      program.build
      @queue = OpenCL::CommandQueue::new(context, device, OpenCL::CommandQueue::PROFILING_ENABLE)
      @kernel = OpenCL::Kernel::new( program, @procedure.name)
#      run_method = "def run("
#      @procedure.parameters.each_index do |i|
    end

    def build(options = {})
      return build_opencl(options) if @lang == ConvolutionGenerator::OpenCL
      ldflags = self.setup_compiler(options)
      extension = ".c" if @lang == ConvolutionGenerator::C
      extension = ".f90" if @lang == ConvolutionGenerator::FORTRAN
#temporary
      c_compiler = options[:CC]
      c_compiler = "cc" if not c_compiler
#end temporary
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
#include <time.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif
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
  } else if ( IsNArray(rb_ptr) ) {
    struct NARRAY *n_ary;
    Data_Get_Struct(rb_ptr, struct NARRAY, n_ary);
    #{param.name} = (void *) n_ary->ptr;
  } else
    rb_raise(rb_eArgError, "wrong type of argument %d", #{i});
EOF
        end
      end
      module_file.print "  VALUE stats = rb_hash_new();\n"
      module_file.print "  struct timespec start, stop;\n"
      module_file.print "  unsigned long long int duration;\n"
      module_file.print "  clock_gettime(CLOCK_REALTIME, &start);\n"
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
      module_file.print "  clock_gettime(CLOCK_REALTIME, &stop);\n"
      module_file.print "  duration = (unsigned long long int)stop.tv_sec * (unsigned long long int)1000000000 + stop.tv_nsec;\n"
      module_file.print "  duration -= (unsigned long long int)start.tv_sec * (unsigned long long int)1000000000 + start.tv_nsec;\n"
      module_file.print "  rb_hash_aset(stats,rb_str_new2(\"duration\"),rb_float_new((double)duration*(double)1e-9));\n"
      module_file.print "  return stats;\n"
      module_file.print  "}"
      module_file.rewind
#     puts module_file.read
      module_file.close
      ConvolutionGenerator::set_lang(previous_lang)
      ConvolutionGenerator::set_output(previous_output)
      module_target = module_file_name.chomp(File::extname(module_file_name))+".o"
      module_final = module_file_name.chomp(File::extname(module_file_name))+".so"
      file module_final => [module_target, target] do
        sh "#{c_compiler} -shared -o #{module_final} #{module_target} #{target} -Wl,-Bsymbolic-functions -Wl,-z,relro -rdynamic -Wl,-export-dynamic #{ldflags}" 
      end
      Rake::Task[module_final].invoke
      require(module_final)
      eval "self.extend(#{module_name})"
      f = File::open(target,"rb")
      @binary = StringIO::new
      @binary.write( f.read )
      f.close
      File.unlink(target)
      File.unlink(module_target)
      File.unlink(module_file_name)
      File.unlink(module_final)
      return self
    end
    def method_missing(meth, *args, &block)
     if meth.to_s == "run" then
       self.build
       self.run(*args,&block)
     else
       super
     end
    end
  end
end
