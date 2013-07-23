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
      cuda_compiler = options[:NVCC]
      cuda_compiler = "nvcc"if not cuda_compiler
      f_flags = options[:FCFLAGS]
      f_flags = "-O2 -Wall" if not f_flags
      f_flags += " -fPIC"
      f_flags += " -fno-second-underscore" if f_compiler == 'g95'
      ld_flags = options[:LDFLAGS]
      ld_flags = "" if not ld_flags
      cuda_flags = options[:NVCCFLAGS]
      cuda_flags = "-O2 -Wall" if not cuda_flags
      cuda_flags += " -fPIC"


      includes = "-I#{RbConfig::CONFIG["archdir"]}"
      includes += " -I#{RbConfig::CONFIG["rubyhdrdir"]} -I#{RbConfig::CONFIG["rubyhdrdir"]}/#{RbConfig::CONFIG["arch"]}"
      ld_flags += " -L#{RbConfig::CONFIG["libdir"]} #{RbConfig::CONFIG["LIBRUBYARG"]} -lrt"
      ld_flags += " -lcuda" if @lang == ConvolutionGenerator::CUDA
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

      rule '.o' => '.cu' do |t|
        cuda_call_string = "#{cuda_compiler} #{cudaflags} -c -o #{t.name} #{t.source}"
        runner.call(t, cuda_call_string)
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
      devices = platform.devices(type)
      if options[:device_name] then
        devices.each{ |d|
          device = d if d.name.match(options[:device_name])
        }
      else
        device = devices.first
      end
      @context = OpenCL::Context::new(nil,[device])
      program = OpenCL::Program::create_with_source(@context, [@code.string])
      opts = options[:CLFLAGS]
      program.build(:options => options[:CLFLAGS])
      @queue = OpenCL::CommandQueue::new(@context, device, OpenCL::CommandQueue::PROFILING_ENABLE)
      @kernel = OpenCL::Kernel::new( program, @procedure.name)
      run_method = <<EOF
def self.run(*args)
  raise "Wrong number of arguments \#{args.length} for #{@procedure.parameters.length}" if args.length > #{@procedure.parameters.length+1} or args.length < #{@procedure.parameters.length}
  params = []
  opts = {}
  opts = args.pop if args.length == #{@procedure.parameters.length+1}
  @procedure.parameters.each_index { |i|
    if @procedure.parameters[i].dimension then
      if @procedure.parameters[i].direction == :in and @procedure.parameters[i].direction == :out then
        params[i] = OpenCL::Buffer::new(@context, OpenCL::Mem::READ_WRITE, :size => args[i].size * args[i].element_size )
        @queue.enqueue_write_buffer(params[i], OpenCL::TRUE, args[i], :cb => args[i].size * args[i].element_size )
      elsif @procedure.parameters[i].direction == :in then
        params[i] = OpenCL::Buffer::new(@context, OpenCL::Mem::READ_ONLY, :size => args[i].size * args[i].element_size )
        @queue.enqueue_write_buffer(params[i], OpenCL::TRUE, args[i], :cb => args[i].size * args[i].element_size )
      elsif @procedure.parameters[i].direction == :out then
        params[i] = OpenCL::Buffer::new(@context, OpenCL::Mem::WRITE_ONLY, :size => args[i].size * args[i].element_size )
      else
        params[i] = OpenCL::Buffer::new(@context, OpenCL::Mem::READ_WRITE, :size => args[i].size * args[i].element_size )
      end
    else
      if @procedure.parameters[i].type.is_a?(Real) then
        params[i] = OpenCL::Half::new(args[i]) if @procedure.parameters[i].type.size == 2
        params[i] = OpenCL::Float::new(args[i]) if @procedure.parameters[i].type.size == 4
        params[i] = OpenCL::Double::new(args[i]) if @procedure.parameters[i].type.size == 8
      elsif @procedure.parameters[i].type.is_a?(Int) then
        if @procedure.parameters[i].type.signed
          params[i] = OpenCL::Char::new(args[i]) if @procedure.parameters[i].type.size == 1
          params[i] = OpenCL::Short::new(args[i]) if @procedure.parameters[i].type.size == 2
          params[i] = OpenCL::Int::new(args[i]) if @procedure.parameters[i].type.size == 4
          params[i] = OpenCL::Long::new(args[i]) if @procedure.parameters[i].type.size == 8
        else
          params[i] = OpenCL::UChar::new(args[i]) if @procedure.parameters[i].type.size == 1
          params[i] = OpenCL::UShort::new(args[i]) if @procedure.parameters[i].type.size == 2
          params[i] = OpenCL::UInt::new(args[i]) if @procedure.parameters[i].type.size == 4
          params[i] = OpenCL::ULong::new(args[i]) if @procedure.parameters[i].type.size == 8
        end
      else
        params[i] = args[i]
      end
    end
  }
  params.each_index{ |i|
    @kernel.set_arg(i, params[i])
  }
  event = @queue.enqueue_NDrange_kernel(@kernel, opts[:global_work_size], opts[:local_work_size])
  @procedure.parameters.each_index { |i|
    if @procedure.parameters[i].dimension then
      if @procedure.parameters[i].direction == :in and @procedure.parameters[i].direction == :out then
        @queue.enqueue_read_buffer(params[i], OpenCL::TRUE, :ptr => args[i], :cb => args[i].size * args[i].element_size )
      elsif @procedure.parameters[i].direction == :out then
        @queue.enqueue_read_buffer(params[i], OpenCL::TRUE, :ptr => args[i], :cb => args[i].size * args[i].element_size )
      end
    end
  }
  result = {}
  result[:start] = event.get_profiling_info(OpenCL::PROFILING_COMMAND_START).unpack("L").first
  result[:end] = event.get_profiling_info(OpenCL::PROFILING_COMMAND_END).unpack("L").first
  result[:duration] = (result[:end] - result[:start])/1000000000.0
  return result
end
EOF
    eval run_method
    return self
    end

    def build(options = {})
      return build_opencl(options) if @lang == ConvolutionGenerator::CL
      ldflags = self.setup_compiler(options)
      extension = ".c" if @lang == ConvolutionGenerator::C
      extension = ".cu" if @lang == ConvolutionGenerator::CUDA
      extension = ".f90" if @lang == ConvolutionGenerator::FORTRAN
#temporary
      c_compiler = options[:CC]
      c_compiler = "cc" if not c_compiler
      linker = options[:LD]
      linker = c_compiler if not linker
#end temporary
      source_file = Tempfile::new([@procedure.name,extension])
      path = source_file.path
      target = path.chomp(File::extname(path))+".o"
      fill_code(source_file)
      source_file.close

      previous_lang = ConvolutionGenerator::get_lang
      previous_output = $output
      ConvolutionGenerator::set_lang(ConvolutionGenerator::C)
      module_file_name = File::split(path.chomp(File::extname(path)))[0] + "/Mod_" + File::split(path.chomp(File::extname(path)))[1].gsub("-","_") + ".c"
      module_name = File::split(module_file_name.chomp(File::extname(module_file_name)))[1]
      module_file = File::open(module_file_name,"w+")
      ConvolutionGenerator::set_output(module_file)
      fill_module(module_file, module_name)
      module_file.rewind
#     puts module_file.read
      module_file.close
      ConvolutionGenerator::set_lang(previous_lang)
      ConvolutionGenerator::set_output(previous_output)
      module_target = module_file_name.chomp(File::extname(module_file_name))+".o"
      module_final = module_file_name.chomp(File::extname(module_file_name))+".so"
      file module_final => [module_target, target] do
        #puts "#{linker} -shared -o #{module_final} #{module_target} #{target} -Wl,-Bsymbolic-functions -Wl,-z,relro -rdynamic -Wl,-export-dynamic #{ldflags}"
        sh "#{linker} -shared -o #{module_final} #{module_target} #{target} -Wl,-Bsymbolic-functions -Wl,-z,relro -rdynamic -Wl,-export-dynamic #{ldflags}" 
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

    def fill_code(source_file)
      @code.rewind
      source_file.puts "#include <inttypes.h>" if @lang == ConvolutionGenerator::C or @lang == ConvolutionGenerator::CUDA
      source_file.puts "#include <cuda.h>" if @lang == ConvolutionGenerator::CUDA
      source_file.write @code.read
      if @lang == ConvolutionGenerator::CUDA then
      end
      @code.rewind
    end

    def fill_module(module_file, module_name)
      module_file.write <<EOF
#include "ruby.h"
#include <inttypes.h>
#include <time.h>
#ifdef HAVE_NARRAY_H
#include "narray.h"
#endif
EOF
      @procedure.header(@lang)
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
      module_file.print "  #{@procedure.properties[:return].type.decl} ret;\n" if @procedure.properties[:return]
      module_file.print "  VALUE stats = rb_hash_new();\n"
      module_file.print "  struct timespec start, stop;\n"
      module_file.print "  unsigned long long int duration;\n"
      module_file.print "  clock_gettime(CLOCK_REALTIME, &start);\n"
      if @lang == ConvolutionGenerator::CUDA then
        module_file.print "  duration = "
      elsif @procedure.properties[:return] then
        module_file.print "  ret = "
      end
      module_file.print "  #{@procedure.name}"
      module_file.print "_" if @lang == ConvolutionGenerator::FORTRAN
      module_file.print "("
      if(@lang == ConvolutionGenerator::FORTRAN) then
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
      if @lang == ConvolutionGenerator::CUDA then
        module_file.print "  rb_hash_aset(stats,ID2SYM(rb_intern(\"duration\")),rb_float_new((double)duration*(double)1e-3));\n"
      else
        module_file.print "  duration = (unsigned long long int)stop.tv_sec * (unsigned long long int)1000000000 + stop.tv_nsec;\n"
        module_file.print "  duration -= (unsigned long long int)start.tv_sec * (unsigned long long int)1000000000 + start.tv_nsec;\n"
        module_file.print "  rb_hash_aset(stats,ID2SYM(rb_intern(\"duration\")),rb_float_new((double)duration*(double)1e-9));\n"
      end
      if @procedure.properties[:return] then
        type_ret = @procedure.properties[:return].type
        module_file.print "  rb_hash_aset(stats,ID2SYM(rb_intern(\"return\")),rb_int_new((long long)ret));\n" if type_ret.kind_of?(Int) and type_ret.signed
        module_file.print "  rb_hash_aset(stats,ID2SYM(rb_intern(\"return\")),rb_int_new((unsigned long long)ret));\n" if type_ret.kind_of?(Int) and not type_ret.signed
        module_file.print "  rb_hash_aset(stats,ID2SYM(rb_intern(\"return\")),rb_float_new((double)ret));\n" if type_ret.kind_of?(Real)
      end
      module_file.print "  return stats;\n"
      module_file.print  "}"
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
