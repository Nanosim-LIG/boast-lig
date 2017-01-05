
require_relative './Parameters.rb'
class Subroutine
  include Parameters
  attr_accessor :code
  def initialize(name,options,functions)
    opts = {:vector_length => 1, :unroll => false, :inline => :included}
    opts.update(options)
    @name = name
    @functions = functions
    @args = []
    @code = nil
    @vector_length = options[:vector_length]
    @usage = opts[:inline]
    @unroll = opts[:unroll]
    @ndime = opts[:ndime]
    @nb_original_param = instance_variables.length + 1
  end
  def construct(block)
    if @usage == :included then
      @code = block
    else
      functions = @functions.nil? ? nil : @functions.values
      inlined = @usage == :inlined ? true : false
      @code = Procedure(@name, @args, :inline => inlined, :functions => functions, &block)
    end
  end
  def call
    if @usage == :included then
      @code.call
    else
      eval "pr @code.call( #{instance_variables[@nb_original_param...instance_variables.length].join(',')} )"
    end
  end
end
