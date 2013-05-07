require "./BOAST.rb"
require 'narray'
module ConvolutionGenerator
  def ConvolutionGenerator::magicfilter_per_ref
    lang = ConvolutionGenerator::get_lang
    ConvolutionGenerator::set_lang(ConvolutionGenerator::FORTRAN)
    kernel = CKernel::new
    kernel.lang = ConvolutionGenerator::FORTRAN
    function_name = "magicfilter_per_ref"
    n = Variable::new("n",Int,{:direction => :in, :signed => false})
    ndat = Variable::new("ndat",Int,{:direction => :in, :signed => false})
    dim_in_min = 0
    dim_in_max = n-1
    dim_out_min = 0
    dim_out_max = n-1
    x = Variable::new("x",Real,{:direction => :in, :dimension => [ Dimension::new(dim_in_min, dim_in_max), Dimension::new(ndat) ] })
    y = Variable::new("y",Real,{:direction => :out, :dimension => [ Dimension::new(ndat), Dimension::new(dim_out_min, dim_out_max) ] })
    p = Procedure::new(function_name, [n,ndat,x,y])
    kernel.code.print(File::read("magicfilter_refs.f90"))

    kernel.procedure = p
    ConvolutionGenerator::set_lang(lang)
    return kernel
  end
