class KOssRef
 attr_reader :kernel

 def initialize(options)

  @opts = {:preprocessor => false, :unroll => false, :inline => :included}
  @opts.update(options)

  # Parameters 
  @Nflu_inner         = Int("Nflu_inner", :dir => :in )
  @Nflusol_inner      = Int("Nflusol_inner", :dir => :in )
  @nb_rhs             = Int("nb_rhs", :dir => :in )
  @afluSize           = Int("afluSize", :dir => :in )
  @pSize              = Int("pSize", :dir => :in )
  @vecSize_1          = Int("vecSize_1", :dir => :in )
  @vecSize_2          = Int("vecSize_2", :dir => :in )
  @matSize            = Int("matSize", :dir => :in )

 @idx_vec_flu        = Int("idx_vec_flu", :dir => :in, :dim => [@vecSize_1,@vecSize_2] )
 @idx_mat_flu        = Int("idx_mat_flu", :dir => :in, :dim => [@matSize] )
 @P_new              = Real("P_new",:dir => :inout, :dim => [@pSize,@nb_rhs])
 @P_inter            = Real("P_inter",:dir => :in, :dim => [@pSize,@nb_rhs])
 @P_old              = Real("P_old",:dir => :in, :dim => [@pSize,@nb_rhs])
 @A_flu              = Real("A_flu",:dir => :in, :dim => [@afluSize])

 end
 
 def generate
  codeF = File::read("./kernel.f90")

   push_env(:lang => FORTRAN) {
     @kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_ref",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @idx_vec_flu, @idx_mat_flu, @P_new, @P_inter, @P_old, @A_flu, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize] ,:functions => nil)
     get_output.print codeF
   }
   return @kernel
 end
end
