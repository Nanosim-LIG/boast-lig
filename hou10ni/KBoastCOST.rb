class KBoastCOST
 attr_reader :kernel


 def initialize(options)


  # Parameters 
  @Nflu_inner         = Int("Nflu_inner", :dir => :in )
  @Nflusol_inner      = Int("Nflusol_inner", :dir => :in )
  @nb_rhs             = Int("nb_rhs", :dir => :in )
  @afluSize           = Int("afluSize", :dir => :in )
  @pSize              = Int("pSize", :dir => :in )
  @vecSize_1          = Int("vecSize_1", :dir => :in )
  @vecSize_2          = Int("vecSize_2", :dir => :in )
  @matSize            = Int("matSize", :dir => :in )
  @idx_vec_flu        = Int("idx_vec_flu", :dir => :in, :dim => [Dim(@vecSize_1),Dim(@vecSize_2)] )
  @idx_mat_flu        = Int("idx_mat_flu", :dir => :in, :dim => [Dim(@matSize)] )
  @cost            		= Int("cost", :dir => :out )
 end


 ##################################################
 ###  Version of the kernel without computation ###
 ##################################################

 def generate

	push_env(:lang => FORTRAN) {
		@kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_boast_cost",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize, @idx_vec_flu, @idx_mat_flu, @cost] ,:functions => nil)

    #..  FUNCTIONS ..#
    register_funccall("matmul") if get_lang == FORTRAN
    call_matmul = lambda{|x,y|
      if get_lang == FORTRAN
        return matmul(x,y)
      end
    }
		
		#.. KERNEL ..#
    opn @kernel.procedure

    # local variables
    i = Variable::new('I', Int)
    j = Variable::new('J', Int)
    i_tmp1 = Variable::new('I_tmp1', Int)
    i_tmp2 = Variable::new('I_tmp2', Int)
    i1 = Variable::new('I1', Int)
    i2 = Variable::new('I2', Int)

    decl i,j,i_tmp1,i_tmp2, i1,i2

		pr @cost === 0

		pr For(i,1,@Nflu_inner+@Nflusol_inner){
      pr i1 === @idx_vec_flu[2,i]
      pr i2 === @idx_vec_flu[3,i]
			pr @cost === @cost + 2*(i2-i1+1)*20*@nb_rhs + (i2-i1+1)*@nb_rhs
		}
		close @kernel.procedure
  }

   return @kernel

 end
 
end
