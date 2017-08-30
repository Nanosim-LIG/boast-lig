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
  #@cost_mem           = Int("cost_mem", :dir => :out )

  #@P_new              = Real("P_new",:dir => :inout, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  #@P_inter            = Real("P_inter",:dir => :in, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  #@P_old              = Real("P_old",:dir => :in, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  #@A_flu              = Real("A_flu",:dir => :in, :dim => [Dim(@afluSize)])

 end


 ###########################################
 ### Non-optimized version of the kernel ###
 ###########################################

 def generate

	push_env(:lang => FORTRAN) {
		@kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_boast",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize, @idx_vec_flu, @idx_mat_flu, @cost] ,:functions => nil)
     #@kernel.procedure = Procedure("fluid_inner_elt_boast",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize, @idx_vec_flu, @idx_mat_flu, @P_new, @P_inter, @P_old, @A_flu] ,:functions => nil)

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
    #i3 = Variable::new('I3', Int)
    #i4 = Variable::new('I4', Int)
    #p_aux = Variable::new("P_aux",Real, :dimension => [Dim(20),Dim(@nb_rhs)])
    #cost = Variable::new('cost', Int)

    decl i,j,i_tmp1,i_tmp2, i1,i2 #,i3,i4,p_aux

		#decl cost
		#pr cost === 0

		pr For(i,1,@Nflu_inner+@Nflusol_inner){
			pr i_tmp1 === 1
#			pr For(j,1,@idx_vec_flu[1,i]){
#				pr i1 === @idx_vec_flu[2*j,i]
#        pr i2 === @idx_vec_flu[2*j+1,i]
#        pr i_tmp2 === i_tmp1 + i2 - i1
#        #pr p_aux.slice(i_tmp1..i_tmp2,1..@nb_rhs) === @P_inter.slice(i1..i2,1..@nb_rhs)
#        pr i_tmp1 === i_tmp2 +1
#				@cost_mem === @cost_mem + 
#			}
      pr i1 === @idx_vec_flu[2,i]
      pr i2 === @idx_vec_flu[3,i]
      #pr i3 === @idx_mat_flu[i]
      #pr i4 === @idx_mat_flu[i+1] - 1
      #reshape_code ="RESHAPE(#{@A_flu.slice(i3..i4)}, (/#{i2}-#{i1}+1, #{i_tmp2}/))"
      #pr @P_new.slice(i1..i2,1..@nb_rhs) === call_matmul.call(reshape_code , p_aux.slice(1..i_tmp2,1..@nb_rhs) ) - @P_old.slice(i1..i2,1..@nb_rhs)		
			pr @cost === @cost + 2*(i2-i1+1)*20*@nb_rhs + (i2-i1+1)*@nb_rhs
		}
		close @kernel.procedure
  }


	#@kernel.cost_function = lambda { |*args| matmul_operation.cost(*args) }

   return @kernel

 end
 
end
