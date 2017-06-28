class KBoast
 attr_reader :kernel


 def initialize(options)

  @opts = {:optim_nested => 1, :optim_main => 1}
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
  @idx_vec_flu        = Int("idx_vec_flu", :dir => :in, :dim => [Dim(@vecSize_1),Dim(@vecSize_2)] )
  @idx_mat_flu        = Int("idx_mat_flu", :dir => :in, :dim => [Dim(@matSize)] )
  @P_new              = Real("P_new",:dir => :inout, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  @P_inter            = Real("P_inter",:dir => :in, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  @P_old              = Real("P_old",:dir => :in, :dim => [Dim(@pSize),Dim(@nb_rhs)])
  @A_flu              = Real("A_flu",:dir => :in, :dim => [Dim(@afluSize)])

 end


 ###########################################
 ### Non-optimized version of the kernel ###
 ###########################################

 def generate_orig

	push_env(:lang => FORTRAN) {
		@kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_boast",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @idx_vec_flu, @idx_mat_flu, @P_new, @P_inter, @P_old, @A_flu, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize] ,:functions => nil)

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
    i3 = Variable::new('I3', Int)
    i4 = Variable::new('I4', Int)
    p_aux = Variable::new("P_aux",Real, :dimension => [Dim(20),Dim(@nb_rhs)])

    decl i,j,i_tmp1,i_tmp2, i1,i2,i3,i4,p_aux

		pr For(i,1,@Nflu_inner+@Nflusol_inner){
			pr i_tmp1 === 1
			For(j,1,@idx_vec_flu[1,i]){
				pr i1 === @idx_vec_flu[2*j,i]
        pr i2 === @idx_vec_flu[2*j+1,i]
        pr i_tmp2 === i_tmp1 + i2 - i1
        pr p_aux.slice(i_tmp1..i_tmp2,1..@nb_rhs) === @P_inter.slice(i1..i2,1..@nb_rhs)
        pr i_tmp1 === i_tmp2 +1
			}
      pr i1 === @idx_vec_flu[2,i]
      pr i2 === @idx_vec_flu[3,i]
      pr i3 === @idx_mat_flu[i]
      pr i4 === @idx_mat_flu[i+1] - 1
      reshape_code ="RESHAPE(#{@A_flu.slice(i3..i4)}, (/#{i2}-#{i1}+1, #{i_tmp2}/))"
      pr @P_new.slice(i1..i2,1..@nb_rhs) === call_matmul.call(reshape_code , p_aux.slice(1..i_tmp2,1..@nb_rhs) ) - @P_old.slice(i1..i2,1..@nb_rhs)		

		}
		close @kernel.procedure
  }
   return @kernel

 end
 


 #######################################
 ### Optimized version of the kernel ###
 #######################################

 def generate

   push_env(:lang => FORTRAN) {

     @kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_boast",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @idx_vec_flu, @idx_mat_flu, @P_new, @P_inter, @P_old, @A_flu, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize] ,:functions => nil)


		#..  FUNCTIONS ..#

		register_funccall("matmul") if get_lang == FORTRAN
		#reecrire matmul pour le C @matmul = Procedure("matmul",[..,..], :return => ..){ ... }
		call_matmul = lambda{|x,y|
			if get_lang == FORTRAN
				return matmul(x,y)
	    #else
			#	return @functions[:matmul].call(x,y)
 			end
		}

		register_funccall(:modulo) if get_lang == FORTRAN


		#.. KERNEL ..#

		opn @kernel.procedure

		# local variables
 		i = Variable::new('I', Int)
 		j = Variable::new('J', Int)
 		i_tmp1 = Variable::new('I_tmp1', Int)
 		i_tmp2 = Variable::new('I_tmp2', Int)
 		i1 = Variable::new('I1', Int)
		i2 = Variable::new('I2', Int)
 		i3 = Variable::new('I3', Int)
 		i4 = Variable::new('I4', Int)
 		p_aux = Variable::new("P_aux",Real, :dimension => [Dim(20),Dim(@nb_rhs)])

 		jmax = Variable::new('jmax', Int)
 		imax = Variable::new('imax', Int)

		decl i,j,i_tmp1,i_tmp2, i1,i2,i3,i4,p_aux
		decl jmax,imax


		##.. Nested loop optimization

   	optim_1 = lambda { |unroll,i|
			loop_1 = lambda { |indice|
				pr i1 === @idx_vec_flu[2*indice,i]
				pr i2 === @idx_vec_flu[2*indice+1,i]
				pr i_tmp2 === i_tmp1 + i2 - i1
				pr p_aux.slice(i_tmp1..i_tmp2,1..@nb_rhs) === @P_inter.slice(i1..i2,1..@nb_rhs)
				pr i_tmp1 === i_tmp2 + 1.to_var
			}
		
			pr jmax === @idx_vec_flu[1,i]-(unroll-1)
			f_nested = For(j,1,jmax, step: unroll){
				unroll.times { |k|
					loop_1[j+k]
      	}
			}
			f_arr_nested = [f_nested]
			if unroll > 1 then
					f2_nested = For(j,(@idx_vec_flu[1,i]-modulo(@idx_vec_flu[1,i],unroll))+1, @idx_vec_flu[1,i]){
						loop_1[j]
					}
					f_arr_nested.push f2_nested
			end
			f_arr_nested
	 }


		##.. Main loop optimization

		optim_2 = lambda { |unroll|
			loop_2 = lambda { |indice|
				pr i_tmp1 === 1
				# Nested loop optimized - No optimization -> = optim_1[1,indice].each{ |f|  pr f }
				optim_1[@opts[:optim_nested],indice].each { |f_nested|
					pr OpenMP::Parallel(:num_threads => @opts[:omp_num_threads]){
					pr OpenMP::For()
          pr f_nested
					}
        }
				pr i1 === @idx_vec_flu[2,indice]
      	pr i2 === @idx_vec_flu[3,indice]
      	pr i3 === @idx_mat_flu[indice]
      	pr i4 === @idx_mat_flu[indice+1] - 1
				reshape_code ="RESHAPE(#{@A_flu.slice(i3..i4)}, (/#{i2}-#{i1}+1, #{i_tmp2}/))"
      	pr @P_new.slice(i1..i2,1..@nb_rhs) === call_matmul.call(reshape_code , p_aux.slice(1..i_tmp2,1..@nb_rhs) ) - @P_old.slice(i1..i2,1..@nb_rhs)
			}
			
			pr imax === @Nflu_inner+@Nflusol_inner-(unroll-1)
			f_main = For(i,1,imax, step: unroll){
				unroll.times { |k|
					loop_2[i+k]
				}
			}
			f_arr_main = [f_main]
			if unroll > 1 then
				f2_main = For(i,(@Nflu_inner+@Nflusol_inner-modulo(@Nflu_inner+@Nflusol_inner,unroll))+1, @Nflu_inner+@Nflusol_inner){
					loop_2[i]
				}
				f_arr_main.push f2_main
			end
			f_arr_main
		}



		## Procedure
		optim_2[@opts[:optim_main]].each { |f_main|
				#pr OpenMP::Parallel(:num_threads => @opts[:omp_num_threads]){
				#pr OpenMP::For()
        pr f_main
				#}
    }

		close @kernel.procedure
   }

   return @kernel
 end

###
# TODO: inverser les boucles du kernel pour voir ce que ca donne
###


end
