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
  codeF = ""
  codeF =<<EOF
 subroutine fluid_inner_elt_ref(Nflu_inner, Nflusol_inner, nb_rhs, idx_vec_flu, idx_mat_flu, P_new, P_inter, P_old, A_flu, afluSize,pSize, vecSize_1,vecSize_2,matSize)
 
 integer,     parameter  :: dp    = 4               ! 4-byte integer
 integer,     parameter  :: dq    = 8               ! Double precision

 integer, intent(in) :: Nflu_inner
 integer, intent(in) :: Nflusol_inner
 integer, intent(in) :: nb_rhs
 integer, intent(in) :: afluSize
 integer, intent(in) :: pSize
 integer, intent(in) :: vecSize_1
 integer, intent(in) :: vecSize_2
 integer, intent(in) :: matSize
 integer, intent(in) :: idx_vec_flu(vecSize_1,vecSize_2) 
 integer, intent(in) :: idx_mat_flu(matSize) 
 real(dq), intent(inout) :: P_new(pSize,nb_rhs) 
 real(dq), intent(in) :: P_inter(pSize,nb_rhs) 
 real(dq), intent(in) :: P_old(pSize,nb_rhs) 
 real(dq), intent(in) :: A_flu(afluSize)  

 integer(dp) :: I
 integer(dp) :: J
 integer(dp) :: I1
 integer(dp) :: I2
 integer(dp) :: I3
 integer(dp) :: I4
 integer(dp) :: I_tmp1
 integer(dp) :: I_tmp2
 real(dq) :: P_aux(20,nb_rhs)


  do I=1,Nflu_inner+Nflusol_inner
   I_tmp1=1
   do J=1, idx_vec_flu(1,I)
    I1=idx_vec_flu(2*J,I)
    I2=idx_vec_flu(2*J+1,I)
    I_tmp2=I_tmp1+I2-I1
    P_aux(I_tmp1:I_tmp2,1:nb_rhs)=P_inter(I1:I2,1:nb_rhs)
    I_tmp1=I_tmp2+1
   end do
   I1=idx_vec_flu(2,I)
   print *,'I1 = ',I1
   I2=idx_vec_flu(3,I)
   print *,'I2 = ',I2
   I3=idx_mat_flu(I)
   print *,'I3 = ',I3
   I4=idx_mat_flu(I+1)-1
   print *,'I4 = ',I4
   print *, 'A_flu = ', RESHAPE(A_flu(I3:I4), (/I2-I1+1, I_tmp2/))
   print *, 'P_aux = ', P_aux(1:I_tmp2,1:nb_rhs)
   print *, 'P_old = ', P_old(I1:I2,1:nb_rhs) 
   P_new(I1:I2,1:nb_rhs)=MATMUL(RESHAPE(A_flu(I3:I4), (/I2-I1+1, I_tmp2/)),P_aux(1:I_tmp2,1:nb_rhs))-P_old(I1:I2,1:nb_rhs)
  end do
end subroutine fluid_inner_elt_ref
EOF

   push_env(:lang => FORTRAN) {
     @kernel = CKernel:: new()
     @kernel.procedure = Procedure("fluid_inner_elt_ref",[@Nflu_inner, @Nflusol_inner, @nb_rhs, @idx_vec_flu, @idx_mat_flu, @P_new, @P_inter, @P_old, @A_flu, @afluSize, @pSize,@vecSize_1,@vecSize_2,@matSize] ,:functions => nil)
     get_output.print codeF
   }
   return @kernel
 end
end
