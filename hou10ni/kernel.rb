SUBROUTINE fluid_inner_elt(Nflu_inner, Nflusol_inner, nb_rhs, idx_vec_flu, idx_mat_flu, P_new, P_inter, P_old, A_flu, afluSize, pSize, vecSize_1, vecSize_2, matSize)
  INTEGER, PARAMETER :: dq=8 !< quadruple precision
  INTEGER, PARAMETER :: dp=4  !< double precision

 !! Parameters
 INTEGER, INTENT(in) :: Nflu_inner
 INTEGER, INTENT(in) :: Nflusol_inner
 INTEGER, INTENT(in) :: nb_rhs
 INTEGER, INTENT(in) :: idx_vec_flu(vecSize_1,vecSize_2)
 INTEGER, INTENT(in) :: idx_mat_flu(matSize)
 REAL(kind=dp), INTENT(inout) :: P_new(pSize,nb_rhs)
 REAL(kind=dp), INTENT(in) :: P_inter(pSize,nb_rhs)
 REAL(kind=dp), INTENT(in) :: P_old(pSize,nb_rhs)
 REAL(kind=dp), INTENT(in) :: A_flu(afluSize)

 !! Variables declaration
 INTEGER :: I,J,I_tmp1,I_tmp2, I1, I2, I3, I4
 !! P_aux allocation
 !REAL(kind=dp),ALLOCATABLE :: P_aux(:,:)
 !ALLOCATE(P_aux((Nfacesperelem+1)*Nphi(MAXVAL(degree_tri_loc)),nb_rhs))
 !! To evaluate P_aux allocation
 REAL(kind=dp) :: P_aux(20,1)

 DO I=1,Nflu_inner+Nflusol_inner
  I_tmp1=1
  DO J=1,idx_vec_flu(1,I)
    I1=idx_vec_flu(2*J,I)
    I2=idx_vec_flu(2*J+1,I)
    I_tmp2=I_tmp1+I2-I1
    P_aux(I_tmp1:I_tmp2,1:nb_rhs)=P_inter(I1:I2,1:nb_rhs)
    I_tmp1=I_tmp2+1
  END DO
  I1=idx_vec_flu(2,I)
  I2=idx_vec_flu(3,I)
  I3=idx_mat_flu(I)
  I4=idx_mat_flu(I+1)-1

  P_new(I1:I2,1:nb_rhs)=MATMUL(RESHAPE(A_flu(I3:I4), (/I2-I1+1, I_tmp2/)),P_aux(1:I_tmp2,1:nb_rhs))-P_old(I1:I2,1:nb_rhs)
 END DO

END SUBROUTINE fluid_inner_elt
