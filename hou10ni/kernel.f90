SUBROUTINE fluid_inner_elt_ref(Nflu_inner, Nflusol_inner, nb_rhs,afluSize, pSize, vecSize_1, vecSize_2, matSize, idx_vec_flu, idx_mat_flu, P_new, P_inter, P_old, A_flu)
  implicit none
  INTEGER, PARAMETER :: dq=8 !< quadruple precision
  INTEGER, PARAMETER :: dp=4  !< double precision

  !! Parameters
  integer, intent(in) :: Nflu_inner
  integer, intent(in) :: Nflusol_inner
  integer, intent(in) :: nb_rhs
  integer, intent(in) :: afluSize
  integer, intent(in) :: pSize
  integer, intent(in) :: vecSize_1
  integer, intent(in) :: vecSize_2
  integer, intent(in) :: matSize
  INTEGER, INTENT(in), DIMENSION(vecSize_1,vecSize_2) :: idx_vec_flu
  INTEGER, INTENT(in), DIMENSION(matSize) :: idx_mat_flu
  REAL(kind=dp), INTENT(inout), DIMENSION(pSize,nb_rhs) :: P_new
  REAL(kind=dp), INTENT(in), DIMENSION(pSize,nb_rhs) :: P_inter
  REAL(kind=dp), INTENT(in), DIMENSION(pSize,nb_rhs) :: P_old
  REAL(kind=dp), INTENT(in), DIMENSION(afluSize) :: A_flu

  !! Variables declaration
  INTEGER :: I,J,I_tmp1,I_tmp2, I1, I2, I3, I4
  REAL(kind=dp), DIMENSION(20,nb_rhs) :: P_aux

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

END SUBROUTINE fluid_inner_elt_ref
