module libconw
  
  !> fake type to copy the documentation from 
  type, private :: doc
     integer :: fil_id
     integer :: m
     integer :: d
     integer, dimension(d) :: n
     integer, dimension(d) :: ld_sw
     integer, dimension(d) :: ld_sf
     integer, dimension(d) :: bc !< array of bc_id for any dimension
     !> array of two resolution levels, of dimension 2**d*product(nsw)
     real :: sw 
     real :: sf !< array of one resolution level doubles by wavelet transform
     real :: fr !< function in real space
     real :: vr !< real space potential
     real :: vf !< multiplication of Vr and fr
     real :: f2 !< square of the function in real space (used for the density) 
  end type doc

  contains

    !> 
    subroutine get_filter_id(filt_name,filt_id)
    end subroutine get_filter_id

    subroutine get_filter_name(filt_id,filt_name)
    end subroutine get_filter_name

    subroutine get_bc_id(bc_name,bc_id)
    end subroutine get_bc_id

    subroutine get_bc_name(bc_id,bc_name)
    end subroutine get_bc_name

    subroutine get_s0s1_input_size(d,fil_id,bc,n)
    end subroutine get_s0s1_input_size

    subroutine get_s0s1_output_size(d,fil_id,bc,n)
    end subroutine get_s0s1_output_size

    subroutine get_s0s1_work_size(d,fil_id,bc,n)
    end subroutine get_s0s1_work_size

!-----------wavelet transforms
    !> normal wavelet transform in only one direction
    !! y := a * dwt_fil |X_i| x(x_1,...,x_i,...,x_d,j)  + a_y * y
    subroutine s0s1_1d(fil_id,tp,d,idim,n,bc_id,nx,ny,narr,a,a_y,x,y)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, intent(in) :: idim
      integer, dimension(d), intent(in) :: n
      integer, intent(in) :: bc_id 
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::sf
      real(tp), dimension(*), intent(inout)  :: y !< @copydoc doc::sw
    end subroutine s0s1_1d

    !> inverse wavelet transform in only one direction
    !! y := a idwt_fil |X_i| x(x_1,...,x_i,...,x_d,j)  + a_y * y
    subroutine s1s0_1d(fil_id,tp,d,idim,n,bc_id,nx,ny,narr,a,a_y,x,y)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc_id
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::sw
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::sf
    end subroutine s1s0_1d

    !> normal wavelet transform
    !! y := a * dwt_fil (x) x + a_y * y
    subroutine s0s1_md(fil_id,tp,d,n,bc,nx,ny,narr,a,a_y,x,y,work)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc  
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::sf
      real(tp), dimension(*), intent(inout)  :: y !< @copydoc doc::sw
      real(tp), dimension(*), intent(inout) :: work
    end subroutine s0s1_md

    !> inverse wavelet transform
    !! sf := a * idwt_fil (x) sw + a_y * sf 
    subroutine s1s0_md(fil_id,tp,d,n,bc,nx,ny,narr,a,a_y,x,y,work)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc  
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::sw
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::sf
      real(tp), dimension(*), intent(inout) :: work
    end subroutine s1s0_md

!------ end of wavelet transforms

!-------filters in one direction
    !> One dimensional convolution of the filter in the direction specified by the user
    !! useful to calculate for example the vector component of the gradient
    !! y := a * fil |X_i| x(x_1,...,x_i,...,x_d,i) + a_x * x + a_y * y
    subroutine s0s0_1d(fil_id,tp,d,idim,n,bc_id,nx,ny,narr,a,a_x,a_y,x,y)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, intent(in) :: idim
      integer, dimension(d), intent(in) :: n
      integer, intent(in) :: bc_id  !<only of the convolved dimension
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_x
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
    end subroutine s0s0_1d

    !> sum of multi dimensional convolutions in the different direction
    !! y = sum_i (a_i * fil |X_i| x(x_1,...,x_i,...,x_d,i))+ a_x * x + a_y * y
    subroutine s0s0_1ds(fil_id,tp,d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,work)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, intent(in) :: idim
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), dimension(d), intent(in) :: a
      real(tp), intent(in) :: a_x
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
      real(tp), dimension(*), intent(inout) :: work !< @copydoc doc::work
    end subroutine s0s0_1ds

    subroutine s0s0_1ds_dot(fil_id,tp,d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,dot_in,work)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, intent(in) :: idim
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), dimension(d), intent(in) :: a
      real(tp), intent(in) :: a_x
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
      real(tp), dimension(d), intent(out) :: dot_in !<array of scalar products <x |out_i>
      real(tp), dimension(*), intent(inout) :: work !< @copydoc doc::work
    end subroutine s0s0_1ds_dot

!-------end filters in one direction


!---- onelevel separable kernel

    !> application of aq separable multi-dimensional convolution on the input array
    !! y := a * fil |X| x + a_x * x + a_y * y
    subroutine s0s0_md(fil_id,tp,d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,work)
      implicit none
      integer, intent(in) :: fil_id
      integer, intent(in) :: tp
      integer, intent(in) :: d
      integer, intent(in) :: idim
      integer, dimension(d), intent(in) :: n
      integer, dimension(d), intent(in) :: bc
      integer, dimension(d), intent(in) :: nx
      integer, dimension(d), intent(in) :: ny
      integer, intent(in) :: narr
      real(tp), intent(in) :: a
      real(tp), intent(in) :: a_x
      real(tp), intent(in) :: a_y
      real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
      real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
      real(tp), dimension(*), intent(inout) :: work !< @copydoc doc::work
    end subroutine s0s0_md

!------------- end onelevel separable kernel


!------ dotproduct to be discussed
!!$    call s0s0_md_x('I',fil_id,tp,d,n,bc,nx,ny,narr,v,a_x,a_y,x,y,work,dotp)
!!$
!!$    call s0s0_md(fil_id,tp,d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,work)
!!$
!!$    subroutine s0vs0_1d_x('N','N',,tp,d,idim,n,bc_id,nx,ny,narr,a,a_y,x,y,v_x,v_y)
!!$
!!$    y = V(x) MF x
!!$
!!$    y = V(x) MF x ;  dotp_out =  (x MF^t) V (MF x) = \int dx psir(x) V(x) psir(x) = \int MFx(x) V(x) MFx(x) = <out| V out> = <MF|V| MF in> -> Specialized
!!$!!                  ;  dotp_in =  x MF x = \int dx \psi(x) psir(x) = <in|out> NO -> different basis
!!$
!!$    y = K x + V x ;   dotp_in =  x ( K x + V x) = \int dx \psi(x) ( K + V) psi(x) = <in|out>  -> OK
!!$ !!                 ; ? dotp_out = <out|out> = < xK | Kx> NO -> ddot

    !get the id of the wavelet transform
    call get_filter_id('SYM8_LP',dwt_fil)
    call get_filter_id('SYM8_D2',kin_fil)
    call get_filter_id('SYM8_MD',md_fil)
    call get_filter_id('SYM8_MI',mi_fil)

    call get_bc_id('G',grow) !grow
    call get_bc_id('S',shrink) !shrink
    call get_bc_id('F',nogrow) !free, nogrow


    call get_s0s1_input_size(3,[n1,n2,n3],dwt_fil,3*[grow],nx)

    !call get_s0s1_output_size(3,[n1,n2,n3],dwt_fil,3*[fbc_id],ny)

    !for the workspace
    call get_s1s0_work_size(3,[n1,n2,n3],dwt_fil,3*[grow],nw_w)
    call get_s1s0_input_size(3,[n1,n2,n3],dwt_fil,3*[grow],ns1)
    call get_s1s0_output_size(3,[n1,n2,n3],dwt_fil,3*[grow],ns0)
    !call get_s0s0_1ds_output_size(3,2*ns0,kwt_fil,3*[nogrow],ntpsi) !nogrow
    call get_s0s0_md_output_size(3,2*ns0,md_fil,3*[grow],npsir)
    call get_s0s0_1ds_work_size(3,2*ns0,kwt_fil,3*[nogrow],nw_k)
    call get_s0s0_md_work_size(3,2*ns0,md_fil,3*[grow],nw_mf)
    nwork=max(nw_w8,nw_k,mw_mf)
    
    !npsifscf >= ns0
    !npsig >= ns1

    !idwt
    call s1s0_md(dwt_fil,8,3,[n1,n2,n3],3*[grow],ns1,ns0,1,&
         1.0_8,0.0_8,&
         psig,psi_fscf,work)

    call s0s0_md(md_fil,8,3,2*ns0,3*[grow],2*ns0,npsir,1,&
         1.0_8,0.0_8,0.0_8,&
         psi_fscf,psir,work)

    !density construction
    !rho=rho+psir**2
    
    !potential multiplication
    !psir=v*psir
    !and other treatments

    !inverse magic filter
    call s0s0_md(mi_fil,8,3,2*ns0,3*[shrink],npsir,2*ns0,1,&
         1.0_8,0.0_8,0.0_8,&
         psir,hpsi,work)

    call s0s0_1ds_dot(kin_fil,8,3,2*ns0,3*[nogrow],2*ns0,2*ns0,1,&
         -0.5_8/[hx,hy,hz],0.0_8,1.0_8,&
         psi_fscf,hpsi,ekin,work)


    !dwt 
    call s0s1_md(dwt_fil,8,3,[n1,n2,n3],3*[shrink],ns0,ns1,1,&
         1.0_8,0.0_8,&
         hpsi,psig,work)


end module libconvw

!->>>>>>>>>>>>lighter version for alpha release
!-----------wavelet transforms
!> normal wavelet transform in only one direction
!! y := a * dwt_fil |X_i| x(x_1,...,x_i,...,x_d,j)  + a_y * y
subroutine s0s1_1d_sym8_lp_d(d,idim,n,bc_id,nx,ny,narr,a,a_y,x,y)
  implicit none
  integer, intent(in) :: d
  integer, intent(in) :: idim
  integer, dimension(d), intent(in) :: n
  integer, intent(in) :: bc_id 
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::sf
  real(kind=8), dimension(*), intent(inout)  :: y !< @copydoc doc::sw
end subroutine s0s1_1d_sym8_lp_d

!> inverse wavelet transform in only one direction
!! y := a idwt_fil |X_i| x(x_1,...,x_i,...,x_d,j)  + a_y * y
subroutine s1s0_1d_sym8_lp_d(d,idim,n,bc_id,nx,ny,narr,a,a_y,x,y)
  implicit none
  integer, intent(in) :: d
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc_id
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::sw
  real(kind=8), dimension(*), intent(inout) :: y !< @copydoc doc::sf
end subroutine s1s0_1d_sym8_lp_d

!> normal wavelet transform
!! y := a * dwt_fil (x) x + a_y * y
subroutine s0s1_sym8_lp_d(d,n,bc,nx,ny,narr,a,a_y,x,y,work)
  implicit none
  integer, intent(in) :: d
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc  
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::sf
  real(kind=8), dimension(*), intent(inout)  :: y !< @copydoc doc::sw
  real(kind=8), dimension(*), intent(inout) :: work
end subroutine s0s1_sym8_lp_d

!> inverse wavelet transform
!! sf := a * idwt_fil (x) sw + a_y * sf 
subroutine s1s0_sym8_lp_d(d,n,bc,nx,ny,narr,a,a_y,x,y,work)
  implicit none
  integer, intent(in) :: d
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc  
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::sw
  real(kind=8), dimension(*), intent(inout) :: y !< @copydoc doc::sf
  real(kind=8), dimension(*), intent(inout) :: work
end subroutine s1s0_sym8_lp_d

!------ end of wavelet transforms

!-------filters in one direction
!> One dimensional convolution of the filter in the direction specified by the user
!! useful to calculate for example the vector component of the gradient
!! y := a * fil |X_i| x(x_1,...,x_i,...,x_d,i) + a_x * x + a_y * y
subroutine s0s0_1d_sym8_md_d(d,idim,n,bc_id,nx,ny,narr,a,a_x,a_y,x,y)
  implicit none
  integer, intent(in) :: d
  integer, intent(in) :: idim
  integer, dimension(d), intent(in) :: n
  integer, intent(in) :: bc_id  !<only of the convolved dimension
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(tp), intent(in) :: a
  real(tp), intent(in) :: a_x
  real(tp), intent(in) :: a_y
  real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
  real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
end subroutine s0s0_1d_sym8_md_d

subroutine s0s0_1d_sym8_mi_d
end subroutine s0s0_1d_sym8_mi_d
!> sum of multi dimensional convolutions in the different direction
!! y = sum_i (a_i * fil |X_i| x(x_1,...,x_i,...,x_d,i))+ a_x * x + a_y * y
subroutine s0s0_1ds_sym8_d2_d(d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,work)
  implicit none
  integer, intent(in) :: d
  integer, intent(in) :: idim
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), dimension(d), intent(in) :: a
  real(kind=8), intent(in) :: a_x
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::x
  real(kind=8), dimension(*), intent(inout) :: y !< @copydoc doc::y
  real(kind=8), dimension(*), intent(inout) :: work !< @copydoc doc::work
end subroutine s0s0_1ds_sym8_d2_d

subroutine s0s0_1ds_sym8_d1_d
end subroutine s0s0_1ds_sym8_d1_d

subroutine s0s0_1ds_dot_sym8_d2_d(d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,dot_in,work)
  implicit none
  integer, intent(in) :: d
  integer, intent(in) :: idim
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(tp), dimension(d), intent(in) :: a
  real(tp), intent(in) :: a_x
  real(tp), intent(in) :: a_y
  real(tp), dimension(*), intent(in) :: x !< @copydoc doc::x
  real(tp), dimension(*), intent(inout) :: y !< @copydoc doc::y
  real(tp), dimension(d), intent(out) :: dot_in !<array of scalar products <x |out_i>
  real(tp), dimension(*), intent(inout) :: work !< @copydoc doc::work
end subroutine s0s0_1ds_dot_sym8_d2_d

!-------end filters in one direction


!---- onelevel separable kernel

!> application of aq separable multi-dimensional convolution on the input array
!! y := a * fil |X| x + a_x * x + a_y * y
subroutine s0s0_sym8_md_d(d,n,bc,nx,ny,narr,a,a_x,a_y,x,y,work)
  implicit none
  integer, intent(in) :: fil_id
  integer, intent(in) :: tp
  integer, intent(in) :: d
  integer, intent(in) :: idim
  integer, dimension(d), intent(in) :: n
  integer, dimension(d), intent(in) :: bc
  integer, dimension(d), intent(in) :: nx
  integer, dimension(d), intent(in) :: ny
  integer, intent(in) :: narr
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: a_x
  real(kind=8), intent(in) :: a_y
  real(kind=8), dimension(*), intent(in) :: x !< @copydoc doc::x
  real(kind=8), dimension(*), intent(inout) :: y !< @copydoc doc::y
  real(kind=8), dimension(*), intent(inout) :: work !< @copydoc doc::work
end subroutine s0s0_sym8_md_d

subroutine s0s0_sym8_mi_d
end subroutine s0s0_sym8_mi_d
  

!------------- end onelevel separable kernel
