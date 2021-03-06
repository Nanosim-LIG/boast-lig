require 'rubygems'
require 'BOAST'
require 'narray'
require "./GenericConvolution2.rb" 

def BOAST::kinetic_per_ref_optim_ekin
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "kinetic_per_ref_optim_ekin"
  n1 = BOAST::Int("n1", :dir => :in)
  n2 = BOAST::Int("n2", :dir => :in)
  n3 = BOAST::Int("n3", :dir => :in)
  hgrid = BOAST::Real("hgrid", :dir => :in, :dim => [BOAST::Dim(3)] )
  kstrten = BOAST::Real("kstrten", :dir => :out, :dim => [BOAST::Dim(3)] )
  x = BOAST::Real("x", :dir => :in,  :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3)] )
  y = BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3)] )
  
  p = BOAST::Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,kstrten])
  kernel.code.print <<EOF
subroutine kinetic_per_ref_optim_ekin(n1,n2,n3,hgrid,x,y,kstrten)
  !   applies the kinetic energy operator onto x to get y. Works for periodic BC
  !   y:=y-1/2Delta x
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)

  integer, intent(in) :: n1,n2,n3
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y
  !real(wp),intent(out)::ekin_out
  real(wp), dimension(6), intent(out) :: kstrten
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: k!,i !for non OMP case
  integer, dimension(lowfil:n1+lupfil) :: mod_arr1
  integer, dimension(lowfil:n2+lupfil) :: mod_arr2
  integer, dimension(lowfil:n3+lupfil) :: mod_arr3
  real(wp) :: ekin1,ekin2,ekin3
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  !real(wp), dimension(8,3) :: ekin_array
!$ integer :: ithread=0
!$ integer :: omp_get_thread_num

  !ekin_out=0._wp
  kstrten=0.0_gp

  !do i=1,8
  !ekin_array(i,1)=10._wp
  !ekin_array(i,2)=10._wp
  !ekin_array(i,3)=10._wp
  !end do
 !$omp parallel default(private) shared(x,y,n1,n2,n3,hgrid,kstrten,fil,mod_arr1,mod_arr2,mod_arr3)
!$  ithread = omp_get_thread_num()
  call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do k=1,14
     fil(-k,:)=fil(k,:)
  end do

  ekin2=0.0_wp
  ekin3=0.0_wp
  ekin1=0.0_wp
  call conv_kin_x(x,y,n1,n2,n3,ekin1,fil,mod_arr1)
  call conv_kin_y(x,y,n1,n2,n3,ekin2,fil,mod_arr2)
  call conv_kin_z(x,y,n1,n2,n3,ekin3,fil,mod_arr3)
  !ekin_array(ithread+1,1)=ekin1
  !ekin_array(ithread+1,2)=ekin2
  !ekin_array(ithread+1,3)=ekin3
  !$omp critical 

!yk
kstrten(1)=kstrten(1)+ekin1
kstrten(2)=kstrten(2)+ekin2
kstrten(3)=kstrten(3)+ekin3

!     ekin_out=ekin_out+ekin1+ekin2+ekin3
  !$omp end critical
  !$omp end parallel


!dee
!!$!  open(unit=97,file='check_ekin3',status='unknown')
!!$    write(197,*) '-------------------------------------------------------------------'
!!$do i=1,8
!!$    write(197,'(3(1X,e24.17))') ekin_array(i,1),ekin_array(i,2),ekin_array(i,3)
!!$end do
!  close(97)

END SUBROUTINE kinetic_per_ref_optim_ekin


subroutine conv_kin_x(x,y,n1,n2,n3,ekin,fil,mod_arr1)
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mod_arr1(lowfil:n1+lupfil)
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  integer :: ndat
  real(wp),intent(in) :: x(0:n1,(n2+1)*(n3+1))
  real(wp),intent(inout) :: y(0:n1,(n2+1)*(n3+1))
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
  integer :: i,l,j,i1

  ndat=(n2+1)*(n3+1)
  !$omp do
  do i=0,ndat/12-1
     do i1=0,n1
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           j=mod_arr1(i1+l)

           tt1=tt1+x(j,i*12+1)*fil(l,1)
           tt2=tt2+x(j,i*12+2)*fil(l,1)
           tt3=tt3+x(j,i*12+3)*fil(l,1)
           tt4=tt4+x(j,i*12+4)*fil(l,1)
           tt5=tt5+x(j,i*12+5)*fil(l,1)
           tt6=tt6+x(j,i*12+6)*fil(l,1)
           tt7=tt7+x(j,i*12+7)*fil(l,1)
           tt8=tt8+x(j,i*12+8)*fil(l,1)
           tt9 =tt9 +x(j,i*12+9 )*fil(l,1)
           tt10=tt10+x(j,i*12+10)*fil(l,1)
           tt11=tt11+x(j,i*12+11)*fil(l,1)
           tt12=tt12+x(j,i*12+12)*fil(l,1)
        end do
        y(i1,i*12+1)=y(i1,i*12+1)+tt1;    ekin=ekin+tt1*x(i1,i*12+1)
        y(i1,i*12+2)=y(i1,i*12+2)+tt2;    ekin=ekin+tt2*x(i1,i*12+2)
        y(i1,i*12+3)=y(i1,i*12+3)+tt3;    ekin=ekin+tt3*x(i1,i*12+3)
        y(i1,i*12+4)=y(i1,i*12+4)+tt4;    ekin=ekin+tt4*x(i1,i*12+4)
        y(i1,i*12+5)=y(i1,i*12+5)+tt5;    ekin=ekin+tt5*x(i1,i*12+5)
        y(i1,i*12+6)=y(i1,i*12+6)+tt6;    ekin=ekin+tt6*x(i1,i*12+6)
        y(i1,i*12+7)=y(i1,i*12+7)+tt7;    ekin=ekin+tt7*x(i1,i*12+7)
        y(i1,i*12+8)=y(i1,i*12+8)+tt8;    ekin=ekin+tt8*x(i1,i*12+8)
        y(i1,i*12+9 )=y(i1,i*12+9 )+tt9 ;    ekin=ekin+tt9 *x(i1,i*12+9 )
        y(i1,i*12+10)=y(i1,i*12+10)+tt10;    ekin=ekin+tt10*x(i1,i*12+10)
        y(i1,i*12+11)=y(i1,i*12+11)+tt11;    ekin=ekin+tt11*x(i1,i*12+11)
        y(i1,i*12+12)=y(i1,i*12+12)+tt12;    ekin=ekin+tt12*x(i1,i*12+12)
       end do
    end do
   !$omp end do

   !$omp do
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr1(i1+l)
             tt=tt+x(j   ,i)*fil(l,1)
          end do
          y(i1,i)=y(i1,i)+tt ; ekin=ekin+tt*x(i1,i)
       end do
    end do
   !$omp end do
END SUBROUTINE conv_kin_x

subroutine conv_kin_y(x,y,n1,n2,n3,ekin,fil,mod_arr2)
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  integer, intent(in) :: mod_arr2(lowfil:n2+lupfil)
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  real(wp),intent(in) :: x(0:n1,0:n2,0:n3)
  real(wp),intent(inout) :: y(0:n1,0:n2,0:n3)
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7
  integer :: l,j,i1,i2,i3

  !$omp do
  do i3=0,n3/8-1
     do i1=0,n1
        do i2=0,n2
           tt0=0.e0_wp
           tt1=0.e0_wp
           tt2=0.e0_wp
           tt3=0.e0_wp
           tt4=0.e0_wp
           tt5=0.e0_wp
           tt6=0.e0_wp
           tt7=0.e0_wp

           do l=lowfil,lupfil
              j=mod_arr2(i2+l)

              tt0=tt0+x(i1,j,i3*8+0)*fil(l,2)
              tt1=tt1+x(i1,j,i3*8+1)*fil(l,2)
              tt2=tt2+x(i1,j,i3*8+2)*fil(l,2)
              tt3=tt3+x(i1,j,i3*8+3)*fil(l,2)
              tt4=tt4+x(i1,j,i3*8+4)*fil(l,2)
              tt5=tt5+x(i1,j,i3*8+5)*fil(l,2)
              tt6=tt6+x(i1,j,i3*8+6)*fil(l,2)
              tt7=tt7+x(i1,j,i3*8+7)*fil(l,2)
           end do
           y(i1,i2,i3*8+0)=y(i1,i2,i3*8+0)+tt0;    ekin=ekin+tt0*x(i1,i2,i3*8+0)
           y(i1,i2,i3*8+1)=y(i1,i2,i3*8+1)+tt1;    ekin=ekin+tt1*x(i1,i2,i3*8+1)
           y(i1,i2,i3*8+2)=y(i1,i2,i3*8+2)+tt2;    ekin=ekin+tt2*x(i1,i2,i3*8+2)
           y(i1,i2,i3*8+3)=y(i1,i2,i3*8+3)+tt3;    ekin=ekin+tt3*x(i1,i2,i3*8+3)
           y(i1,i2,i3*8+4)=y(i1,i2,i3*8+4)+tt4;    ekin=ekin+tt4*x(i1,i2,i3*8+4)
           y(i1,i2,i3*8+5)=y(i1,i2,i3*8+5)+tt5;    ekin=ekin+tt5*x(i1,i2,i3*8+5)
           y(i1,i2,i3*8+6)=y(i1,i2,i3*8+6)+tt6;    ekin=ekin+tt6*x(i1,i2,i3*8+6)
           y(i1,i2,i3*8+7)=y(i1,i2,i3*8+7)+tt7;    ekin=ekin+tt7*x(i1,i2,i3*8+7)
        end do
     end do
  end do
   
  !$omp end do

  !$omp do
  do i3=(n3/8)*8,n3
     do i1=0,n1
        do i2=0,n2
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=mod_arr2(i2+l)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           end do
           y(i1,i2,i3)=y(i1,i2,i3)+tt
           ekin=ekin+tt*x(i1,i2,i3)
        end do
     end do
  end do
  !$omp end do
END SUBROUTINE conv_kin_y


subroutine conv_kin_z(x,y,n1,n2,n3,ekin,fil,mod_arr3)
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, parameter :: lowfil=-14,lupfil=14
  integer, intent(in) :: n1,n2,n3
  real(wp), intent(in) :: fil(lowfil:lupfil,3) 
  integer, intent(in) :: mod_arr3(lowfil:n3+lupfil)
  real(wp),intent(in) :: x((n2+1)*(n1+1),0:n3)
  real(wp),intent(inout) :: y((n2+1)*(n1+1),0:n3)
  real(wp),intent(inout) :: ekin
  real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
  integer :: ndat
  integer ::i,l,j,i3

  ndat=(n2+1)*(n1+1)
  !$omp do
  do i=0,ndat/12-1
     do i3=0,n3
        tt1=0.e0_wp
        tt2=0.e0_wp
        tt3=0.e0_wp
        tt4=0.e0_wp
        tt5=0.e0_wp
        tt6=0.e0_wp
        tt7=0.e0_wp
        tt8=0.e0_wp
        tt9 =0.e0_wp
        tt10=0.e0_wp
        tt11=0.e0_wp
        tt12=0.e0_wp

        do l=lowfil,lupfil
           j=mod_arr3(i3+l)

           !print *,'j,mod_arr',j,'aa',mod_arr3(:)

           tt1=tt1+x(i*12+1,j)*fil(l,3)
           tt2=tt2+x(i*12+2,j)*fil(l,3)
           tt3=tt3+x(i*12+3,j)*fil(l,3)
           tt4=tt4+x(i*12+4,j)*fil(l,3)
           tt5=tt5+x(i*12+5,j)*fil(l,3)
           tt6=tt6+x(i*12+6,j)*fil(l,3)
           tt7=tt7+x(i*12+7,j)*fil(l,3)
           tt8=tt8+x(i*12+8,j)*fil(l,3)
           tt9 =tt9 +x(i*12+9 ,j)*fil(l,3)
           tt10=tt10+x(i*12+10,j)*fil(l,3)
           tt11=tt11+x(i*12+11,j)*fil(l,3)
           tt12=tt12+x(i*12+12,j)*fil(l,3)
        end do

        y(i*12+1,i3)=y(i*12+1,i3)+tt1;    ekin=ekin+tt1*x(i*12+1,i3)
        y(i*12+2,i3)=y(i*12+2,i3)+tt2;    ekin=ekin+tt2*x(i*12+2,i3)
        y(i*12+3,i3)=y(i*12+3,i3)+tt3;    ekin=ekin+tt3*x(i*12+3,i3)
        y(i*12+4,i3)=y(i*12+4,i3)+tt4;    ekin=ekin+tt4*x(i*12+4,i3)
        y(i*12+5,i3)=y(i*12+5,i3)+tt5;    ekin=ekin+tt5*x(i*12+5,i3)
        y(i*12+6,i3)=y(i*12+6,i3)+tt6;    ekin=ekin+tt6*x(i*12+6,i3)
        y(i*12+7,i3)=y(i*12+7,i3)+tt7;    ekin=ekin+tt7*x(i*12+7,i3)
        y(i*12+8,i3)=y(i*12+8,i3)+tt8;    ekin=ekin+tt8*x(i*12+8,i3)
        y(i*12+9 ,i3)=y(i*12+9 ,i3)+tt9 ;    ekin=ekin+tt9*x(i*12+9 ,i3)
        y(i*12+10,i3)=y(i*12+10,i3)+tt10;    ekin=ekin+tt10*x(i*12+10,i3)
        y(i*12+11,i3)=y(i*12+11,i3)+tt11;    ekin=ekin+tt11*x(i*12+11,i3)
        y(i*12+12,i3)=y(i*12+12,i3)+tt12;    ekin=ekin+tt12*x(i*12+12,i3)
     end do
  end do
  !$omp end do

  !$omp do
  do i=(ndat/12)*12+1,ndat
     do i3=0,n3
        tt=0.e0_wp
        do l=lowfil,lupfil
           j=mod_arr3(i3+l)
           tt=tt+x(i,j)*fil(l,3)
        end do
        y(i,i3)=y(i,i3)+tt; ekin=ekin+tt*x(i,i3)
     end do
  end do
  !$omp end do
END SUBROUTINE conv_kin_z


subroutine fill_mod_arr(arr,nleft,nright,n)
  implicit none
  integer,intent(in) :: nleft,nright,n
  integer,intent(out) :: arr(nleft:nright)
  integer :: i
  
  if (nleft >= -n) then
     do i=nleft,-1
        arr(i)=n+i
     end do
  else
     do i=nleft,-1
        arr(i)=modulo(i,n)
     end do
  endif
  
  do i=max(0,nleft),min(n-1,nright)
     arr(i)=i
  end do
  
  if (nright < 2*n) then
     do i=n,nright
        arr(i)=i-n
     end do
  else
     do i=n,nright
        arr(i)=modulo(i,n)
     end do
  endif
END SUBROUTINE fill_mod_arr
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def kinetic_per_ref_optim
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "kinetic_per_ref_optim"
  n1 = BOAST::Int("n1", :dir => :in)
  n2 = BOAST::Int("n2", :dir => :in)
  n3 = BOAST::Int("n3", :dir => :in)
  hgrid = BOAST::Real("hgrid", :dir => :in, :dim => [BOAST::Dim(3)] )
  c = BOAST::Real("c", :dir => :in)
  x = BOAST::Real("x", :dir => :in,  :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3)] )
  y = BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3)] )
  
  p = BOAST::Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,c])
  kernel.code.print <<EOF
subroutine kinetic_per_ref_optim(n1,n2,n3,hgrid,x,y,c)
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)

  integer, intent(in) :: n1,n2,n3
  real(gp), intent(in) :: c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y

!  stop 'convolut_kinetic_per_c should never be called'

  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer, dimension(lowfil:n1+lupfil) :: mod_arr1   
  integer, dimension(lowfil:n2+lupfil) :: mod_arr2   
  integer, dimension(lowfil:n3+lupfil) :: mod_arr3   
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  integer :: k

  !$omp parallel default(private) shared(x,y,n1,n2,n3,c,hgrid,fil,mod_arr1,mod_arr2,mod_arr3)
  call fill_mod_arr(mod_arr1,lowfil,n1+lupfil,n1+1)
  call fill_mod_arr(mod_arr2,lowfil,n2+lupfil,n2+1)
  call fill_mod_arr(mod_arr3,lowfil,n3+lupfil,n3+1)

  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do k=1,14
     fil(-k,:)=fil(k,:)
  end do

  call conv_kin_x(mod_arr1,fil,c,x,y,n1,(n2+1)*(n3+1))   
  call conv_kin_y_new(mod_arr2,fil,x,y,n1,n2,n3)
  call conv_kin_z_new(mod_arr3,fil,x,y,n3,(n1+1)*(n2+1))
  !$omp end parallel


  contains


  subroutine conv_kin_x(mod_arr,fil0,c0,x0,y0,n,ndat)
    implicit none
    integer, intent(in) :: n,ndat
    integer, dimension(lowfil:n+lupfil) :: mod_arr   
    real(wp), dimension(lowfil:lupfil,3) :: fil0
    real(gp), intent(in) :: c0
    real(wp), intent(in) :: x0(0:n,ndat)
    real(wp), intent(out) :: y0(0:n,ndat)
    integer :: i1,i,l,j
    real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12

    !$omp do 
    do i=0,ndat/12-1
       do i1=0,n1
          tt1 =x0(i1,i*12+1)*c0
          tt2 =x0(i1,i*12+2)*c0
          tt3 =x0(i1,i*12+3)*c0
          tt4 =x0(i1,i*12+4)*c0
          tt5 =x0(i1,i*12+5)*c0
          tt6 =x0(i1,i*12+6)*c0
          tt7 =x0(i1,i*12+7)*c0
          tt8 =x0(i1,i*12+8)*c0
          tt9 =x0(i1,i*12+9 )*c0
          tt10=x0(i1,i*12+10)*c0
          tt11=x0(i1,i*12+11)*c0
          tt12=x0(i1,i*12+12)*c0

          do l=lowfil,lupfil
             j=mod_arr(i1+l)

             tt1=tt1+x0(j,i*12+1)*fil0(l,1)
             tt2=tt2+x0(j,i*12+2)*fil0(l,1)
             tt3=tt3+x0(j,i*12+3)*fil0(l,1)
             tt4=tt4+x0(j,i*12+4)*fil0(l,1)
             tt5=tt5+x0(j,i*12+5)*fil0(l,1)
             tt6=tt6+x0(j,i*12+6)*fil0(l,1)
             tt7=tt7+x0(j,i*12+7)*fil0(l,1)
             tt8=tt8+x0(j,i*12+8)*fil0(l,1)
             tt9 =tt9 +x0(j,i*12+9 )*fil0(l,1)
             tt10=tt10+x0(j,i*12+10)*fil0(l,1)
             tt11=tt11+x0(j,i*12+11)*fil0(l,1)
             tt12=tt12+x0(j,i*12+12)*fil0(l,1)
          end do
          y0(i1,i*12+1)=tt1
          y0(i1,i*12+2)=tt2
          y0(i1,i*12+3)=tt3
          y0(i1,i*12+4)=tt4
          y0(i1,i*12+5)=tt5
          y0(i1,i*12+6)=tt6
          y0(i1,i*12+7)=tt7
          y0(i1,i*12+8)=tt8
          y0(i1,i*12+9 )=tt9 
          y0(i1,i*12+10)=tt10
          y0(i1,i*12+11)=tt11
          y0(i1,i*12+12)=tt12
       end do
    end do
    !$omp end do

    !$omp do 
    do i=(ndat/12)*12+1,ndat
       do i1=0,n1
          tt=x0(i1,i)*c0
          do l=lowfil,lupfil
             j=mod_arr(i1+l)
             tt=tt+x0(j   ,i)*fil0(l,1)
          end do
          y0(i1,i)=tt
       end do
    end do
    !$omp end do
  END SUBROUTINE conv_kin_x
subroutine conv_kin_y_new(mod_arr,fil0,x0,y0,n1,n2,n3)
  implicit none
  !Arguments
  integer, parameter :: wp=kind(1.0d0)
  integer, parameter :: lowfil=-14,lupfil=14,ncache=4*1024
  integer, parameter :: unrol=4
  integer, intent(in) :: n1,n2,n3
  integer, dimension(lowfil:n2+lupfil), intent(in) :: mod_arr
  real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
  real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x0
  real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y0
  integer :: i1,i2,i3,l,j,j3,nol,ndim,k
  real(wp) :: tt,tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7
  real(wp), dimension(ncache) :: zcache

  !number of lines in n1 to be treated simultaneously 
  ndim=(n2-lowfil+lupfil+1)
  nol=min(ncache/(unrol*ndim),n1) 

  if (nol <1) stop 'nol'

  
  !$omp do
  do i3=0,n3/unrol-1
     do i1=0,(n1/nol)*nol-1,nol
        k=1
        !first boundary
        do j=0,nol-1
           do i2=1,-lowfil
              do j3=0,unrol-1
                 !zcache(i2+j*ndim+j3*ndim*nol)=&
                 zcache(k)=&
                      x0(i1+j,i2+lowfil-1+n2+1,i3*unrol+j3)
                 k=k+1
              end do
           end do
           do i2=-lowfil+1,n2-lowfil+1
              do j3=0,unrol-1
                 !zcache(i2+j*ndim+j3*ndim*nol)=&
                 zcache(k)=&
                      x0(i1+j,i2+lowfil-1,i3*unrol+j3)
                 k=k+1
              end do
           end do
           do i2=n2-lowfil+2,n2-lowfil+1+lupfil
              do j3=0,unrol-1
                 !zcache(i2+j*ndim+j3*ndim*nol)=&
                 zcache(k)=&
                      x0(i1+j,i2+lowfil-1-n2-1,i3*unrol+j3)
                 k=k+1
              end do
           end do
        end do
!zcache=zcache(unrol,ndim,nol)
        do j=0,nol-1
           k=1+j*unrol*ndim
           do i2=0,n2
              tt0=0.e0_wp
              tt1=0.e0_wp
              tt2=0.e0_wp
              tt3=0.e0_wp
!!$              tt4=0.e0_wp
!!$              tt5=0.e0_wp
!!$              tt6=0.e0_wp
!!$              tt7=0.e0_wp

              do l=lowfil,lupfil
                 tt0=tt0+zcache(k+0)*fil0(l,2)
                 tt1=tt1+zcache(k+1)*fil0(l,2)
                 tt2=tt2+zcache(k+2)*fil0(l,2)
                 tt3=tt3+zcache(k+3)*fil0(l,2)
!!$                 tt4=tt4+zcache(k+4)*fil0(l,2)
!!$                 tt5=tt5+zcache(k+5)*fil0(l,2)
!!$                 tt6=tt6+zcache(k+6)*fil0(l,2)
!!$                 tt7=tt7+zcache(k+7)*fil0(l,2)

                 k=k+unrol
              end do
              k=k-(-lowfil+lupfil)*unrol
              y0(i1+j,i2,i3*unrol+0)=y0(i1+j,i2,i3*unrol+0)+tt0
              y0(i1+j,i2,i3*unrol+1)=y0(i1+j,i2,i3*unrol+1)+tt1
              y0(i1+j,i2,i3*unrol+2)=y0(i1+j,i2,i3*unrol+2)+tt2
              y0(i1+j,i2,i3*unrol+3)=y0(i1+j,i2,i3*unrol+3)+tt3
!!$              y0(i1+j,i2,i3*unrol+4)=y0(i1+j,i2,i3*unrol+4)+tt4
!!$              y0(i1+j,i2,i3*unrol+5)=y0(i1+j,i2,i3*unrol+5)+tt5
!!$              y0(i1+j,i2,i3*unrol+6)=y0(i1+j,i2,i3*unrol+6)+tt6
!!$              y0(i1+j,i2,i3*unrol+7)=y0(i1+j,i2,i3*unrol+7)+tt7
           end do
        end do
     end do


     do i1=(n1/nol)*nol,n1
        do i2=0,n2
           tt0=0.e0_wp
           tt1=0.e0_wp
           tt2=0.e0_wp
           tt3=0.e0_wp
!!$           tt4=0.e0_wp
!!$           tt5=0.e0_wp
!!$           tt6=0.e0_wp
!!$           tt7=0.e0_wp

           do l=lowfil,lupfil
              j=mod_arr(i2+l)

              tt0=tt0+x0(i1,j,i3*unrol+0)*fil0(l,2)
              tt1=tt1+x0(i1,j,i3*unrol+1)*fil0(l,2)
              tt2=tt2+x0(i1,j,i3*unrol+2)*fil0(l,2)
              tt3=tt3+x0(i1,j,i3*unrol+3)*fil0(l,2)
!!$              tt4=tt4+x0(i1,j,i3*unrol+4)*fil0(l,2)
!!$              tt5=tt5+x0(i1,j,i3*unrol+5)*fil0(l,2)
!!$              tt6=tt6+x0(i1,j,i3*unrol+6)*fil0(l,2)
!!$              tt7=tt7+x0(i1,j,i3*unrol+7)*fil0(l,2)
           end do
           y0(i1,i2,i3*unrol+0)=y0(i1,i2,i3*unrol+0)+tt0
           y0(i1,i2,i3*unrol+1)=y0(i1,i2,i3*unrol+1)+tt1
           y0(i1,i2,i3*unrol+2)=y0(i1,i2,i3*unrol+2)+tt2
           y0(i1,i2,i3*unrol+3)=y0(i1,i2,i3*unrol+3)+tt3
!!$           y0(i1,i2,i3*unrol+4)=y0(i1,i2,i3*unrol+4)+tt4
!!$           y0(i1,i2,i3*unrol+5)=y0(i1,i2,i3*unrol+5)+tt5
!!$           y0(i1,i2,i3*unrol+6)=y0(i1,i2,i3*unrol+6)+tt6
!!$           y0(i1,i2,i3*unrol+7)=y0(i1,i2,i3*unrol+7)+tt7
        end do
     end do

  end do
  !$omp end do

  !$omp do 
  do i3=(n3/unrol)*unrol,n3
     do i1=0,n1
        do i2=0,n2
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=mod_arr(i2+l)
              tt=tt+x0(i1,j   ,i3)*fil0(l,2)
           end do
           y0(i1,i2,i3)=y0(i1,i2,i3)+tt
        end do
     end do
  end do
  !$omp end do
END SUBROUTINE conv_kin_y_new

  subroutine conv_kin_z_new(mod_arr,fil0,x0,y0,n,ndat)
    implicit none
    integer, parameter :: wp=kind(1.0d0)
    integer, parameter :: lowfil=-14,lupfil=14,ncache=2*1024
    integer, parameter :: unrol=4,nol=unrol
    integer, intent(in) :: n,ndat
    integer, dimension(lowfil:n+lupfil), intent(in) :: mod_arr   
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
    real(wp),intent(in):: x0(ndat,0:n)
    real(wp),intent(inout)::y0(ndat,0:n)
    integer :: i3,i,l,j,k,ndim!,nol
    real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12
    real(wp), dimension(ncache) :: zcache

    !number of lines in ndat to be treated simultaneously 
    ndim=(n-lowfil+lupfil+1)
!!$    nol=min(ncache/(unrol*ndim),unrol)
!!$
!!$    !for the moment the number of lines equals the unrolling pattern
!!$    if (nol <unrol) stop 'nol'

    !$omp do 
    do i=0,ndat/unrol-1
       !copy buffer
       l=1
       do i3=1,-lowfil
          do k=1,nol
             zcache(l)=&
                  !zcache(i3+k*ndim)=&
                  x0(i*unrol+k,i3+lowfil-1+n+1)
             l=l+1
          end do
       end do
       do i3=-lowfil+1,n-lowfil+1
          do k=1,nol
             zcache(l)=&
                  x0(i*unrol+k,i3+lowfil-1)
             l=l+1
          end do
       end do
       do i3=n-lowfil+2,n-lowfil+1+lupfil
          do k=1,nol
             zcache(l)=&
                  x0(i*unrol+k,i3+lowfil-1-n-1)
             l=l+1
          end do
       end do

!zcache=zcache(unrol,ndim)

       k=0
       do i3=0,n
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
!!$          tt5=0.e0_wp
!!$          tt6=0.e0_wp
!!$          tt7=0.e0_wp
!!$          tt8=0.e0_wp
!!$          tt9 =0.e0_wp
!!$          tt10=0.e0_wp
!!$          tt11=0.e0_wp
!!$          tt12=0.e0_wp

          do l=lowfil,lupfil

             tt1=tt1+  zcache(k+1 )*fil0(l,3)
             tt2=tt2+  zcache(k+2 )*fil0(l,3)
             tt3=tt3+  zcache(k+3 )*fil0(l,3)
             tt4=tt4+  zcache(k+4 )*fil0(l,3)
!!$             tt5=tt5+  zcache(k+5 )*fil0(l,3)
!!$             tt6=tt6+  zcache(k+6 )*fil0(l,3)
!!$             tt7=tt7+  zcache(k+7 )*fil0(l,3)
!!$             tt8=tt8+  zcache(k+8 )*fil0(l,3)
!!$             tt9 =tt9 +zcache(k+9 )*fil0(l,3)
!!$             tt10=tt10+zcache(k+10)*fil0(l,3)
!!$             tt11=tt11+zcache(k+11)*fil0(l,3)
!!$             tt12=tt12+zcache(k+12)*fil0(l,3)

             k=k+unrol

          end do
          k=k-(-lowfil+lupfil)*unrol

          y0(i*unrol+1,i3)=y0(i*unrol+1,i3)+tt1
          y0(i*unrol+2,i3)=y0(i*unrol+2,i3)+tt2
          y0(i*unrol+3,i3)=y0(i*unrol+3,i3)+tt3
          y0(i*unrol+4,i3)=y0(i*unrol+4,i3)+tt4
!!$          y0(i*unrol+5,i3)=y0(i*unrol+5,i3)+tt5
!!$          y0(i*unrol+6,i3)=y0(i*unrol+6,i3)+tt6
!!$          y0(i*unrol+7,i3)=y0(i*unrol+7,i3)+tt7
!!$          y0(i*unrol+8,i3)=y0(i*unrol+8,i3)+tt8
!!$          y0(i*unrol+9 ,i3)=y0(i*unrol+9 ,i3)+tt9 
!!$          y0(i*unrol+10,i3)=y0(i*unrol+10,i3)+tt10
!!$          y0(i*unrol+11,i3)=y0(i*unrol+11,i3)+tt11
!!$          y0(i*unrol+12,i3)=y0(i*unrol+12,i3)+tt12
       end do
    end do
    !$omp end do

    !$omp do 
    do i=(ndat/unrol)*unrol+1,ndat
       do i3=0,n
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr(i3+l)
             tt=tt+x0(i,j)*fil0(l,3)
          end do
          y0(i,i3)=y0(i,i3)+tt
       end do
    end do
    !$omp end do
  END SUBROUTINE conv_kin_z_new


  subroutine conv_kin_y(mod_arr,fil0,x0,y0,n1,n2,n3)
    implicit none
    !Arguments
    integer, intent(in) :: n1,n2,n3
    integer, dimension(lowfil:n2+lupfil), intent(in) :: mod_arr
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
    real(wp), dimension(0:n1,0:n2,0:n3), intent(in) :: x0
    real(wp), dimension(0:n1,0:n2,0:n3), intent(inout) :: y0
    integer :: i1,i2,i3,l,j
    real(wp) :: tt,tt0,tt1,tt2,tt3,tt4,tt5,tt6,tt7

    !$omp do
    do i3=0,n3/8-1
        do i1=0,n1
           do i2=0,n2
              tt0=0.e0_wp
              tt1=0.e0_wp
              tt2=0.e0_wp
              tt3=0.e0_wp
              tt4=0.e0_wp
              tt5=0.e0_wp
              tt6=0.e0_wp
              tt7=0.e0_wp
    
              do l=lowfil,lupfil
                 j=mod_arr(i2+l)
    
                 tt0=tt0+x0(i1,j,i3*8+0)*fil0(l,2)
                 tt1=tt1+x0(i1,j,i3*8+1)*fil0(l,2)
                 tt2=tt2+x0(i1,j,i3*8+2)*fil0(l,2)
                 tt3=tt3+x0(i1,j,i3*8+3)*fil0(l,2)
                 tt4=tt4+x0(i1,j,i3*8+4)*fil0(l,2)
                 tt5=tt5+x0(i1,j,i3*8+5)*fil0(l,2)
                 tt6=tt6+x0(i1,j,i3*8+6)*fil0(l,2)
                 tt7=tt7+x0(i1,j,i3*8+7)*fil0(l,2)
              end do
              y0(i1,i2,i3*8+0)=y0(i1,i2,i3*8+0)+tt0
              y0(i1,i2,i3*8+1)=y0(i1,i2,i3*8+1)+tt1
              y0(i1,i2,i3*8+2)=y0(i1,i2,i3*8+2)+tt2
              y0(i1,i2,i3*8+3)=y0(i1,i2,i3*8+3)+tt3
              y0(i1,i2,i3*8+4)=y0(i1,i2,i3*8+4)+tt4
              y0(i1,i2,i3*8+5)=y0(i1,i2,i3*8+5)+tt5
              y0(i1,i2,i3*8+6)=y0(i1,i2,i3*8+6)+tt6
              y0(i1,i2,i3*8+7)=y0(i1,i2,i3*8+7)+tt7
           end do
        end do
     end do
     !$omp end do

     !$omp do 
     do i3=(n3/8)*8,n3
        do i1=0,n1
           do i2=0,n2
              tt=0.e0_wp
              do l=lowfil,lupfil
                 j=mod_arr(i2+l)
                 tt=tt+x0(i1,j   ,i3)*fil0(l,2)
              end do
              y0(i1,i2,i3)=y0(i1,i2,i3)+tt
           end do
        end do
     end do
     !$omp end do
  END SUBROUTINE conv_kin_y


  subroutine conv_kin_z(mod_arr,fil0,x0,y0,n,ndat)
    implicit none
    integer, intent(in) :: n,ndat
    integer, dimension(lowfil:n+lupfil), intent(in) :: mod_arr   
    real(wp), dimension(lowfil:lupfil,3), intent(in) :: fil0
    real(wp),intent(in):: x0(ndat,0:n)
    real(wp),intent(inout)::y0(ndat,0:n)
    integer :: i3,i,l,j
    real(wp) :: tt,tt1,tt2,tt3,tt4,tt5,tt6,tt7,tt8,tt9,tt10,tt11,tt12

    !$omp do 
    do i=0,ndat/12-1
       do i3=0,n3
          tt1=0.e0_wp
          tt2=0.e0_wp
          tt3=0.e0_wp
          tt4=0.e0_wp
          tt5=0.e0_wp
          tt6=0.e0_wp
          tt7=0.e0_wp
          tt8=0.e0_wp
          tt9 =0.e0_wp
          tt10=0.e0_wp
          tt11=0.e0_wp
          tt12=0.e0_wp

          do l=lowfil,lupfil
             j=mod_arr(i3+l)

             tt1=tt1+x0(i*12+1,j)*fil0(l,3)
             tt2=tt2+x0(i*12+2,j)*fil0(l,3)
             tt3=tt3+x0(i*12+3,j)*fil0(l,3)
             tt4=tt4+x0(i*12+4,j)*fil0(l,3)
             tt5=tt5+x0(i*12+5,j)*fil0(l,3)
             tt6=tt6+x0(i*12+6,j)*fil0(l,3)
             tt7=tt7+x0(i*12+7,j)*fil0(l,3)
             tt8=tt8+x0(i*12+8,j)*fil0(l,3)
             tt9 =tt9 +x0(i*12+9 ,j)*fil0(l,3)
             tt10=tt10+x0(i*12+10,j)*fil0(l,3)
             tt11=tt11+x0(i*12+11,j)*fil0(l,3)
             tt12=tt12+x0(i*12+12,j)*fil0(l,3)
          end do

          y0(i*12+1,i3)=y0(i*12+1,i3)+tt1
          y0(i*12+2,i3)=y0(i*12+2,i3)+tt2
          y0(i*12+3,i3)=y0(i*12+3,i3)+tt3
          y0(i*12+4,i3)=y0(i*12+4,i3)+tt4
          y0(i*12+5,i3)=y0(i*12+5,i3)+tt5
          y0(i*12+6,i3)=y0(i*12+6,i3)+tt6
          y0(i*12+7,i3)=y0(i*12+7,i3)+tt7
          y0(i*12+8,i3)=y0(i*12+8,i3)+tt8
          y0(i*12+9 ,i3)=y0(i*12+9 ,i3)+tt9 
          y0(i*12+10,i3)=y0(i*12+10,i3)+tt10
          y0(i*12+11,i3)=y0(i*12+11,i3)+tt11
          y0(i*12+12,i3)=y0(i*12+12,i3)+tt12
       end do
    end do
    !$omp end do

    !$omp do 
    do i=(ndat/12)*12+1,ndat
       do i3=0,n3
          tt=0.e0_wp
          do l=lowfil,lupfil
             j=mod_arr(i3+l)
             tt=tt+x0(i,j)*fil0(l,3)
          end do
          y0(i,i3)=y0(i,i3)+tt
       end do
    end do
    !$omp end do
  END SUBROUTINE conv_kin_z


END SUBROUTINE kinetic_per_ref_optim

subroutine fill_mod_arr(arr,nleft,nright,n)
  implicit none
  integer,intent(in) :: nleft,nright,n
  integer,intent(out) :: arr(nleft:nright)
  integer :: i
  
  if (nleft >= -n) then
     do i=nleft,-1
        arr(i)=n+i
     end do
  else
     do i=nleft,-1
        arr(i)=modulo(i,n)
     end do
  endif
  
  do i=max(0,nleft),min(n-1,nright)
     arr(i)=i
  end do
  
  if (nright < 2*n) then
     do i=n,nright
        arr(i)=i-n
     end do
  else
     do i=n,nright
        arr(i)=modulo(i,n)
     end do
  endif
END SUBROUTINE fill_mod_arr

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def kinetic_per_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "kinetic_per_ref"
  n1 = BOAST::Int("n1", :dir => :in)
  n2 = BOAST::Int("n2", :dir => :in)
  n3 = BOAST::Int("n3", :dir => :in)
  hgrid = BOAST::Real("hgrid", :dir => :in, :dim => [BOAST::Dim(3)] )
  c = BOAST::Real("c", :dir => :in)
  x = BOAST::Real("x", :dir => :in,  :dim => [ BOAST::Dim(0, n1 - 1),  BOAST::Dim(0, n2 - 1), BOAST::Dim(0, n3 - 1)] )
  y = BOAST::Real("y", :dir => :out, :dim => [ BOAST::Dim(0, n1 - 1),  BOAST::Dim(0, n2 - 1), BOAST::Dim(0, n3 - 1)] )
  
  p = BOAST::Procedure::new(function_name, [n1,n2,n3,hgrid,x,y,c])
  kernel.code.print <<EOF
subroutine kinetic_per_ref(n1,n2,n3,hgrid,x,y,c)
!   applies the kinetic energy operator onto x to get y. Works for periodic BC
  implicit none
  integer, parameter :: wp=kind(1.0d0), gp=kind(1.0d0)
  integer, intent(in) :: n1,n2,n3
  real(gp),intent(in)::c
  real(gp), dimension(3), intent(in) :: hgrid
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(in) :: x
  real(wp), dimension(0:n1-1,0:n2-1,0:n3-1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-14,lupfil=14
  integer :: i1,i2,i3,i,l,j
  real(wp) :: tt
  real(wp), dimension(3) :: scale
  real(wp), dimension(lowfil:lupfil,3) :: fil
  
  scale(:)=real(-.5_gp/hgrid(:)**2,wp)

  ! second derivative filters for Daubechies 16
  fil(0,:)=   -3.5536922899131901941296809374e0_wp*scale(:)
  fil(1,:)=    2.2191465938911163898794546405e0_wp*scale(:)
  fil(2,:)=   -0.6156141465570069496314853949e0_wp*scale(:)
  fil(3,:)=    0.2371780582153805636239247476e0_wp*scale(:)
  fil(4,:)=   -0.0822663999742123340987663521e0_wp*scale(:)
  fil(5,:)=    0.02207029188482255523789911295638968409e0_wp*scale(:)
  fil(6,:)=   -0.409765689342633823899327051188315485e-2_wp*scale(:)
  fil(7,:)=    0.45167920287502235349480037639758496e-3_wp*scale(:)
  fil(8,:)=   -0.2398228524507599670405555359023135e-4_wp*scale(:)
  fil(9,:)=    2.0904234952920365957922889447361e-6_wp*scale(:)
  fil(10,:)=  -3.7230763047369275848791496973044e-7_wp*scale(:)
  fil(11,:)=  -1.05857055496741470373494132287e-8_wp*scale(:)
  fil(12,:)=  -5.813879830282540547959250667e-11_wp*scale(:)
  fil(13,:)=   2.70800493626319438269856689037647576e-13_wp*scale(:)
  fil(14,:)=  -6.924474940639200152025730585882e-18_wp*scale(:)

  do i=1,14
     fil(-i,:)=fil(i,:)
  enddo
  
  do i3=0,n3-1
     ! (1/2) d^2/dx^2
     do i2=0,n2-1
        do i1=0,n1-1
           tt=x(i1,i2,i3)*c
           do l=lowfil,lupfil
              j=modulo(i1+l,n1)
              tt=tt+x(j   ,i2,i3)*fil(l,1)
           enddo
           y(i1,i2,i3)=tt
        enddo
     enddo
     
     ! + (1/2) d^2/dy^2
     do i1=0,n1-1
        do i2=0,n2-1
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i2+l,n2)
              tt=tt+x(i1,j   ,i3)*fil(l,2)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
     
  enddo

  ! + (1/2) d^2/dz^2
  do i2=0,n2-1
     do i1=0,n1-1
        do i3=0,n3-1
           tt=0.e0_wp
           do l=lowfil,lupfil
              j=modulo(i3+l,n3)
              tt=tt+x(i1,i2,   j)*fil(l,3)
           enddo
           y(i1,i2,i3)=y(i1,i2,i3)+tt
        enddo
     enddo
  enddo
  
END SUBROUTINE kinetic_per_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end

def Kinetic(conv_filter, optims = GenericOptimization::new, ekin = false)


  kinetic_operation = GenericConvolutionOperator::new(conv_filter,  :transpose => 0, :a_x =>(not ekin), :a_y => true, :eks => ekin, :a => true, :work => true, :kinetic => true, :zero_out_work => true)
  kinetic_operation.optimize(optims)

  p, subops= kinetic_operation.procedure()
  
  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op| 
    BOAST::pr op 
    puts "chosen:"+ op.name
  }
  BOAST::pr p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| kinetic_operation.cost(*args) }
  return kernel

end

