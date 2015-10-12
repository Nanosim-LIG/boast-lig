require 'rubygems'
require 'BOAST'
require 'narray'
require "./GenericConvolution.rb" 

def generate_filter
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "generate_filter"
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :inout, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  p = BOAST::Procedure::new(function_name, [nord, c1D])

kernel.code.print <<EOF
subroutine generate_filter (nord, c1D)
  integer, intent(in) :: nord
  integer :: m
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2),&
     intent(inout) :: c1D
!> @file
!! @author
!!    Copyright (C) 2014 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


m=nord/2
  !Only nord=2,4,6,8,16
   ! Coefficients from Matematica program.
select case(nord)
case(2)

     c1D(-m:m,0) = (/&
     -1.d0/2.d0,0.d0,1.d0/2.d0&
     /)

     c1D(-m:m,-1) = (/&
     -3.d0/2.d0, 2.d0,-1.d0/2.d0&
     /)
     c1D(-m:m,1) = (/&
     1.d0/2.d0,-2.d0,3.d0/2.d0&
     /)

   case(4)

     c1D(-m:m,0) = (/&
     1.d0/12.d0,-2.d0/3.d0,0.d0,2.d0/3.d0,-1.d0/12.d0&
     /)
     c1D(-m:m,-1) = (/&
     -1.d0/4.d0,-5.d0/6.d0,3.d0/2.d0,-1.d0/2.d0,1.d0/12.d0&
     /)
     c1D(-m:m,-2) = (/&
     -25.d0/12.0,4.d0,-3.d0,4.d0/3.d0,-1.d0/4.d0&
     /)
     c1D(-m:m,1) = (/&
     -1.d0/12.d0,+1.d0/2.d0,-3.d0/2.d0,+5.d0/6.d0,+1.d0/4.d0&
     /)
     c1D(-m:m,2) = (/&
     1.d0/4.d0,-4.d0/3.d0,3.d0,-4.d0,25.d0/12.d0&
     /)

   case(6)

     c1D(-m:m,0) = (/&
     -1.d0/60.d0,3.d0/20.d0,-3.d0/4.d0,0.d0,&
3.d0/4.d0,-3.d0/20.d0,1.d0/60.d0&
     /)
     c1D(-m:m,-1) = (/&
     1.d0/30.d0,-2/5.d0,-7.d0/12.d0,4.d0/3.d0,&
-1.d0/2.d0,2.d0/15.d0,-1.d0/60.d0&
     /)
     c1D(-m:m,1) = (/&
     1.d0/60.d0,-2.d0/15.d0,1.d0/2.d0,-4.d0/3.d0,&
7.d0/12.d0,2.d0/5.d0,-1.d0/30.d0&
     /)
     c1D(-m:m,-2) = (/&
     -1.d0/6.d0,-77.d0/60.d0,5.d0/2.d0,-5.d0/3.d0,&
5.d0/6.d0,-1.d0/4.d0,1.d0/30.d0&
     /)
     c1D(-m:m,2) = (/&
     -1.d0/30.d0,1.d0/4.d0,-5/6.d0,5.d0/3.d0,&
-5.d0/2.d0,77.d0/60.d0,1.d0/6.d0&
     /)
     c1D(-m:m,-3) = (/&
     -49.d0/20.d0,6.d0,-15.d0/2.d0,20.d0/3.d0,&
-15.d0/4.d0,6.d0/5.d0,-1.d0/6.d0&
     /)
     c1D(-m:m,3) = (/&
     1.d0/6.d0,-6.d0/5.d0,15.d0/4.d0,-20.d0/3.d0,&
15.d0/2.d0,-6.d0,49.d0/20.d0&
     /)

   case(8)

     c1D(-m:m,0) = (/&
     1.d0/280.d0,-4.d0/105.d0,1.d0/5.d0,-4.d0/5.d0,0.d0,&
4.d0/5.d0,-1.d0/5.d0,4.d0/105.d0,-1.d0/280.d0&
     /)
     c1D(-m:m,-1) = (/&
     -1.d0/168.d0,1.d0/14.d0,-1.d0/2.d0,-9.d0/20.d0,5.d0/4.d0,&
-1.d0/2.d0,1.d0/6.d0,-1.d0/28.d0,1.d0/280.d0&
     /)
     c1D(-m:m,1) = (/&
     -1.d0/280.d0,1.d0/28.d0,-1.d0/6.d0,1.d0/2.d0,-5.d0/4.d0,&
9.d0/20.d0,1.d0/2.d0,-1.d0/14.d0,1.d0/168.d0&
     /)
     c1D(-m:m,-2) = (/&
     1.d0/56.d0,-2.d0/7.d0,-19.d0/20.d0,2.d0,-5.d0/4.d0,&
2.d0/3.d0,-1.d0/4.d0,2.d0/35.d0,-1.d0/168.d0&
     /)
     c1D(-m:m,2) = (/&
      1.d0/168.d0,-2.d0/35.d0,1.d0/4.d0,-2.d0/3.d0,&
5.d0/4.d0,-2.d0,19.d0/20.d0,2.d0/7.d0,-1.d0/56.d0&
     /)
     c1D(-m:m,-3) = (/&
     -1.d0/8.d0,-223.d0/140.d0,7.d0/2.d0,-7.d0/2.d0,&
35.d0/12.d0,-7.d0/4.d0,7.d0/10.d0,& 
     -1.d0/6.d0,1.d0/56.d0&
     /)
     c1D(-m:m,3) = (/&
     -1.d0/56.d0,1.d0/6.d0,-7.d0/10.d0,7.d0/4.d0,&
-35.d0/12.d0,7.d0/2.d0,-7.d0/2.d0,&
     223.d0/140.d0,1.d0/8.d0&
     /)
     c1D(-m:m,-4) = (/&
     -761.d0/280.d0,8.d0,-14.d0,56.d0/3.d0,&
35.d0/2.d0,56.d0/5.d0,-14.d0/3.d0,& 
     8.d0/7.d0,-1.d0/8.d0&
     /)
     c1D(-m:m,4) = (/&
     1.d0/8.d0,-8.d0/7.d0,14.d0/3.d0,-56.d0/5.d0,&
35.d0/2.d0,-56.d0/3.d0,14.d0,-8.d0,& 
     761.d0/280.d0&
     /)

   case(16)

     c1D(-m:m,0) = (/&
     1.d0/102960.d0,-8.d0/45045.d0,2.d0/1287.d0,-56.d0/6435.d0,7.d0/198.d0,&
     -56.d0/495.d0,14.d0/45.d0,&
     -8.d0/9.d0,0.d0,8.d0/9.d0,-14.d0/45.d0,56.d0/495.d0,-7.d0/198.d0,&
     56.d0/6435.d0,-2.d0/1287.d0,8.d0/45045.d0,-1.d0/102960.d0&
     /)
     c1D(-m:m,-1) = (/&
     -1.d0/80080.d0,1.d0/4290.d0,-3.d0/1430.d0,7.d0/572.d0,-7.d0/132.d0,&
     21.d0/110.d0,-7.d0/10.d0,-17.d0/72.d0,&
     9.d0/8.d0,-1.d0/2.d0,7.d0/30.d0,-21.d0/220.d0,7.d0/220.d0,&
     -7.d0/858.d0,3.d0/2002.d0,-1.d0/5720.d0,1.d0/102960.d0&
     /)
     c1D(-m:m,1) = (/&
     -1.d0/102960.d0,1.d0/5720.d0,-3.d0/2002.d0,7.d0/858.d0,-7.d0/220.d0,&
     21.d0/220.d0,-7.d0/30.d0,1.d0/2.d0,&
     -9.d0/8.d0,17.d0/72.d0,7.d0/10.d0,-21.d0/110.d0,&
     7.d0/132.d0,-7.d0/572.d0,3.d0/1430.d0,-1.d0/4290.d0,1.d0/80080.d0&
     /)
     c1D(-m:m,-2) = (/&
      1.d0/48048.d0,-2.d0/5005.d0,15.d0/4004.d0,-10.d0/429.d0,5.d0/44.d0,&
      -6.d0/11.d0,-1207.d0/2520.d0,10.d0/7.d0,&
      -45.d0/56.d0,10.d0/21.d0,-1.d0/4.d0,6.d0/55.d0,-5.d0/132.d0,&
      10.d0/1001.d0,-15.d0/8008.d0,2.d0/9009.d0,-1.d0/80080.d0&
      /)
     c1D(-m:m,2) = (/&
     1.d0/80080.d0,-2.d0/9009.d0,15.d0/8008.d0,-10.d0/1001.d0,5.d0/132.d0,&
     -6.d0/55.d0,1.d0/4.d0,-10.d0/21.d0,45.d0/56.d0,&      
     -10.d0/7.d0,1207.d0/2520.d0,6.d0/11.d0,-5.d0/44.d0,10.d0/429.d0,&
     -15.d0/4004.d0,2.d0/5005.d0,-1.d0/48048.d0&
     /)
     c1D(-m:m,-3) = (/&
     -1.d0/21840.d0,1.d0/1092.d0,-5.d0/546.d0,5.d0/78.d0,-5.d0/12.d0,&
     -20417.d0/27720.d0,11.d0/6.d0,-55.d0/42.d0,55.d0/56.d0,&
     -55.d0/84.d0,11.d0/30.d0,-1.d0/6.d0,5.d0/84.d0,-5.d0/312.d0,&
     5.d0/1638.d0,-1.d0/2730.d0,1.d0/48048.d0&
     /)
     c1D(-m:m,3) = (/&
     -1.d0/48048.d0,1.d0/2730.d0,-5.d0/1638.d0,5.d0/312.d0,-5.d0/84.d0,&
      1.d0/6.d0,-11.d0/30.d0,55.d0/84.d0,&
     -55.d0/56.d0,55.d0/42.d0,-11.d0/6.d0,20417.d0/27720.d0,5.d0/12.d0,&
     -5.d0/78.d0,5.d0/546.d0,-1.d0/1092.d0,1.d0/21840.d0&
     /)
     c1D(-m:m,-4) = (/&
     1.d0/7280.d0,-4.d0/1365.d0,3.d0/91.d0,-4.d0/13.d0,&
     -28271.d0/27720.d0,12.d0/5.d0,-11.d0/5.d0,+44.d0/21.d0,&
     -99.d0/56.d0,44.d0/35.d0,-11.d0/15.d0,12.d0/35.d0,-1.d0/8.d0,&
     4.d0/117.d0,-3.d0/455.d0,4.d0/5005.d0,-1.d0/21840.d0&
     /)
     c1D(-m:m,4) = (/&
     1.d0/21840.d0,-4.d0/5005.d0,3.d0/455.d0,-4.d0/117.d0,1.d0/8.d0,&
     -12.d0/35.d0,11.d0/15.d0,-44.d0/35.d0,&
     99.d0/56.d0,-44.d0/21.d0,11.d0/5.d0,-12.d0/5.d0,&
     28271.d0/27720.d0,4.d0/13.d0,-3.d0/91.d0,4.d0/1365.d0,-1.d0/7280.d0&
     /)
     c1D(-m:m,-5) = (/&
     -1.d0/1680.d0,1.d0/70.d0,-3.d0/14.d0,-485333.d0/360360.d0,13.d0/4.d0,&
     -39.d0/10.d0,143.d0/30.d0,-143.d0/28.d0,&
     1287.d0/280.d0,-143.d0/42.d0,143.d0/70.d0,-39.d0/40.d0,13.d0/36.d0,&
     -1.d0/10.d0,3.d0/154.d0,-1.d0/420.d0,1.d0/7280.d0&
     /)
    c1D(-m:m,5) = (/&
    -1.d0/7280.d0,1.d0/420.d0,-3.d0/154.d0,1.d0/10.d0,-13.d0/36.d0,&
     39.d0/40.d0,-143.d0/70.d0,143.d0/42.d0,&
    -1287.d0/280.d0,143.d0/28.d0,-143.d0/30.d0,39.d0/10.d0,-13.d0/4.d0,&
     485333.d0/360360.d0,3.d0/14.d0,-1.d0/70.d0,1.d0/1680.d0&
    /)
    c1D(-m:m,-6) = (/&
    1.d0/240.d0,-2.d0/15.d0,-631193.d0/360360.d0,14.d0/3.d0,-91.d0/12.d0,&
    182.d0/15.d0,-1001.d0/60.d0,286.d0/15.d0,&
    -143.d0/8.d0,286.d0/21.d0,-1001.d0/120.d0,182.d0/45.d0,&
    -91.d0/60.d0,14.d0/33.d0,-1.d0/12.d0,2.d0/195.d0,-1.d0/1680.d0&
    /)
    c1D(-m:m,6) = (/&
    1.d0/1680.d0,-2.d0/195.d0,1.d0/12.d0,-14.d0/33.d0,91.d0/60.d0,&
    -182.d0/45.d0,1001.d0/120.d0,-286.d0/21.d0,&
    143.d0/8.d0,-286.d0/15.d0,1001.d0/60.d0,-182.d0/15.d0,91.d0/12.d0,&
    -14.d0/3.d0,631193.d0/360360.d0,2.d0/15.d0,-1.d0/240.d0&
    /)
    c1D(-m:m,-7) = (/&
    -1.d0/16.d0,-835397.d0/360360.d0,15.d0/2.d0,-35.d0/2.d0,455.d0/12.d0,&
    -273.d0/4.d0,1001.d0/10.d0,-715.d0/6.d0,&
    6435.d0/56.d0,-715.d0/8.d0,1001.d0/18.d0,-273.d0/10.d0,455.d0/44.d0,&
    -35.d0/12.d0,15.d0/26.d0,-1.d0/14.d0,1.d0/240.d0&
    /)
    c1D(-m:m,7) = (/&
    -1.d0/240.d0,1.d0/14.d0,-15.d0/26.d0,35.d0/12.d0,-455.d0/44.d0,&
    273.d0/10.d0,-1001.d0/18.d0,715.d0/8.d0,&
    -6435.d0/56.d0,715.d0/6.d0,-1001.d0/10.d0,273.d0/4.d0,-455.d0/12.d0,&
    35.d0/2.d0,-15.d0/2.d0,835397.d0/360360.d0,1.d0/16.d0&
    /)
    c1D(-m:m,-8) = (/&
    -2436559.d0/720720.d0,16.d0,-60.d0,560.d0/3.d0,-455.d0,4368.d0/5.d0,&
    -4004.d0/3.d0,11440.d0/7.d0,&
    -6435.d0/4.d0,11440.d0/9.d0,-4004.d0/5.d0,4368.d0/11.d0,-455.d0/3.d0,&
    560.d0/13.d0,-60.d0/7.d0,16.d0/15.d0,-1.d0/16.d0&
    /)
    c1D(-m:m,8) = (/&
    1.d0/16.d0,-16.d0/15.d0,60.d0/7.d0,-560.d0/13.d0,455.d0/3.d0,&
    -4368.d0/11.d0,4004.d0/5.d0,-11440.d0/9.d0,&
    6435.d0/4.d0,-11440.d0/7.d0,4004.d0/3.d0,-4368.d0/5.d0,455.d0,&
    -560.d0/3.d0,60.d0,-16.d0,2436559.d0/720720.d0&
    /)

   end select
end subroutine generate_filter

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def div_u_i_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "div_u_i_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03), BOAST::Dim(3)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine div_u_i_ref(geocode,n01,n02,n03,u,du,nord,hgrids,c1D)
  implicit none

!  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,nord,geocode
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03,3) :: u
  real(kind=8), dimension(n01,n02,n03) :: du

  !c..local variables
  integer :: n,m,n_cell,ii
  integer :: i,j,ib,i1,i2,i3
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2), intent(in)  :: c1D
  real(kind=8) :: hx,hy,hz, d1, d2, d3
  logical :: perx,pery,perz

  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
!  perx=(geocode /= 'F')
!  pery=(geocode == 'P')
!  perz=(geocode /= 'F')

  perx=(geocode /= 1)
  pery=(geocode == 2)
  perz=(geocode /= 1)

  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if

  !Only nord=2,4,6,8,16

  select case(nord)
  case(2,4,6,8,16)
     !O.K.
  case default
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end select

  !call set_diff_corff (nord, m, c1D)
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01

           du(i1,i2,i3) = 0.0d0

           d1 = 0.d0
           if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,1)!/hx
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,1)!/hx
               end do
              end if
           else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,1)!/hx
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)!/hx
               end do
              end if
           else
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,1)!/hx
              end do
           end if
           d1=d1/hx

           d2 = 0.d0
           if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,2)!/hy
               end do
              else
               do j=-m,m
                 d2 = d2 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,2)!/hy
               end do
              end if
           else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                d2 = d2 + c1D(j,0)*u(i1,ii,i3,2)!/hy
               end do
              else
               do j=-m,m
                 d2 = d2 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)!/hy
               end do
              end if
           else
              do j=-m,m
                 d2 = d2 + c1D(j,0)*u(i1,i2 + j,i3,2)!/hy
              end do
           end if
           d2=d2/hy

           d3 = 0.d0
           if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,3)!/hz
               end do
              else
               do j=-m,m
                 d3 = d3 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,3)!/hz
               end do
              end if
           else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                d3 = d3 + c1D(j,0)*u(i1,i2,ii,3)!/hz
               end do
              else
               do j=-m,m
                d3 = d3 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)!/hz
               end do
              end if
           else
              do j=-m,m
                 d3 = d3 + c1D(j,0)*u(i1,i2,i3 + j,3)!/hz
              end do
           end if
           d3=d3/hz

           du(i1,i2,i3) = d1+d2+d3

        end do
     end do
  end do

end subroutine div_u_i_ref

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def nabla_u_and_square_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "nabla_u_and_square_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03), BOAST::Dim(3)] )
  du2 = BOAST::Real("du2", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,du2,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine nabla_u_and_square_ref(geocode,n01,n02,n03,u,du,du2,nord,hgrids,c1D)
  implicit none
  !c..declare the pass
  integer, intent(in) :: geocode 
  integer, intent(in) :: n01,n02,n03,nord
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03) :: u
  real(kind=8), dimension(n01,n02,n03,3) :: du
  real(kind=8), dimension(n01,n02,n03) :: du2

  !c..local variables
  integer :: n,m,n_cell
  integer :: i,j,ib,i1,i2,i3,ii
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
  real(kind=8) :: hx,hy,hz
  logical :: perx,pery,perz

  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions

  perx=(geocode /= 1)
  pery=(geocode == 2)
  perz=(geocode /= 1)


  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     call f_err_throw('Ngrid in has to be setted > than n=nord + 1')
     !stop
  end if

  !Only nord=2,4,6,8,16

  select case(nord)
  case(2,4,6,8,16)
     !O.K.
  case default
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end select

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01

           du(i1,i2,i3,1) = 0.0d0
           du2(i1,i2,i3) = 0.0d0

           if (i1.le.m) then
            if (perx) then
             do j=-m,m
              ii=modulo(i1 + j + n01 - 1, n01 ) + 1
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(ii,i2,i3)/hx
             end do
            else
             do j=-m,m
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
             end do
            end if
           else if (i1.gt.n01-m) then
            if (perx) then
             do j=-m,m
              ii=modulo(i1 + j - 1, n01 ) + 1
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(ii,i2,i3)/hx
             end do
            else
             do j=-m,m
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
             end do
            end if
           else
              do j=-m,m
                 du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
              end do
           end if

           du2(i1,i2,i3) = du(i1,i2,i3,1)*du(i1,i2,i3,1)
           du(i1,i2,i3,2) = 0.0d0

           if (i2.le.m) then
            if (pery) then
             do j=-m,m
              ii=modulo(i2 + j + n02 - 1, n02 ) + 1
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,ii,i3)/hy
             end do
            else
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do
            end if
           else if (i2.gt.n02-m) then
            if (pery) then
             do j=-m,m
              ii=modulo(i2 + j - 1, n02 ) + 1
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,ii,i3)/hy
             end do
            else
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
             end do
            end if
           else
              do j=-m,m
                 du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
              end do
           end if

           du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,2)*du(i1,i2,i3,2)

           du(i1,i2,i3,3) = 0.0d0

           if (i3.le.m) then
            if (perz) then
             do j=-m,m
              ii=modulo(i3 + j + n03 - 1, n03 ) + 1
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,ii)/hz
             end do
            else
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
            end if
           else if (i3.gt.n03-m) then
            if (perz) then
             do j=-m,m
              ii=modulo(i3 + j - 1, n03 ) + 1
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,ii)/hz
             end do
            else
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
             end do
            end if
           else
              do j=-m,m
                 du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
              end do
           end if

           du2(i1,i2,i3) = du2(i1,i2,i3) + du(i1,i2,i3,3)*du(i1,i2,i3,3)

        end do
     end do
  end do
end subroutine nabla_u_and_square_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def update_rhopol_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "update_rhopol_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
dlogeps = BOAST::Real("dlogeps", :dir => :in, :dim => [BOAST::Dim(3), BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)] )
  rhopol = BOAST::Real("rhopol", :dir => :inout, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)] )
  rhores2 = BOAST::Real("rhores2", :dir => :out)
  eta = BOAST::Real("eta", :dir => :in)
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,nord,hgrids,eta,dlogeps,rhopol,rhores2,c1D])
  kernel.code.print <<EOF
subroutine update_rhopol_ref(geocode,n01,n02,n03,u,nord,hgrids,eta,dlogeps,rhopol,rhores2,c1D)

  implicit none

  integer, intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,nord
  real(kind=8), intent(in) :: eta
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03), intent(in) :: u
  real(kind=8), dimension(3,n01,n02,n03), intent(in) :: dlogeps
  real(kind=8), dimension(n01,n02,n03), intent(inout) :: rhopol
  real(kind=8), intent(out) :: rhores2

  !c..local variables
  integer :: n,m,n_cell
  integer :: i,j,ib,i1,i2,i3,isp,i1_max,i2_max,ii
  !real(kind=8), parameter :: oneo4pi=0.25d0/pi_param
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
  real(kind=8) :: hx,hy,hz,max_diff,fact,dx,dy,dz,res,rho
  real(kind=8) :: oneo4pi
  logical :: perx,pery,perz

  oneo4pi=0.0079577471545947673d0
  !1.0d0/(16.d0*atan(1.d0))
  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 1)
  pery=(geocode == 2)
  perz=(geocode /= 1)

  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if


  !Only nord=2,4,6,8,16
  if (all(nord /=[2,4,6,8,16])) then
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end if

  rhores2=0.d0
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01

           dx=0.d0

           if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                dx = dx + c1D(j,0)*u(ii,i2,i3)
               end do
              else
               do j=-m,m
                dx = dx + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
               end do
              end if
           else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                dx = dx + c1D(j,0)*u(ii,i2,i3)
               end do
              else
               do j=-m,m
                dx = dx + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
               end do
              end if
           else
              do j=-m,m
                 dx = dx + c1D(j,0)*u(i1 + j,i2,i3)
              end do
           end if
           dx=dx/hx

           dy = 0.0d0
           if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                dy = dy + c1D(j,0)*u(i1,ii,i3)
               end do
              else
               do j=-m,m
                dy = dy + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
               end do
              end if
           else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                dy = dy + c1D(j,0)*u(i1,ii,i3)
               end do
              else
               do j=-m,m
                dy = dy + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
               end do
              end if
           else
              do j=-m,m
                 dy = dy + c1D(j,0)*u(i1,i2 + j,i3)
              end do
           end if
           dy=dy/hy

           dz = 0.0d0
           if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                dz = dz + c1D(j,0)*u(i1,i2,ii)
               end do
              else
               do j=-m,m
                dz = dz + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
               end do
              end if
           else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                dz = dz + c1D(j,0)*u(i1,i2,ii)
               end do
              else
               do j=-m,m
                dz = dz + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
               end do
              end if
           else
              do j=-m,m
                 dz = dz + c1D(j,0)*u(i1,i2,i3 + j)
              end do
           end if
           dz=dz/hz

           !retrieve the previous treatment
           res = dlogeps(1,i1,i2,i3)*dx + &
                dlogeps(2,i1,i2,i3)*dy + dlogeps(3,i1,i2,i3)*dz
           res = res*oneo4pi
           rho=rhopol(i1,i2,i3)
           res=res-rho
           res=eta*res
           rhores2=rhores2+res*res
           rhopol(i1,i2,i3)=res+rho

        end do
     end do
  end do

end subroutine update_rhopol_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end

def nabla_u_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "nabla_u_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),BOAST::Dim(3)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine nabla_u_ref(geocode,n01,n02,n03,u,du,nord,hgrids,c1D)
  implicit none

  !c..this routine computes 'nord' order accurate first derivatives 

  !c..input:
  !c..ngrid       = number of points in the grid, 
  !c..u(ngrid)    = function values at the grid points

  !c..output:
  !c..du(ngrid)   = first derivative values at the grid points

  !c..declare the pass
  integer, intent(in) :: geocode 
  integer, intent(in) :: n01,n02,n03,nord
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03) :: u
  real(kind=8), dimension(n01,n02,n03,3) :: du

  !c..local variables
  integer :: n,m,n_cell,ii
  integer :: i,j,ib,i1,i2,i3
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D
  real(kind=8) :: hx,hy,hz
  logical :: perx,pery,perz

  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
  perx=(geocode /= 1)
  pery=(geocode == 2)
  perz=(geocode /= 1)


  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if

  !Only nord=2,4,6,8,16

  select case(nord)
  case(2,4,6,8,16)
     !O.K.
  case default
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end select

  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           du(i1,i2,i3,1) = 0.0d0

           if (i1.le.m) then
            if (perx) then
             do j=-m,m
              ii=modulo(i1 + j + n01 - 1, n01 ) + 1
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(ii,i2,i3)/hx
             end do
            else
             do j=-m,m
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-m-1)*u(j+m+1,i2,i3)/hx
             end do
            end if
           else if (i1.gt.n01-m) then
            if (perx) then
             do j=-m,m
              ii=modulo(i1 + j - 1, n01 ) + 1
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(ii,i2,i3)/hx
             end do
            else
             do j=-m,m
              du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)/hx
             end do
            end if
           else
            do j=-m,m
             du(i1,i2,i3,1) = du(i1,i2,i3,1) + c1D(j,0)*u(i1 + j,i2,i3)/hx
            end do
           end if

           du(i1,i2,i3,2) = 0.0d0

           if (i2.le.m) then
            if (pery) then
             do j=-m,m
              ii=modulo(i2 + j + n02 - 1, n02 ) + 1
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,ii,i3)/hy
             end do
            else
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-m-1)*u(i1,j+m+1,i3)/hy
             end do
            end if
           else if (i2.gt.n02-m) then
            if (pery) then
             do j=-m,m
              ii=modulo(i2 + j - 1, n02 ) + 1
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,ii,i3)/hy
             end do
            else
             do j=-m,m
              du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)/hy
             end do
            end if
           else
            do j=-m,m
             du(i1,i2,i3,2) = du(i1,i2,i3,2) + c1D(j,0)*u(i1,i2 + j,i3)/hy
            end do
           end if

           du(i1,i2,i3,3) = 0.0d0

           if (i3.le.m) then
            if (perz) then
             do j=-m,m
              ii=modulo(i3 + j + n03 - 1, n03 ) + 1
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,ii)/hz
             end do
            else
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-m-1)*u(i1,i2,j+m+1)/hz
             end do
            end if
           else if (i3.gt.n03-m) then
            if (perz) then
             do j=-m,m
              ii=modulo(i3 + j - 1, n03 ) + 1
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,ii)/hz
             end do
            else
             do j=-m,m
              du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)/hz
             end do
            end if
           else
            do j=-m,m
             du(i1,i2,i3,3) = du(i1,i2,i3,3) + c1D(j,0)*u(i1,i2,i3 + j)/hz
            end do
           end if

        end do
     end do
  end do
end subroutine nabla_u_ref
EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end



#version without all the *** divisions, separating the 3 convolutions, and with openmp
def div_u_i_opt
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "div_u_i_opt"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),BOAST::Dim(3)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine div_u_i_opt(geocode,n01,n02,n03,u,du,nord,hgrids,c1D)
  implicit none
!  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,nord,geocode
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03,3) :: u
  real(kind=8), dimension(n01,n02,n03) :: du

  !c..local variables
  integer :: n,m,n_cell,ii,iii
  integer :: i,j,ib,i1,i2,i3
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2), intent(in)  :: c1D
  real(kind=8) :: hx,hy,hz, d1
  logical :: perx,pery,perz

  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
!  perx=(geocode /= 'F')
!  pery=(geocode == 'P')
!  perz=(geocode /= 'F')

  perx=(geocode /= 1)
  pery=(geocode == 2)
  perz=(geocode /= 1)

  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if

  !Only nord=2,4,6,8,16

  select case(nord)
  case(2,4,6,8,16)
     !O.K.
  case default
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end select

  !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d1) 
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01

           du(i1,i2,i3) = 0.0d0

           d1 = 0.d0
           if (i1.le.m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,1)
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i1-m-1)*u(j+m+1,i2,i3,1)
               end do
              end if
           else if (i1.gt.n01-m) then
              if (perx) then
               do j=-m,m
                ii=modulo(i1 + j - 1, n01 ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,1)
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)
               end do
              end if
           else
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1 + j,i2,i3,1)
              end do
           end if

        du(i1,i2,i3) =d1/hx
        end do
     end do
  end do
  !$omp end parallel do

  !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d1) 
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           d1=0.d0
           if (i2.le.m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                d1 = d1 + c1D(j,0)*u(i1,ii,i3,2)
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i2-m-1)*u(i1,j+m+1,i3,2)
               end do
              end if
           else if (i2.gt.n02-m) then
              if (pery) then
               do j=-m,m
                ii=modulo(i2 + j - 1, n02 ) + 1
                d1 = d1 + c1D(j,0)*u(i1,ii,i3,2)
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)
               end do
              end if
           else
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1,i2 + j,i3,2)
              end do
           end if
          du(i1,i2,i3) = du(i1,i2,i3) + d1/hy
        end do
     end do
  end do
  !$omp end parallel do

  !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d1) 
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           d1=0.d0
           if (i3.le.m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                d1 = d1 + c1D(j,0)*u(i1,i2,ii,3)
               end do
              else
               do j=-m,m
                 d1 = d1 + c1D(j,i3-m-1)*u(i1,i2,j+m+1,3)
               end do
              end if
           else if (i3.gt.n03-m) then
              if (perz) then
               do j=-m,m
                ii=modulo(i3 + j - 1, n03 ) + 1
                d1 = d1 + c1D(j,0)*u(i1,i2,ii,3)
               end do
              else
               do j=-m,m
                d1 = d1 + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)
               end do
              end if
           else
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1,i2,i3 + j,3)
              end do
           end if
           du(i1,i2,i3) = du(i1,i2,i3)+d1/hz
        end do
     end do
  end do
  !$omp end parallel do


end subroutine div_u_i_opt

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


#ugly one, just to check it works this way, and help debug

def div_u_i_1d
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "div_u_i_1d"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),BOAST::Dim(3)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  dim = BOAST::Int("dim", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,dim,nord,hgrids,c1D])

  kernel.code.print <<EOF
subroutine div_u_i_1d(geocode,n01,n02,n03,u,du,dim,nord,hgrids,c1D)
  implicit none
!  character(len=1), intent(in) :: geocode
  integer, intent(in) :: n01,n02,n03,nord,geocode
  real(kind=8), dimension(3), intent(in) :: hgrids
  real(kind=8), dimension(n01,n02,n03,3) :: u
  real(kind=8), dimension(n01,n02,n03) :: du

  !c..local variables
  integer :: n,m,n_cell,ii,iii
  integer :: i,j,i1,i2,i3,dim, it, dimit
  real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2), intent(in)  :: c1D
  real(kind=8) :: hx,hy,hz, d1,h
  logical :: perx

  n = nord+1
  m = nord/2
  hx = hgrids(1)!acell/real(n01,kind=8)
  hy = hgrids(2)!acell/real(n02,kind=8)
  hz = hgrids(3)!acell/real(n03,kind=8)
  n_cell = max(n01,n02,n03)

  !buffers associated to the geocode
  !conditions for periodicity in the three directions
!  perx=(geocode /= 'F')
!  pery=(geocode == 'P')
!  perz=(geocode /= 'F')

  ! Beware that n_cell has to be > than n.
  if (n_cell.lt.n) then
     write(*,*)'ngrid in has to be setted > than n=nord + 1'
     stop
  end if

  !Only nord=2,4,6,8,16

  select case(nord)
  case(2,4,6,8,16)
     !O.K.
  case default
     write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
     stop
  end select


  select case(dim)
  case(1)
    h=hx
    perx=(geocode /= 1)
  case(2)
    h=hy
    perx=(geocode == 2)
  case(3)
    h=hz
    perx=(geocode /= 1)
  case default
    write(*,*)'Only dim 1,2,3 allowed'
  stop
  end select

  !$omp parallel do default(shared) private(i1,i2,i3,j,ii,it,dimit, d1) 
  do i3=1,n03
     do i2=1,n02
        do i1=1,n01
           select case(dim)
           case(1)
             it=i1
             dimit=n01
           case(2)
             it=i2
             dimit=n02
           case(3)
             it=i3
             dimit=n03
           case default
             write(*,*)'Only dim 1,2,3 allowed'
           stop
           end select

           d1 = 0.d0
           if (it.le.m) then
              if (perx) then
                select case(dim)
                case(1)
               do j=-m,m
                ii=modulo(it + j + dimit - 1, dimit ) + 1
                    d1 = d1 + c1D(j,0)*u(ii,i2,i3,dim)
               end do
                case(2)
               do j=-m,m
                ii=modulo(it + j + dimit - 1, dimit ) + 1
                    d1 = d1 + c1D(j,0)*u(i1,ii,i3,dim)
               end do
                case(3)
               do j=-m,m
                ii=modulo(it + j + dimit - 1, dimit ) + 1
                    d1 = d1 + c1D(j,0)*u(i1,i2,ii,dim)
               end do
                case default
                    write(*,*)'Only dim 1,2,3 allowed'
                stop
                end select

              else

                select case(dim)
                case(1)
               do j=-m,m
                 d1 = d1 + c1D(j,it-m-1)*u(j+m+1,i2,i3,dim)
               end do
                case(2)
               do j=-m,m
                 d1 = d1 + c1D(j,it-m-1)*u(i1,j+m+1,i3,dim)
               end do
                case(3)
               do j=-m,m
                 d1 = d1 + c1D(j,it-m-1)*u(i1,i2,j+m+1,dim)
               end do
                case default
                    write(*,*)'Only dim 1,2,3 allowed'
                stop
                end select
    
              end if
           else if (it.gt.dimit-m) then
              if (perx) then

                select case(dim)
                case(1)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                d1 = d1 + c1D(j,0)*u(ii,i2,i3,dim)
               end do
                case(2)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                d1 = d1 + c1D(j,0)*u(i1,ii,i3,dim)
               end do
                case(3)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                d1 = d1 + c1D(j,0)*u(i1,i2,ii,dim)
               end do
                case default
                    write(*,*)'Only dim 1,2,3 allowed'
                stop
                end select

              else

                select case(dim)
                case(1)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                 d1 = d1 + c1D(j,it-dimit+m)*u(dimit + j - m,i2,i3,dim)
               end do
                case(2)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                 d1 = d1 + c1D(j,it-dimit+m)*u(i1,dimit + j - m,i3,dim)
               end do
                case(3)
               do j=-m,m
                ii=modulo(it + j - 1, dimit ) + 1
                 d1 = d1 + c1D(j,it-dimit+m)*u(i1,i2,dimit + j - m,dim)
               end do
                case default
                    write(*,*)'Only dim 1,2,3 allowed'
                stop
                end select


              end if
           else

                select case(dim)
                case(1)
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(it + j,i2,i3,dim)
              end do
                case(2)
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1,it + j,i3,dim)
              end do
                case(3)
              do j=-m,m
                 d1 = d1 + c1D(j,0)*u(i1,i2,it + j,dim)
              end do
                case default
                    write(*,*)'Only dim 1,2,3 allowed'
                stop
                end select


           end if

        du(i1,i2,i3) = du(i1,i2,i3) + d1/h
        end do
     end do
  end do
  !$omp end parallel do


end subroutine div_u_i_1d

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


def nabla_u_epsilon_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "nabla_u_epsilon_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03), BOAST::Dim(3)] )
  eps = BOAST::Real("eps", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,eps,c1D])
  kernel.code.print <<EOF

    subroutine nabla_u_epsilon_ref(geocode,n01,n02,n03,u,du,nord,hgrids,eps,c1D)
      implicit none

      integer, intent(in) :: geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u 
      real(kind=8), dimension(n01,n02,n03,3), intent(out) :: du 
      real(kind=8), dimension(n01,n02,n03), intent(in) :: eps 

      !c..local variables
      integer :: n,m,n_cell,ii
      integer :: i,j,ib,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz, d,e
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      !perx=(geocode /= 'F')
      !pery=(geocode == 'P')
      !perz=(geocode /= 'F')
      perx=(geocode /= 1)
      pery=(geocode == 2)
      perz=(geocode /= 1)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         write(*,*)'ngrid in has to be setted > than n=nord + 1'
         stop
      end if

      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select


      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e) 
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d= 0.0d0
               e=eps(i1,i2,i3)


               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               du(i1,i2,i3,1) = d*e
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e) 
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0
               e=eps(i1,i2,i3)
               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               du(i1,i2,i3,2)=d*e
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e) 
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0
               e=eps(i1,i2,i3)
               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               du(i1,i2,i3,3)=d*e
            end do
         end do
      end do
      !$omp end parallel do

    end subroutine nabla_u_epsilon_ref

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end



def nabla_u_square_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "nabla_u_square_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  du2 = BOAST::Real("du2", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [BOAST::Dim(3)] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du2,nord,hgrids,c1D])
  kernel.code.print <<EOF

    subroutine nabla_u_square_ref(geocode,n01,n02,n03,u,du2,nord,hgrids,c1D)
      implicit none


      integer, intent(in) :: geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u
      real(kind=8), dimension(n01,n02,n03), intent(out) :: du2 

      !c..local variables
      integer :: n,m,n_cell
      integer :: i,j,ib,i1,i2,i3,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz,d
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 1)
      pery=(geocode == 2)
      perz=(geocode /= 1)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         write(*,*)'ngrid in has to be setted > than n=nord + 1'
         stop
      end if

      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select


      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz

      !$omp parallel do default(none) &
      !$omp private(i3,i2,i1,d,ii,j) &
      !$omp shared(du2,perx,m,n01,n02,n03,c1D_1,u)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0
               !du2(i1,i2,i3) = 0.0d0

               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               !du(i1,i2,i3,1)= d
               du2(i1,i2,i3) = d*d

            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0

               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               !du(i1,i2,i3,2)= d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d) 
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0

               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               !du(i1,i2,i3,3)=d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d

            end do
         end do
      end do
      !$omp end parallel do
    end subroutine nabla_u_square_ref


EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


#  def div_u_i_3d(kernel_1D)
#    kernel = BOAST::CKernel::new(:kernels => [kernel_1D])
#    function_name = "div_u_i_3d"
#    n01 = BOAST::Int("n01", :dir => :in)
#    n02 = BOAST::Int("n02", :dir => :in)
#    n03 = BOAST::Int("n03", :dir => :in)
#    dim = BOAST::Int("dim", :dir => :in)
##    input = Variable::new("input",Real,{:direction => :in, :dimension => [ Dimension::new(n1 * n2 * n3) ] })
##    output = Variable::new("output",Real,{:direction => :out, :dimension => [ Dimension::new(n1 * n2 * n3) ] })
##    temp = Variable::new("temp",Real,{:direction => :out, :dimension => [ Dimension::new(n1 * n2 * n3) ] })

#  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),3] )
#  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
#  nord = BOAST::Int("nord", :dir => :in)
#  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [3] )

#    p = BOAST::Procedure::new(function_name,[n01,n02,n03,u,du, dim,nord, hgrids]) {
#      kernel_1D.procedure.call(n01,n02,n03, u, du, 1 , nord, hgrids).print
#      kernel_1D.procedure.call(n01,n02,n03, u, du, 2 , nord, hgrids).print
#      kernel_1D.procedure.call(n01,n02,n03, u, du, 3 , nord, hgrids).print
#    }
#    kernel.procedure = p
#    return kernel
#  end


def Poisson_conv(conv_filter, optims=GenericOptimization::new, accum=1)
    
if (accum !=0) then
  conv_operation = GenericConvolutionOperator1d::new(conv_filter, :transpose => 0, :work => false, :a_y => accum, :ld => false, :narr => false, :a => true, :poisson => true)
else
  conv_operation = GenericConvolutionOperator1d::new(conv_filter, :transpose => 0, :work => false, :ld => false, :narr => false, :a => true, :poisson => true)
end 


  #test of 1d kernels optimizations in view of many-d
  conv_operation.optimize(optims)

  p, subops= conv_operation.procedure()

  kernel = BOAST::CKernel::new

  print_header

  subops.each_value { |op| 
    BOAST::pr op 
    puts "chosen:"+ op.name
  }
  BOAST::pr p

  kernel.procedure = p
  kernel.cost_function = lambda { |*args| conv_operation.cost(*args) }
  return kernel

end


def Poisson_brocker(optims,n1,n2,n3)

    function_name = "Poisson_brocker"
#    nords = BOAST::Int("nords", :dim => [BOAST::Dim(5)], :local =>true)
    kernels=[]
    j = BOAST::Int "j"
    nord = BOAST::Int("nord", :dir => :in)
    idim = BOAST::Int("idim", :dir => :in)
    a = BOAST::Real("a", :dir => :in)
    bc = BOAST::Real("bc", :dir => :in)
    u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1), BOAST::Dim(0, 2)] )
    du=BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
    nn = BOAST::Int("nn", :dir => :in, :dim => [BOAST::Dim(3)])

    suffix = ".c" if BOAST::get_lang == BOAST::C
    suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
    BOAST::decl j
    nords=[2,4,6,8,16]
    
    nords.each{ |nord_n|
        filter= NArray.float(nord_n+1, nord_n+1)
        generate_filter.run(nord_n, filter)
        conv_filter = PoissonFilter::new('poisson'+nord_n.to_s,filter.to_a,nord_n)
        k = Poisson_conv(conv_filter, optims)
        k.build(:openmp => true)
        kernels.push k
    }



    kernel = BOAST::CKernel::new(:kernels => [kernels[0],kernels[1],kernels[2],kernels[3],kernels[4]])
    p = BOAST::Procedure(function_name, [nord,idim, nn, bc, u, du, a]){
        BOAST::pr BOAST::Case( nord, 2, lambda {
          BOAST::pr  kernels[0].procedure.call(3,idim,nn, bc, u, du, a)
        },4,lambda {
          BOAST::pr  kernels[1].procedure.call(3,idim,nn, bc, u, du, a)
        },6,lambda {
          BOAST::pr  kernels[2].procedure.call(3,idim,nn, bc, u, du, a)
        },8,lambda {
          BOAST::pr  kernels[3].procedure.call(3,idim,nn, bc, u, du, a)
        },16,lambda {
          BOAST::pr  kernels[4].procedure.call(3,idim,nn, bc, u, du, a)
        })
    }

    BOAST::pr p
    kernel.procedure = p


    c = lambda { |*args| 3 }

#BOAST::Procedure(function_name, [nord,*args]){
#        BOAST::pr BOAST::Case( nord, 2, lambda {
#          kernels[0].cost(*args)
#        },4,lambda {
#          kernels[1].cost(*args)
#        },6,lambda {
#          kernels[2].cost(*args)
#        },8,lambda {
#          kernels[3].cost(*args)
#        },16,lambda {
#          kernels[4].cost(*args)
#        })
#    }

    kernel.cost_function = c

  File::open("poisson_kernels#{suffix}","w") { |f|
    f.puts kernels[0]
    f.puts kernels[1]
    f.puts kernels[2]
    f.puts kernels[3]
    f.puts kernels[4]
    f.puts kernel
  }

    return kernel

end 

def div_u_i(n1,n2,n3,k2)

  function_name = "div_u_i"
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1), BOAST::Dim(0, 2)] )
  du=BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)] )
  cc=BOAST::Real("cc", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :optional => true )
  hgrids = BOAST::Real("hgrids",:dir => :in, :dim => [ BOAST::Dim(0, 2)])
  nord = BOAST::Int("nord", :dir => :in)
  geocode = BOAST::Int("geocode", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03,u,du,nord,hgrids,cc]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    uxy=BOAST::Real("uxy", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :allocate => :heap )
    uxz=BOAST::Real("uxz", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :allocate => :heap )
    uyz=BOAST::Real("uyz", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :allocate => :heap )
    du1=BOAST::Real("du1", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :allocate => :heap )
    du2=BOAST::Real("du2", :dir => :out, :dim => [ BOAST::Dim(0, n01-1),  BOAST::Dim(0, n02-1), BOAST::Dim(0, n03-1)], :allocate => :heap )
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::decl uxy, uxz, uyz, du1, du2
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03

    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }
    
    BOAST::pr BOAST::If(present(cc), lambda{
        du1.alloc
        du2.alloc
        uxy.alloc
        uyz.alloc
        uxz.alloc
        BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u[0,0,0,0].address, du, a0)
        BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u[0,0,0,1].address, du1, a1)
        BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u[0,0,0,2].address, du2, a2)
        BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u[0,0,0,1].address, uxy, a0)
        BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u[0,0,0,2].address, uxz, a0)
        BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u[0,0,0,2].address, uyz, a1)

    BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: nil, private: [i1,i2,i3]) { 
        BOAST::pr BOAST::For(i3, 0,n3-1,openmp: true){
          BOAST::pr BOAST::For(i2, 0,n2-1){
            BOAST::pr BOAST::For(i1, 0,n1-1){
        cc[i1, i2, i3] === (u[i1,i2,i3,0]*u[i1,i2,i3,0])*du[i1,i2,i3] + 
                       2.d0*u[i1,i2,i3,0]*u[i1,i2,i3,1]*uxy[i1,i2,i3] + 
                       2.d0*u[i1,i2,i3,0]*u[i1,i2,i3,2]*uxz[i1,i2,i3] + 
                         (u[i1,i2,i3,1]*u[i1,i2,i3,1])*du1[i1,i2,i3]+ 
                       2.d0*u[i1,i2,i3,1]*u[i1,i2,i3,2]*uyz[i1,i2,i3]+
                       (u[i1,i2,i3,2]*u[i1,i2,i3,2])*du2[i1,i2,i3]
        du[i1,i2,i3]=du[i1,i2,i3]+du1[i1,i2,i3]+du2[i1,i2,i3]
            }
          }
        }
      }

    },lambda{
        BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u[0,0,0,0].address, du, a0)
        BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u[0,0,0,1].address, du, a1)
        BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u[0,0,0,2].address, du, a2)
    }
    }
    BOAST::pr p
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
    return kernel
end


def nabla_u_and_square(n1,n2,n3,k2)

  function_name = "nabla_u_and_square"
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1-1), BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1), BOAST::Dim(0, 2)] )
  ddu = BOAST::Real("ddu", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  hgrids = BOAST::Real("hgrids",:dir => :in, :dim => [ BOAST::Dim(0, 2)])
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)
  geocode = BOAST::Int("geocode", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03,u,du,ddu,nord,hgrids]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    i1 = BOAST::Int("i1")
    i2 = BOAST::Int("i2")
    i3 = BOAST::Int("i3")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::decl i1
    BOAST::decl i2
    BOAST::decl i3
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03
    BOAST::pr i1 === 0
    BOAST::pr i2 === 0
    BOAST::pr i3 === 0
    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }

    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]
    
    BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u, du[0,0,0,0].address, a0)
    BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u, du[0,0,0,1].address, a1)
    BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u, du[0,0,0,2].address, a2)
    BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: nil, private: [i1,i2,i3]) { 
        BOAST::pr BOAST::For(i3, 0,n3-1,openmp: true){
          BOAST::pr BOAST::For(i2, 0,n2-1){
            BOAST::pr BOAST::For(i1, 0,n1-1){
               BOAST::pr ddu[i1,i2,i3] === (du[i1,i2,i3,0]*du[i1,i2,i3,0]) + (du[i1,i2,i3,1]*du[i1,i2,i3,1]) + (du[i1,i2,i3,2]*du[i1,i2,i3,2])
            }
          }
        }
      }
    }

    
    BOAST::pr p
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
    return kernel
end

def nabla_u_square(n1,n2,n3,k2)

  function_name = "nabla_u_square"
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  ddu = BOAST::Real("ddu", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  hgrids = BOAST::Real("hgrids",:dir => :in, :dim => [ BOAST::Dim(0, 2)])
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)
  geocode = BOAST::Int("geocode", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03,u,ddu,nord,hgrids]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    i1 = BOAST::Int("i1")
    i2 = BOAST::Int("i2")
    i3 = BOAST::Int("i3")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1-1), BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1), BOAST::Dim(0, 2)] , :local =>true, :allocate => :heap)
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::decl du
    BOAST::decl i1
    BOAST::decl i2
    BOAST::decl i3
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03
    BOAST::pr i1 === 0
    BOAST::pr i2 === 0
    BOAST::pr i3 === 0
    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }

    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]
    du.alloc
    BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u, du[0,0,0,0].address, a0)
    BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u, du[0,0,0,1].address, a1)
    BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u, du[0,0,0,2].address, a2)
    BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: nil, private: [i1,i2,i3]) { 
        BOAST::pr BOAST::For(i3, 0,n3-1,openmp: true){
          BOAST::pr BOAST::For(i2, 0,n2-1){
            BOAST::pr BOAST::For(i1, 0,n1-1){
               BOAST::pr ddu[i1,i2,i3] === (du[i1,i2,i3,0]*du[i1,i2,i3,0]) + (du[i1,i2,i3,1]*du[i1,i2,i3,1]) + (du[i1,i2,i3,2]*du[i1,i2,i3,2])
            }
          }
        }
      }
    du.dealloc
    }

    
    BOAST::pr p
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
    return kernel
end


def nabla_u_epsilon(n1,n2,n3,k2)

  function_name = "nabla_u_epsilon"
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1-1), BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1),BOAST::Dim(0, 2)] )
  eps = BOAST::Real("eps", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  hgrids = BOAST::Real("hgrids",:dir => :in, :dim => [ BOAST::Dim(0, 2)])
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)
  geocode = BOAST::Int("geocode", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03,u,du,nord,hgrids,eps]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    i1 = BOAST::Int("i1")
    i2 = BOAST::Int("i2")
    i3 = BOAST::Int("i3")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::decl i1
    BOAST::decl i2
    BOAST::decl i3
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03
    BOAST::pr i1 === 0
    BOAST::pr i2 === 0
    BOAST::pr i3 === 0
    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }

    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]
    
    BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u, du[0,0,0,0].address, a0)
    BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u, du[0,0,0,1].address, a1)
    BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u, du[0,0,0,2].address, a2)
    BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: nil, private: [i1,i2,i3]) { 
        BOAST::pr BOAST::For(i3, 0,n3-1,openmp: true){
          BOAST::pr BOAST::For(i2, 0,n2-1){
            BOAST::pr BOAST::For(i1, 0,n1-1){
              BOAST::pr du[i1,i2,i3,0] === du[i1,i2,i3,0]*eps[i1,i2,i3]
              BOAST::pr du[i1,i2,i3,1] === du[i1,i2,i3,1]*eps[i1,i2,i3]
              BOAST::pr du[i1,i2,i3,2] === du[i1,i2,i3,2]*eps[i1,i2,i3]
            }
          }
        }
      }
    }

    
    BOAST::pr p
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
    return kernel
end


def update_rhopol(n1,n2,n3,k2)

  function_name = "update_rhopol"
  u       = BOAST::Real("u",       :dir => :in,    :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  rhopol  = BOAST::Real("rhopol",  :dir => :inout, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  dlogeps = BOAST::Real("dlogeps", :dir => :in,    :dim => [ BOAST::Dim(3),        BOAST::Dim(0, n1-1), BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  hgrids  = BOAST::Real("hgrids",  :dir => :in,    :dim => [ BOAST::Dim(0, 2)])
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  nord = BOAST::Int("nord")
  rhores2 = BOAST::Real("rhores2", :dir => :out)
  geocode = BOAST::Int("geocode", :dir => :in)
  eta = BOAST::Real("eta", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03, u,nord, hgrids,eta,dlogeps,rhopol,rhores2]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    du      = BOAST::Real("du",      :dir => :out,   :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1), BOAST::Dim(0, 2)],  :local =>true, :allocate =>:heap )
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    i1 = BOAST::Int("i1")
    i2 = BOAST::Int("i2")
    i3 = BOAST::Int("i3")
    res = BOAST::Real("res")
    rho = BOAST::Real("rho")
    oneo4pi = BOAST::Real("oneo4pi")
    tmp_rhores2 = BOAST::Real("tmp_rhores2")
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::decl i1
    BOAST::decl i2
    BOAST::decl i3
    BOAST::decl res
    BOAST::decl rho
    BOAST::decl du
    BOAST::decl oneo4pi
    BOAST::decl tmp_rhores2
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03
    BOAST::pr i1 === 0
    BOAST::pr i2 === 0
    BOAST::pr i3 === 0
    BOAST::pr oneo4pi === 0.0079577471545947673


    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }
    
    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]
    
    du.alloc
    BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u, du[0,0,0,0].address, a0)
    BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u, du[0,0,0,1].address, a1)
    BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u, du[0,0,0,2].address, a2)
    BOAST::pr tmp_rhores2 === 0.0
    BOAST::pr BOAST::OpenMP::Parallel(default: :shared, reduction: {"+" => tmp_rhores2}, private: [i1,i2,i3,res,rho]) { 
      BOAST::pr BOAST::For(i3, 0,n3-1, openmp: true){
        BOAST::pr BOAST::For(i2, 0,n2-1){
          BOAST::pr BOAST::For(i1, 0,n1-1){
            BOAST::pr res === dlogeps[1,i1,i2,i3]*du[i1,i2,i3,0]+dlogeps[2,i1,i2,i3]*du[i1,i2,i3,1]+dlogeps[3,i1,i2,i3]*du[i1,i2,i3,2]
            BOAST::pr res === res*oneo4pi
            BOAST::pr rho === rhopol[i1,i2,i3]
            BOAST::pr res === (res-rho)*eta
            BOAST::pr tmp_rhores2 === tmp_rhores2 + res*res
            BOAST::pr rhopol[i1,i2,i3] === res+rho
          }
        }
      }
    }
    du.dealloc
    BOAST::pr rhores2 === tmp_rhores2
  }

    
  BOAST::pr p
  kernel.procedure = p
  kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
  return kernel
end



def nabla_u(n1,n2,n3,k2)

  function_name = "nabla_u"
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
    du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1),BOAST::Dim(0,2)] )
#  du0 = BOAST::Real("du0", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
#  du1 = BOAST::Real("du1", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
#  du2 = BOAST::Real("du2", :dir => :out, :dim => [ BOAST::Dim(0, n1-1),  BOAST::Dim(0, n2-1), BOAST::Dim(0, n3-1)] )
  hgrids = BOAST::Real("hgrids",:dir => :in, :dim => [ BOAST::Dim(0, 2)])
  geocode = BOAST::Int("geocode", :dir => :in)
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)

  kernel = BOAST::CKernel::new(:kernels => [k2])


#  suffix = ".c" if BOAST::get_lang == BOAST::C
#  suffix = ".f90" if BOAST::get_lang == BOAST::FORTRAN
#  File::open("poisson_kernels#{suffix}","w") { |f|
#    f.puts k2
#  }

#  print_header


  p = BOAST::Procedure::new(function_name,[geocode, n01,n02,n03,u,du,nord,hgrids]){
    nn = BOAST::Int("nn", :dim => [BOAST::Dim(3)], :local =>true)
    bc0 = BOAST::Int("bc0")
    bc1 = BOAST::Int("bc1")
    bc2 = BOAST::Int("bc2")
    a0 = BOAST::Real("a0")
    a1 = BOAST::Real("a1")
    a2 = BOAST::Real("a2")
    BOAST::decl nn, bc0, bc1, bc2, a0, a1, a2
    BOAST::pr nn[1] === n01
    BOAST::pr nn[2] === n02
    BOAST::pr nn[3] === n03
    BOAST::pr bc0===BC::PERIODIC
    BOAST::pr bc1===BC::PERIODIC
    BOAST::pr bc2===BC::PERIODIC

    BOAST::pr BOAST::If(geocode == 1) {
      BOAST::pr bc0 === BC::NPERIODIC
      BOAST::pr bc2 === BC::NPERIODIC
    }
    BOAST::pr BOAST::If(geocode != 2){
      BOAST::pr bc1 === BC::NPERIODIC
    }
    
    BOAST::pr a0 === BOAST::Real(1.0) / hgrids[0]
    BOAST::pr a1 === BOAST::Real(1.0) / hgrids[1]
    BOAST::pr a2 === BOAST::Real(1.0) / hgrids[2]
    
    BOAST::pr k2.procedure.call(nord, 0, nn, bc0, u, du[0,0,0,0].address, a0)
    BOAST::pr k2.procedure.call(nord, 1, nn, bc1, u, du[0,0,0,1].address, a1)
    BOAST::pr k2.procedure.call(nord, 2, nn, bc2, u, du[0,0,0,2].address, a2)

    }

    
    BOAST::pr p
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(nord,*args) }
    return kernel
end



