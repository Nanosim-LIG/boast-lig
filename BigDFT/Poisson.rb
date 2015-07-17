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


def fssnord3dmatnabla3var_lg_ref
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "fssnord3dmatnabla3var_lg_ref"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),3] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [3] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine fssnord3dmatnabla3var_lg_ref(geocode,n01,n02,n03,u,du,nord,hgrids,c1D)
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

end subroutine fssnord3dmatnabla3var_lg_ref

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end

#version without all the *** divisions, separating the 3 convolutions, and with openmp
def fssnord3dmatnabla3var_lg_opt
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "fssnord3dmatnabla3var_lg_opt"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),3] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [3] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,nord,hgrids,c1D])
  kernel.code.print <<EOF
subroutine fssnord3dmatnabla3var_lg_opt(geocode,n01,n02,n03,u,du,nord,hgrids,c1D)
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


end subroutine fssnord3dmatnabla3var_lg_opt

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end


#ugly one, just to check it works this way, and help debug

def fssnord3dmatnabla3var_lg_1d
  lang = BOAST::get_lang
  BOAST::set_lang(BOAST::FORTRAN)
  kernel = BOAST::CKernel::new
  function_name = "fssnord3dmatnabla3var_lg_1d"
  geocode = BOAST::Int("geocode", :dir => :in )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03),3] )
  du = BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n01),  BOAST::Dim(0, n02), BOAST::Dim(0, n03)] )
  dim = BOAST::Int("dim", :dir => :in)
  nord = BOAST::Int("nord", :dir => :in)
  c1D = BOAST::Real("c1D", :dir => :in, :dim =>  [ BOAST::Dim(-nord/2,nord/2), BOAST::Dim(-nord/2,nord/2)])
  hgrids = BOAST::Real("hgrids", :dir => :in, :dim => [3] )
  p = BOAST::Procedure::new(function_name, [geocode,n01,n02,n03,u,du,dim,nord,hgrids,c1D])

  kernel.code.print <<EOF
subroutine fssnord3dmatnabla3var_lg_1d(geocode,n01,n02,n03,u,du,dim,nord,hgrids,c1D)
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


  select case(dim)
  case(1)
    h=hx
  case(2)
    h=hy
  case(3)
    h=hz
  case default
    write(*,*)'Only dim 1,2,3 allowed'
  stop
  end select

  !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d1) 
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


end subroutine fssnord3dmatnabla3var_lg_1d

EOF
  kernel.procedure = p
  BOAST::set_lang(lang)
  return kernel
end



#  def fssnord3dmatnabla3var_lg_3d(kernel_1D)
#    kernel = BOAST::CKernel::new(:kernels => [kernel_1D])
#    function_name = "fssnord3dmatnabla3var_lg_3d"
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


def Poisson_conv(conv_filter, optims=GenericOptimization::new)
    

  conv_operation = GenericConvolutionOperator1d::new(conv_filter, :transpose => 0, :work => false, :a_y => 1, :ld => false, :narr => false)

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

FILTER = [
   9.7125097125097125E-006,
  -1.7760017760017760E-004,
   1.5540015540015540E-003,
  -8.7024087024087024E-003,
  0.035353535353535348     ,
  -0.11313131313131313     ,
   0.31111111111111112     ,
  -0.88888888888888875     ,
   0.0000000000000000     ,
   0.88888888888888875     ,
  -0.31111111111111112     ,
   0.11313131313131313     ,
 -0.035353535353535348     ,
   8.7024087024087024E-003,
  -1.5540015540015540E-003,
   1.7760017760017760E-004,
  -9.7125097125097125E-006
]

def fssnord3dmatnabla3var_lg_full(n1,n2,n3,hgrids)

  function_name = "fssnord3dmatnabla3var_lg_full"
  u = BOAST::Real("u", :dir => :in, :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3),3] )
  du=BOAST::Real("du", :dir => :out, :dim => [ BOAST::Dim(0, n1),  BOAST::Dim(0, n2), BOAST::Dim(0, n3)] )
  n01 = BOAST::Int("n01", :dir => :in)
  n02 = BOAST::Int("n02", :dir => :in)
  n03 = BOAST::Int("n03", :dir => :in)
  n = NArray.int(3)
  n[0] = n1
  n[1] = n2
  n[2] = n3
  #for now let's say that all filters are the same
  filter_x=FILTER.collect { |n| n / hgrids[0] }
  conv_filter = ConvolutionFilter::new('poisson',filter_x,8)
  optims = GenericOptimization::new(:unroll_range => 1, :mod_arr_test => true,:tt_arr_test => true, :dimensions => [n1,n2,n3])

  k2 = Poisson_conv(conv_filter, optims)
  k2.build(:openmp => true)
  kernel = BOAST::CKernel::new(:kernels => [k2])

  print_header

  bc = NArray.int(3)
  bc[0] = BC::PERIODIC
  bc[1] = BC::PERIODIC
  bc[2] = BC::PERIODIC

  p = BOAST::Procedure::new(function_name,[n01,n02,n03,u,du]){
   BOAST::pr k2.procedure.call(3, 0, n, BC::PERIODIC, u, du, 0.5)
   BOAST::pr k2.procedure.call(3, 1, n, BC::PERIODIC, u[0..(n01-1), 0..(n02-1), 0..(n03-1), 1], du, 0.5)
   BOAST::pr k2.procedure.call(3, 2, n, BC::PERIODIC, u[0..(n01-1), 0..(n02-1), 0..(n03-1), 2], du, 0.5)
  }
    kernel.procedure = p
    kernel.cost_function = lambda { |*args| 3*k2.cost(*args)+2 }
    return kernel
end

