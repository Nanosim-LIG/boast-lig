!> Simple non-optimized version of the major convolution routines
subroutine magicfilter_free_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,-lupfil:n1-1-lowfil), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /

  do j=1,ndat
     do i=-lupfil,n1-1-lowfil
        tt=0.e0_wp
        do l=max(-i,lowfil),min(lupfil,n1-1-i)
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_free_ref


!> Simple non-optimized version of the major convolution routines
subroutine magicfilter_free_t_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, parameter :: lowfil=-8,lupfil=7
  integer, intent(in) :: n1,ndat
  real(wp), dimension(lowfil:n1-1+lupfil,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1-1), intent(out) :: y
  !local variables
  integer :: i,j,l
  real(wp) :: tt
  ! the filtered output data structure has grown by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /

  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           tt=tt+x(i+l,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_free_t_ref

subroutine magicfilter_per_ref(n1,ndat,x,y)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n1,ndat
  real(kind=8), dimension(0:(n1-1),ndat), intent(in) :: x
  real(kind=8), dimension(ndat,0:(n1-1)), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-8,lupfil=7
  integer :: i,j,k,l
  real(kind=8) :: tt

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(kind=8) fil(lowfil:lupfil)
  DATA fil / &
       8.4334247333529341094733325815816e-7_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.9940415697834003993178616713_wp,&
       -0.604895289196983516002834636e-1_wp, &
       -0.2103025160930381434955489412839065067e-1_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       2.72734492911979659657715313017228e-6_wp /
  do j=1,ndat
     do i=0,n1-1
        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)  
           tt=tt+x(  k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_per_ref

subroutine magicfilter_per_t_ref(n1,ndat,x,y)
  !use module_base
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  integer, intent(in) :: n1,ndat
  real(wp), dimension(0:n1-1,ndat), intent(in) :: x
  real(wp), dimension(ndat,0:n1-1), intent(out) :: y
  !local variables
  integer, parameter :: lowfil=-7,lupfil=8
  integer :: i,j,k,l
  real(wp) :: tt
  ! the filtered output data structure has shrunk by the filter length

  !          THE MAGIC FILTER FOR DAUBECHIES-16
  real(wp) fil(lowfil:lupfil)
  DATA fil / &
       2.72734492911979659657715313017228e-6_wp,&
       -0.5185986881173432922848639136911487e-4_wp,&
       0.49443227688689919192282259476750972e-3_wp,&
       -0.344128144493493857280881509686821861e-2_wp,&
       0.1337263414854794752733423467013220997e-1_wp,&
       -0.2103025160930381434955489412839065067e-1_wp,&
       -0.604895289196983516002834636e-1_wp,&
       0.9940415697834003993178616713_wp,&
       0.612625895831207982195380597e-1_wp,&
       0.2373821463724942397566389712597274535e-1_wp,&
       -0.942047030201080385922711540948195075e-2_wp,&
       0.174723713672993903449447812749852942e-2_wp,&
       -0.30158038132690463167163703826169879e-3_wp,&
       0.8762984476210559564689161894116397e-4_wp,&
       -0.1290557201342060969516786758559028e-4_wp,&
       8.4334247333529341094733325815816e-7_wp /


  do j=1,ndat
     do i=0,n1-1

        tt=0.e0_wp
        do l=lowfil,lupfil
           k=modulo(i+l,n1)
           tt=tt+x(k,j)*fil(l)
        enddo
        y(j,i)=tt
     enddo
  enddo

END SUBROUTINE magicfilter_per_t_ref
