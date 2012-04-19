subroutine magicfilter_per_1(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
      enddo
      y(j + 0, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_1
subroutine magicfilter_per_2(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 1, 2
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
    enddo
  enddo
  do j=ndat - modulo(ndat, 2) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_2
subroutine magicfilter_per_3(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 2, 3
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
    enddo
  enddo
  do j=ndat - modulo(ndat, 3) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_3
subroutine magicfilter_per_4(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 3, 4
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
    enddo
  enddo
  do j=ndat - modulo(ndat, 4) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_4
subroutine magicfilter_per_5(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 4, 5
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
    enddo
  enddo
  do j=ndat - modulo(ndat, 5) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_5
subroutine magicfilter_per_6(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 5, 6
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
    enddo
  enddo
  do j=ndat - modulo(ndat, 6) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_6
subroutine magicfilter_per_7(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 6, 7
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
    enddo
  enddo
  do j=ndat - modulo(ndat, 7) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_7
subroutine magicfilter_per_8(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8) :: tt8
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 7, 8
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      tt8 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
        tt8 = tt8 + x(k, j + 7) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
      y(j + 7, i) = tt8
    enddo
  enddo
  do j=ndat - modulo(ndat, 8) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_8
subroutine magicfilter_per_9(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8) :: tt8
  real(kind=8) :: tt9
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 8, 9
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      tt8 = 0.0
      tt9 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
        tt8 = tt8 + x(k, j + 7) * fil(l)
        tt9 = tt9 + x(k, j + 8) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
      y(j + 7, i) = tt8
      y(j + 8, i) = tt9
    enddo
  enddo
  do j=ndat - modulo(ndat, 9) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_9
subroutine magicfilter_per_10(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8) :: tt8
  real(kind=8) :: tt9
  real(kind=8) :: tt10
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 9, 10
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      tt8 = 0.0
      tt9 = 0.0
      tt10 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
        tt8 = tt8 + x(k, j + 7) * fil(l)
        tt9 = tt9 + x(k, j + 8) * fil(l)
        tt10 = tt10 + x(k, j + 9) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
      y(j + 7, i) = tt8
      y(j + 8, i) = tt9
      y(j + 9, i) = tt10
    enddo
  enddo
  do j=ndat - modulo(ndat, 10) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_10
subroutine magicfilter_per_11(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8) :: tt8
  real(kind=8) :: tt9
  real(kind=8) :: tt10
  real(kind=8) :: tt11
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 10, 11
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      tt8 = 0.0
      tt9 = 0.0
      tt10 = 0.0
      tt11 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
        tt8 = tt8 + x(k, j + 7) * fil(l)
        tt9 = tt9 + x(k, j + 8) * fil(l)
        tt10 = tt10 + x(k, j + 9) * fil(l)
        tt11 = tt11 + x(k, j + 10) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
      y(j + 7, i) = tt8
      y(j + 8, i) = tt9
      y(j + 9, i) = tt10
      y(j + 10, i) = tt11
    enddo
  enddo
  do j=ndat - modulo(ndat, 11) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_11
subroutine magicfilter_per_12(n, ndat, x, y)
  integer(kind=4), parameter :: lowfil = -8
  integer(kind=4), parameter :: upfil = 7
  integer(kind=4), intent(in) :: n
  integer(kind=4), intent(in) :: ndat
  real(kind=8), intent(in), dimension(0:n-1, ndat) :: x
  real(kind=8), intent(out), dimension(ndat, 0:n-1) :: y
  integer(kind=4) :: i
  integer(kind=4) :: j
  integer(kind=4) :: k
  integer(kind=4) :: l
  real(kind=8) :: tt1
  real(kind=8) :: tt2
  real(kind=8) :: tt3
  real(kind=8) :: tt4
  real(kind=8) :: tt5
  real(kind=8) :: tt6
  real(kind=8) :: tt7
  real(kind=8) :: tt8
  real(kind=8) :: tt9
  real(kind=8) :: tt10
  real(kind=8) :: tt11
  real(kind=8) :: tt12
  real(kind=8), parameter, dimension(lowfil:upfil) :: fil = (/ &
8.4334247333529341094733325815816e-7, &
-0.1290557201342060969516786758559028e-4, &
0.8762984476210559564689161894116397e-4, &
-0.30158038132690463167163703826169879e-3, &
0.174723713672993903449447812749852942e-2, &
-0.942047030201080385922711540948195075e-2, &
0.2373821463724942397566389712597274535e-1, &
0.612625895831207982195380597e-1, &
0.9940415697834003993178616713, &
-0.604895289196983516002834636e-1, &
-0.2103025160930381434955489412839065067e-1, &
0.1337263414854794752733423467013220997e-1, &
-0.344128144493493857280881509686821861e-2, &
0.49443227688689919192282259476750972e-3, &
-0.5185986881173432922848639136911487e-4, &
2.72734492911979659657715313017228e-6 /)
  do j=1, ndat - 11, 12
    do i=0, n-1
      tt1 = 0.0
      tt2 = 0.0
      tt3 = 0.0
      tt4 = 0.0
      tt5 = 0.0
      tt6 = 0.0
      tt7 = 0.0
      tt8 = 0.0
      tt9 = 0.0
      tt10 = 0.0
      tt11 = 0.0
      tt12 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j + 0) * fil(l)
        tt2 = tt2 + x(k, j + 1) * fil(l)
        tt3 = tt3 + x(k, j + 2) * fil(l)
        tt4 = tt4 + x(k, j + 3) * fil(l)
        tt5 = tt5 + x(k, j + 4) * fil(l)
        tt6 = tt6 + x(k, j + 5) * fil(l)
        tt7 = tt7 + x(k, j + 6) * fil(l)
        tt8 = tt8 + x(k, j + 7) * fil(l)
        tt9 = tt9 + x(k, j + 8) * fil(l)
        tt10 = tt10 + x(k, j + 9) * fil(l)
        tt11 = tt11 + x(k, j + 10) * fil(l)
        tt12 = tt12 + x(k, j + 11) * fil(l)
      enddo
      y(j + 0, i) = tt1
      y(j + 1, i) = tt2
      y(j + 2, i) = tt3
      y(j + 3, i) = tt4
      y(j + 4, i) = tt5
      y(j + 5, i) = tt6
      y(j + 6, i) = tt7
      y(j + 7, i) = tt8
      y(j + 8, i) = tt9
      y(j + 9, i) = tt10
      y(j + 10, i) = tt11
      y(j + 11, i) = tt12
    enddo
  enddo
  do j=ndat - modulo(ndat, 12) + 1, ndat
    do i=0, n-1
      tt1 = 0.0
      do l=lowfil, upfil
        k = modulo(i + l, n)
        tt1 = tt1 + x(k, j) * fil(l)
      enddo
      y(j, i) = tt1
    enddo
  enddo
END SUBROUTINE magicfilter_per_12
