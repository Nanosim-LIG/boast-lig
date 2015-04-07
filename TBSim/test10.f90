program test10

  use defs
  use linalg

  interface ! Interface to the C nstime function (returns time in ns as integer*8).
    integer(kind = 8) function nstime()
    end function nstime
  end interface

! Parameters.

  integer, parameter :: ntests = 100000
  character(len = 1), parameter :: trans = "T"

! Global variables.

  integer :: i, j
  integer :: itest
  integer(kind = 8) :: tstart, tstop

  real(dp) :: a(10, 10)
  complex(dp) :: b(10), c(10), cref(10)

! Code.

  write (6, '("Trans = ",a1,".")') trans

! Initialize a & b.

  call random_number(a)

  b = (/(cmplx(.1d0, .2d0)*i, i = 1, 10)/)

! Reference time = matmul.

  tstart = nstime()

  if (trans == "N") then
    do itest = 1, ntests
      cref = matmul(a, b)
    end do
  else
    do itest = 1, ntests
      cref = matmul(transpose(a), b) ! Transpose on the fly.
    end do
  end if

  tstop = nstime()

! NB : 10 real*complex multiplications+9 complex additions = 38 flops per matrix line = 380 flops per operation.

  write (6, '("MATMUL                                    = ",i10," mus (",f6.3," Gflops).")') &
&   (tstop-tstart)/1000, 380*ntests/real(tstop-tstart, dp)

#ifdef MKL

! MKL DZGEMV.

! MKL might allocate workspace memory on first call, so we make a blank call before timing.

!  call DZGEMV(trans, 10, 10, (1.d0, 0.d0), a, 10, b, 1, (0.d0, 0.d0), c, 1)

!  tstart = nstime()

!  do itest = 1, ntests
!    call DZGEMV(trans, 10, 10, (1.d0, 0.d0), a, 10, b, 1, (0.d0, 0.d0), c, 1)
!  end do

!  tstop = nstime()

!  write (6, '("MKL DZGEMV                                = ",i10," mus (",f6.3," Gflops).")') &
!&   (tstop-tstart)/1000, 380*ntests/real(tstop-tstart, dp)

! MKL double DGEMV.

! MKL might allocate workspace memory on first call, so we make a blank call before timing.

  call DGEMV(trans, 10, 10, (1.d0, 0.d0), a, 10, b(1), 2, (0.d0, 0.d0), c(1), 2)

  tstart = nstime()

  do itest = 1, ntests
    call DGEMV(trans, 10, 10, (1.d0, 0.d0), a, 10, b(1), 2, (0.d0, 0.d0), c(1), 2)
    call DGEMV(trans, 10, 10, (1.d0, 0.d0), a, 10, b(2), 2, (0.d0, 0.d0), c(2), 2)
  end do

  tstop = nstime()

  write (6, '("MKL double DGEMV                          = ",i10," mus (",f6.3," Gflops).")') &
&   (tstop-tstart)/1000, 380*ntests/real(tstop-tstart, dp)

#endif

! dmatzvec_10x10 unrolled by the compiler.

  tstart = nstime()

  do itest = 1, ntests
    call dmatzvec_10x10(trans, 10, 10, a, 10, b, c)
  end do

  tstop = nstime()

  write (6, '("DMATZVEC_10X10 (unrolled by the compiler) = ",i10," mus (",f6.3," Gflops).")') &
&   (tstop-tstart)/1000, 380*ntests/real(tstop-tstart, dp)

! dmatzvec_10x10 unrolled by hand.

  tstart = nstime()

  do itest = 1, ntests
    call dmatzvec_10x10_unrolled(trans, 10, 10, a, 10, b, c)
  end do

  tstop = nstime()

  write (6, '("Norm of residual = ",d12.5,".")') sqrt(zsqnorm(c-cref, 10))

  write (6, '("DMATZVEC_10X10 (unrolled by hand)         = ",i10," mus (",f6.3," Gflops).")') &
&   (tstop-tstart)/1000, 380*ntests/real(tstop-tstart, dp)

contains

!
! (n, m) real matrix-complex vector product with n, m = 1 or 10.
!
! Unrolled by the compiler.
!

  subroutine dmatzvec_10x10(trans, n, m, a, lda, b, c)

    implicit none

! Input/output variables.

    integer, intent(in) :: lda
    integer, intent(in) :: n, m
    real(dp), intent(in) :: a(lda, *)
    complex(dp), intent(in) :: b(*)
    complex(dp), intent(out) :: c(*)
    character(len = 1), intent(in) :: trans

! Local variables.

    integer :: i, j
    complex(dp) :: bj

! Subroutine code.

    if (trans == "N") then

      if (n == 1) then

        c(1) = a( 1, 1)*b(1)

        do j = 2, m ! Do not unroll : Seldom used, and unfavorable stride anyway.
          c(1) = c(1)+a(1, j)*b(j)
        end do

      else if (n == 10) then

        bj = b(1)
        do i = 1, 10
          c(i) = a(i, 1)*bj
        end do

        do j = 2, m
          bj = b(j)
          do i = 1, 10
            c(i) = c(i)+a(i, j)*bj
          end do
        end do

      else

        stop "dmatzvec_10x10 : Error, invalid (n, m) for trans = 'N'."

      end if

    else if (trans == "T") then

      if (m == 1) then

        bj = b(1)
        do i = 1, n ! Do not unroll : Seldom used.
          c(i) = a(1, i)*bj
        end do

      else if (m == 10) then

        do i = 1, n
          c(i) = 0.d0
          do j = 1, 10
            c(i) = c(i)+a(j, i)*b(j)
          end do
        end do

      else

        stop "dmatzvec_10x10 : Error, invalid (n, m) for trans = 'T'."

      end if

    else

      stop "dmatzvec_10x10 : Error, invalid trans."

    end if

  end subroutine dmatzvec_10x10

!
! Complex scalar product with n = 1 or 10.
!
! Unrolled by the compiler.
!

  complex(dp) function zdotc_10(n, a, b)

    implicit none

! Input/output variables.

    integer, intent(in) :: n
    complex(dp), intent(in) :: a(*), b(*)

! Local variables.

    integer :: i

! Subroutine code.

    if (n == 1) then

      zdotc_10 = conjg(a(1))*b(1)

    else if (n == 10) then

      zdotc_10 = 0.d0
      do i = 1, 10
        zdotc_10 = zdotc_10+conjg(a(i))*b(i)
      end do

    else

      stop "zdotc_10 : Error, invalid n."

    end if

  end function zdotc_10

!
! (n, m) real matrix-complex vector product with n, m = 1 or 10.
!
! Unrolled by hand.
!

  subroutine dmatzvec_10x10_unrolled(trans, n, m, a, lda, b, c)

    implicit none

! Input/output variables.

    integer, intent(in) :: lda
    integer, intent(in) :: n, m
    real(dp), intent(in) :: a(lda, *)
    complex(dp), intent(in) :: b(*)
    complex(dp), intent(out) :: c(*)
    character(len = 1), intent(in) :: trans

! Local variables.

    integer :: i, j
    complex(dp) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, bj

! Subroutine code.

    if (trans == "N") then

      if (n == 1) then

        c(1) = a( 1, 1)*b(1)

        do j = 2, m ! Do not unroll : Seldom used, and unfavorable stride anyway.
          c(1) = c(1)+a(1, j)*b(j)
        end do

      else if (n == 10) then

!         bj = b(1)
!         c( 1) = a( 1, 1)*bj
!         c( 2) = a( 2, 1)*bj
!         c( 3) = a( 3, 1)*bj
!         c( 4) = a( 4, 1)*bj
!         c( 5) = a( 5, 1)*bj
!         c( 6) = a( 6, 1)*bj
!         c( 7) = a( 7, 1)*bj
!         c( 8) = a( 8, 1)*bj
!         c( 9) = a( 9, 1)*bj
!         c(10) = a(10, 1)*bj
!
!         do j = 2, m
!           bj = b(j)
!           c( 1) = c( 1)+a( 1, j)*bj
!           c( 2) = c( 2)+a( 2, j)*bj
!           c( 3) = c( 3)+a( 3, j)*bj
!           c( 4) = c( 4)+a( 4, j)*bj
!           c( 5) = c( 5)+a( 5, j)*bj
!           c( 6) = c( 6)+a( 6, j)*bj
!           c( 7) = c( 7)+a( 7, j)*bj
!           c( 8) = c( 8)+a( 8, j)*bj
!           c( 9) = c( 9)+a( 9, j)*bj
!           c(10) = c(10)+a(10, j)*bj
!         end do

#ifdef SSE
        call dmatzvec_n_10x10_sse(a, b, c)
#else
        call dmatzvec_n_10x10_avx(a, b, c)
#endif


      else

        stop "dmatzvec_10x10_unrolled : Error, invalid (n, m) for trans = 'N'."

      end if

    else if (trans == "T") then

      if (m == 1) then

        bj = b(1)
        do i = 1, n ! Do not unroll : Seldom used.
          c(i) = a(1, i)*bj
        end do

      else if (m == 10) then

!         b1  = b( 1) ! Does that make sense ? Is there any gain in doing this ?
!         b2  = b( 2) ! Not enough registers anyway (if the compiler ever try to make use of them !)
!         b3  = b( 3) ! Will likely turn into indexed access on the stack instead of indexed access
!         b4  = b( 4) ! on the heap, which costs a copy and shall make no difference in perfs.
!         b5  = b( 5)
!         b6  = b( 6)
!         b7  = b( 7)
!         b8  = b( 8)
!         b9  = b( 9)
!         b10 = b(10)
!         do i = 1, n
!           c(i) = a(1, i)*b1+a(2, i)*b2+a(3, i)*b3+a(4, i)*b4+a( 5, i)*b5 &
!   &            + a(6, i)*b6+a(7, i)*b7+a(8, i)*b8+a(9, i)*b9+a(10, i)*b10
!         end do

#ifdef SSE
        call dmatzvec_t_10x10_sse(a, b, c)
#else
        call dmatzvec_t_10x10_avx(a, b, c)
#endif

      else

        stop "dmatzvec_10x10_unrolled : Error, invalid (n, m) for trans = 'T'."

      end if

    else

      stop "dmatzvec_10x10_unrolled : Error, invalid trans."

    end if

  end subroutine dmatzvec_10x10_unrolled

!
! Complex scalar product with n = 1 or 10.
!
! Unrolled by hand.
!

  complex(dp) function zdotc_10_unrolled(n, a, b)

    implicit none

! Input/output variables.

    integer, intent(in) :: n
    complex(dp), intent(in) :: a(*), b(*)

! Local variables.

    integer :: i

! Subroutine code.

    if (n == 1) then

      zdotc_10_unrolled = conjg(a(1))*b(1)

    else if (n == 10) then

      zdotc_10_unrolled = conjg(a(1))*b(1)+conjg(a(2))*b(2)+conjg(a(3))*b(3)+conjg(a(4))*b(4)+conjg(a( 5))*b( 5) &
  &                     + conjg(a(6))*b(6)+conjg(a(7))*b(7)+conjg(a(8))*b(8)+conjg(a(9))*b(9)+conjg(a(10))*b(10)

    else

      stop "zdotc_10_unrolled : Error, invalid n."

    end if

  end function zdotc_10_unrolled

end program test10
