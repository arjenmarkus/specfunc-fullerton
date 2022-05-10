! fullterton_aux.f90 --
!     Module with auxiliary routines. Probably not useful outside the context
!     of the fullerton library
!
module fullerton_aux
    use iso_fortran_env
    use, intrinsic :: ieee_arithmetic

    implicit none

    private
    public :: r1mach, i1mach, inits, csevl, r9pak, r9upak

    !
    ! Common and less common constants used in the library
    !
    real, parameter, public :: euler  = 0.5772156649015329e0
    real, parameter, public :: pi     = 3.14159265358979324e0
    real, parameter, public :: pi2    = 1.5707963267948966e0


    !
    ! sq2pil = log(sqrt(2.*pi)),  sqpi2l = log (sqrt(pi/2.))
    !
    real, parameter, public :: sq2pil = 0.91893853320467274e0
    real, parameter, public :: sqpi2l = 0.22579135264472743e0

    !
    ! pi2rec = 2/pi - 0.625
    !
    real, parameter, public :: pi2rec = 0.0116197723675813430e0

contains

! r1mach --
!     Rewrite of the classical R1MACH routine, based on the C implementation
!
real function r1mach(i)
    integer, intent(in) :: i

!  Original comments:
!
!  SINGLE-PRECISION MACHINE CONSTANTS
!  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
!  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
!  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
!  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
!  R1MACH(5) = LOG10(B)
!
    select case (i)
        case (1)
            r1mach = tiny(1.0)

        case (2)
            r1mach = huge(1.0)

        case (3)
            r1mach = epsilon(1.0) / radix(1.0)

        case (4)
            r1mach = epsilon(1.0)

        case (5)
            r1mach = log10( real(radix(1.0), kind=kind(1.0)) )

        case default
            write( error_unit, * ) "R1MACH: invalid argument - should be between 1 and 5"
            error stop
    end select
end function r1mach

! i1mach --
!     Rewrite of the classical I1MACH routine, based on the C implementation
!
integer function i1mach(i)
    integer, intent(in) :: i

!   Original comments:
!
!    I1MACH( 1) = THE STANDARD INPUT UNIT.
!    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
!    I1MACH( 3) = THE STANDARD PUNCH UNIT.
!    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
!    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
!    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
!    INTEGERS HAVE FORM SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!    I1MACH( 7) = A, THE BASE.
!    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
!    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
!    FLOATS HAVE FORM  SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!               WHERE  EMIN .LE. E .LE. EMAX.
!    I1MACH(10) = B, THE BASE.
!  SINGLE-PRECISION
!    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
!  DOUBLE-PRECISION
!    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
!    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
!    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
!
    select case (i)
        case (1)
            i1mach = input_unit
        case (2)
            i1mach = output_unit
        case (3)
            i1mach = -1 ! punch cards not supported
        case (4)
            i1mach = error_unit
        case (5)
            i1mach = bit_size(1)
        case (6)
            i1mach = 4 ! Fixed to the value from the classic routine
        case (7)
            i1mach = radix(1)
        case (8)
            i1mach = digits(1)
        case (9)
            i1mach = huge(1)
        case (10)
            i1mach = radix(1.0)
        case (11)
            i1mach = digits(1.0)
        case (12)
            i1mach = exponent(tiny(1.0))
        case (13)
            i1mach = exponent(huge(1.0))
        case (14)
            i1mach = digits(1.0d0)
        case (15)
            i1mach = exponent(tiny(1.0d0))
        case (16)
            i1mach = exponent(huge(1.0d0))
        case default
            write(error_unit, *) "I1MACH: invalid argument - should be between 1 and 16"
            error stop
    end select
end function i1mach

! inits --
!     Original:
!
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     initialize the orthogonal series so that inits is the number of terms
!     needed to insure the error is no larger than eta.  ordinarily, eta
!     will be chosen to be one-tenth machine precision.
!
integer function inits (os, nos, eta)
!
!             input arguments --
! os     array of nos coefficients in an orthogonal series.
! nos    number of coefficients in os.
! eta    requested accuracy of series.

    integer, intent(in) :: nos
    real, intent(in)    :: os(nos)
    real, intent(in)    :: eta

    real                :: err
    integer             :: i, ii

    ! Note
    ! This condition seems superfluous
    ! if (nos.lt.1) call seteru ( "inits   number of coefficients lt 1", 35, 2, 2)

    i   = 0
    err = 0.
    do ii = 1,nos
        i   = nos + 1 - ii
        err = err + abs(os(i))
        if (err > eta) exit
    enddo

    ! Note: useful check
    if (i == nos) then
        write(error_unit, * ) "inits: eta may be too small"
    endif

    inits = i

end function inits

! csevl --
!     Original:
!
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate the n-term chebyshev series cs at x.  adapted from
!     r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
!     and parker, chebyshev polys in numerical analysis, oxford press, p.56.
!
real function csevl (x, cs, n)
!
!             input arguments --
! x      value at which the series is to be evaluated.
! cs     array of n terms of a chebyshev series.  in eval-
!        uating cs, only half the first coef is summed.
! n      number of terms in array cs.
!

    integer, intent(in) :: n
    real, intent(in)    :: x
    real, intent(in)    :: cs(n)

    integer             :: i, ni
    real                :: b0, b1, b2, twox

    ! Note:
    ! The error conditions seem superfluous
    !
    ! if (n.lt.1) call seteru ("csevl   number of terms le 0", 28, 2,2)
    ! if (n.gt.1000) call seteru ("hcsevl   number of terms gt 1000", 31, 3, 2)
    ! if (x.lt.(-1.1) .or. x.gt.1.1) call seteru ("csevl   x outside (-1,+1)", 25, 1, 1)

    b1 = 0.
    b0 = 0.
    twox = 2.*x
    do i = 1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
    enddo

    csevl = 0.5 * (b0-b2)

end function csevl

! r9pak --
!     Original:
!     december 1979 edition. w. fullerton, c3, los alamos scientific lab.
!
!     pack a base 2 exponent into floating point number x.  this routine
!     is almost the inverse of r9upak.  it is not exactly the inverse,
!     because abs(x) need not be between 0.5 and 1.0.  if both r9pak and
!     2.0**n were known to be in range, we could compute
!                 r9pak = x * 2.0**n .
!
real function r9pak (y, n)
    real, intent(in)    :: y
    integer, intent(in) :: n

    integer             :: ny

    call r9upak (y, r9pak, ny)

    r9pak = set_exponent( r9pak, (n+ny) )

end function r9pak

! ru9pak --
!     Original:
!     august 1980 portable edition.  w. fullerton, los alamos scientific lab
!
!     unpack floating point number x so that x = y * 2.0**n, where
!     0.5 .le. abs(y) .lt. 1.0 .
!
subroutine r9upak (x, y, n)
    implicit none

    real, intent(in)     :: x
    real, intent(out)    :: y
    integer, intent(out) :: n

    !
    ! Code replaced by two intrinsic functions
    ! Testing has shown these to be equivalent to the original
    !
    y = fraction( x )
    n = exponent( x )
end subroutine r9upak

end module fullerton_aux
