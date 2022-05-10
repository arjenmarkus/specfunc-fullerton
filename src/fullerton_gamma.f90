! fullerton_gamma.f90 --
!     Module cotaining implementations of the gamma function and various related functions
!
module fullerton_gamma
    use ieee_arithmetic
    use fullerton_aux
    use fullerton_gamma_aux

    implicit none
    private
    public :: fgamma, log_gamma, sign_log_gamma, inc_gamma, inc_gamma_compl, gamma_tricomi, reciprocal_gamma, digamma

    ! Do distinguish between the intrinsic gamma function and this one
    interface fgamma
        module procedure gamma
    end interface

    interface log_gamma
        module procedure alngam
    end interface

    interface sign_log_gamma
        module procedure algams
    end interface

    interface inc_gamma
        module procedure gami
    end interface

    interface inc_gamma_compl
        module procedure gamic
    end interface

    interface gamma_tricomi
        module procedure gamit
    end interface

    interface reciprocal_gamma
        module procedure gamr
    end interface

    interface digamma
        module procedure psi
    end interface

contains

! algams --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate log abs (gamma(x)) and return the sign of gamma(x) in sgngam.
!     sgngam is either +1.0 or -1.0.
!
subroutine algams (x, algam, sgngam)

    real, intent(in) :: x
    real, intent(out) :: algam, sgngam

    integer :: int

    algam = alngam(x)
    sgngam = 1.0

    if ( x <= 0.0 ) then
        int = mod (-aint(x), 2.0) + 0.1
        if ( int == 0 ) sgngam = -1.0
    endif
end subroutine algams

! gamma --
!     Original:
!     jan 1984 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function gamma (x)
    real, intent(in) :: x

    real, save :: gcs(23) = [   &
        .008571195590989331e0,  &
        .004415381324841007e0,  &
        .05685043681599363e0,   &
       -.004219835396418561e0,  &
        .001326808181212460e0,  &
       -.0001893024529798880e0, &
        .0000360692532744124e0, &
       -.0000060567619044608e0, &
        .0000010558295463022e0, &
       -.0000001811967365542e0, &
        .0000000311772496471e0, &
       -.0000000053542196390e0, &
        .0000000009193275519e0, &
       -.0000000001577941280e0, &
        .0000000000270798062e0, &
       -.0000000000046468186e0, &
        .0000000000007973350e0, &
       -.0000000000001368078e0, &
        .0000000000000234731e0, &
       -.0000000000000040274e0, &
        .0000000000000006910e0, &
       -.0000000000000001185e0, &
        .0000000000000000203e0  ]

    integer, save   :: ngcs = 0
    real, save      :: xmin = 0.0
    real, save      :: xmax = 0.0
    real, save      :: xsml = 0.0
    real, save      :: dxrel= 0.0

    real            :: y, sinpiy
    integer         :: i, n

    !
    ! Initialisation:
    ! initialize.  find legal bounds for x, and determine the number of
    ! terms in the series required to attain an accuracy ten times better
    ! than machine precision.
    !
    if ( ngcs == 0 ) then
        ngcs = inits (gcs, 23, 0.1*r1mach(3))

        call r9gaml (xmin, xmax)
        xsml = exp (max (log (r1mach(1)), -log(r1mach(2))) + 0.01)
        dxrel = sqrt (r1mach(4))
    endif

    !
    ! Argument ranges:
    !                x  < -1.0
    !            abs(x) <  10.0
    !            abs(x) >= 10.0
    !
    y = abs(x)
    if ( y <= 10.0 ) then
        !
        ! compute gamma(x) for abs(x) <= 10.0.  reduce interval and
        ! find gamma(1+y) for 0. <= y <= 1. first of all.
        !
        n = x
        if ( x < 0.0 ) n = n - 1
        y = x - float(n)
        n = n - 1
        gamma = 0.9375 + csevl(2.*y-1., gcs, ngcs)

        if ( n == 0 ) then
            return

        elseif (n < 0 ) then
            !
            ! compute gamma(x) for x .lt. 1.
            !
            n = -n
            if ( x == 0.0 ) then
                gamma = ieee_value( x, ieee_positive_inf )

            elseif ( x < 0.0 .and. x+float(n-2) == 0.0 ) then
                gamma = ieee_value( x, ieee_positive_inf )

            else
                do i = 1,n
                    gamma = gamma / (x+float(i-1))
                enddo
            endif
        else
            !
            ! gamma(x) for x >= 2.
            !
            do i = 1,n
                gamma = (y+float(i))*gamma
            enddo
        endif

    else
        !
        ! compute gamma(x) for abs(x) .gt. 10.0.  recall y = abs(x).
        !
        if ( x > xmax ) then
            gamma = ieee_value( x, ieee_positive_inf )
        else

            if ( x < xmin ) then
                gamma = 0.
            else
                gamma = exp((y-0.5)*log(y) - y + sq2pil + r9lgmc(y) )

                if ( x < 0.0 ) then
                    sinpiy = sin (pi*y)

                    if ( sinpiy == 0.0 ) then
                        gamma = ieee_value( x, ieee_positive_inf )
                    else
                        gamma = -pi / (y*sinpiy*gamma)
                    endif
                endif
            endif
       endif
    endif
end function gamma

! gamr --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!     this routine, not gamma(x), should be the fundamental one.
!
real function gamr (x)
    real, intent(in) :: x

    real    :: alngx, sgngx
    integer :: ir, irold

    gamr = 0.0

    if ( x <= 0.0 .and. aint(x) == x ) then
        return

    else
        if ( abs(x) <= 10.0 ) then
            gamr = 1.0/gamma(x)
        else
            call algams (x, alngx, sgngx)
            gamr = sgngx * exp(-alngx)
        endif
    endif
end function gamr

! gami --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate the incomplete gamma function defined by
!
!     gami = integral from t = 0 to x of exp(-t) * t**(a-1.0) .
!
!     gami is evaluated for positive values of a and non-negative values
!     of x.  a slight deterioration of 2 or 3 digits accuracy will occur
!     when gami is very large or very small, because logarithmic variables
!     are used.
!
real function gami (a, x)
    real, intent(in) :: a, x

    real :: factor

    if ( a  <= 0.0 .or. x < 0.0 ) then
        gami = ieee_value( x, ieee_quiet_nan )

    else
        gami = 0.0

        if ( x > 0.0 ) then
            !
            ! the only error possible in the expression below is a fatal overflow.
            !
            factor = exp (alngam(a) + a*log(x) )

            gami   = factor * gamit(a, x)
        endif
    endif
end function gami

! gamic --
!     Original:
!    july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!    evaluate the complementary incomplete gamma function
!
!    gamic = integral from t = x to infinity of exp(-t) * t**(a-1.)  .
!
!    gamic is evaluated for arbitrary real values of a and for non-negative
!    values of x (even though gamic is defined for x .lt. 0.0), except that
!    for x = 0 and a .le. 0.0, gamic is undefined.
!
!         a slight deterioration of 2 or 3 digits accuracy will occur when
!    gamic is very large or very small in absolute value, because log-
!    arithmic variables are used.  also, if the parameter a is very close
!    to a negative integer (but not a negative integer), there is a loss
!    of accuracy, which is reported if the result is less than half
!    machine precision.
!
!    ref. -- w. gautschi, an evaluation procedure for incomplete gamma
!    functions, acm trans. math. software.
!
real function gamic (a, x)
    real, intent(in) :: a, x

    real, save :: eps    = 0.0
    real, save :: sqeps  = 0.0
    real, save :: alneps = 0.0
    real, save :: bot    = 0.0

    real       :: aeps, algap1, alngs, alx, e, gstar, h, fm, sga, sgng, sgngam, sgngs, t
    integer    :: ma, izero

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps    = 0.5*r1mach(3)
        sqeps  = sqrt(r1mach(4))
        alneps = -log(r1mach(3))
        bot    = log(r1mach(1))
    endif

    !
    ! Argument ranges:
    !           x <  0.0
    !           x == 0.0
    !           x >= 1.0
    !
    if ( x <= 0.0 ) then
        gamic = ieee_value( x, ieee_quiet_nan )

    elseif ( x == 0.0 ) then
        if ( a <= 0.0 ) then
            gamic = ieee_value( x, ieee_quiet_nan )
        else
            gamic = exp (alngam(a+1.0) - log(a))
        endif
    else
        alx   = log(x)
        sga   = 1.0
        if (a /= 0.0) sga = sign (1.0, a)

        ma    = a + 0.5*sga
        aeps  = a - float(ma)

        izero = 0

        if ( x < 1.0 ) then
            if ( a <= 0.5 .and. abs(aeps) < 0.001 ) then
                fm = -ma
                e  = 2.0
                if ( fm > 1.0 ) then
                    e = 2.0*(fm+2.0)/(fm*fm-1.0)
                endif

                e = e - alx*x**(-0.001)

                if ( e*abs(aeps) <= eps ) then
                    gamic = r9gmic (a, x, alx)
                    return
                endif
            endif

            !
            ! Fall through for the last part
            !
            call algams (a+1.0, algap1, sgngam)
            gstar = r9gmit (a, x, algap1, sgngam, alx)
            if ( gstar == 0.0 ) then
                izero = 1
            else
                alngs = log (abs(gstar))
                sgngs = sign (1.0, gstar)
            endif
        else
            if ( a < x ) then
                gamic = exp (r9lgic(a, x, alx))
                return
            else
                sgngam = 1.0
                algap1 = alngam (a+1.0)
                sgngs = 1.0
                alngs = r9lgit (a, x, algap1)
            endif
        endif

        !
        ! evaluation of gamic(a,x) in terms of tricomi-s incomplete gamma fn.
        !
        h = 1.0
        if ( izero /= 1 ) then
            t = a*alx + alngs
            if ( t > alneps ) then
                sgng  = -sgngs * sga * sgngam
                t     = t + algap1 - log(abs(a))
                gamic = sgng * exp(t)
                return
            elseif ( t > -alneps ) then
                h     = 1.0 - sgngs*exp(t)
            endif
        endif

        sgng = sign (1.0, h) * sga * sgngam
        t    = log(abs(h)) + algap1 - log(abs(a))
        gamic = sgng * exp(t)
    endif
end function gamic

! gamit --
!     Original:
!    july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!    evaluate tricomi-s incomplete gamma function defined by
!
!    gamit = x**(-a)/gamma(a) * integral t = 0 to x of exp(-t) * t**(a-1.)
!
!    and analytic continuation for a .le. 0.0.  gamma(x) is the complete
!    gamma function of x.  gamit is evaluated for arbitrary real values of
!    a and for non-negative values of x (even though gamit is defined for
!    x .lt. 0.0).
!
!         a slight deterioration of 2 or 3 digits accuracy will occur when
!    gamit is very large or very small in absolute value, because log-
!    arithmic variables are used.  also, if the parameter a is very close
!    to a negative integer (but not a negative integer), there is a loss
!    of accuracy, which is reported if the result is less than half
!    machine precision.
!
!    ref. -- w. gautschi, an evaluation procedure for incomplete gamma
!    functions, acm trans. math. software.
!
real function gamit (a, x)
    real, intent(in) :: a, x

    real, save :: alneps = 0.0
    real, save :: sqeps  = 0.0
    real, save :: bot    = 0.0

    real       :: aeps, ainta, algap1, alng, alx, h, sga, sgngam, t

    !
    ! Initialisation
    !
    if ( alneps == 0.0 ) then
        alneps = -log(r1mach(3))
        sqeps  = sqrt(r1mach(4))
        bot    = log(r1mach(1))
    endif

    !
    ! Argument ranges:
    !     x < 0.0
    !     0.0 <= x <= 1.0
    !     x > 1.0
    !
    if ( x < 0.0 ) then
        gamit = ieee_value( x, ieee_quiet_nan )

    else
        if ( x /= 0.0 ) then
            alx = log(x)
        endif

        sga = 1.0
        if ( a /= 0.0 ) then
            sga = sign (1.0, a)
        endif

        ainta = aint (a+0.5*sga)
        aeps = a - ainta

        if ( x == 0.0 ) then
            gamit = 0.0
            if ( ainta > 0.0 .or. aeps /= 0.0 ) then
                gamit = gamr(a+1.0)
            endif
        elseif ( x <= 1.0 ) then
            if ( a  >= -0.5 .or. aeps /= 0.0) then
                call algams (a+1.0, algap1, sgngam)
            endif
            gamit = r9gmit (a, x, algap1, sgngam, alx)
        else
            if ( a >= x ) then
                t = r9lgit (a, x, alngam(a+1.0))
                gamit = exp(t)

            else
                alng = r9lgic (a, x, alx)
                !
                ! evaluate gamit in terms of log(gamic(a,x))
                !
                h = 1.0
                if ( aeps /= 0.0 .or. ainta > 0.0 ) then
                    call algams (a+1.0, algap1, sgngam)
                    t = log(abs(a)) + alng - algap1

                    if ( t <= alneps ) then
                        if ( t > -alneps ) then
                             h = 1.0 - sga*sgngam*exp(t)
                        endif
                        t = -a*alx + log(abs(h))
                        gamit = sign (exp(t), h)
                    else
                        t = t - a*alx
                        gamit = -sga*sgngam*exp(t)
                    endif
                endif
            endif
        endif
    endif
end function gamit

! psi --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
real function psi (x)
    real, intent(in) :: x

!
! series for psi        on the interval  0.          to  1.00000d+00
!                                        with weighted error   2.03e-17
!                                         log weighted error  16.69
!                               significant figures required  16.39
!                                    decimal places required  17.37
!
    real, save :: psics(23) = [ &
        -.038057080835217922e0, &
         .49141539302938713e0,  &
        -.056815747821244730e0, &
         .008357821225914313e0, &
        -.001333232857994342e0, &
         .000220313287069308e0, &
        -.000037040238178456e0, &
         .000006283793654854e0, &
        -.000001071263908506e0, &
         .000000183128394654e0, &
        -.000000031353509361e0, &
         .000000005372808776e0, &
        -.000000000921168141e0, &
         .000000000157981265e0, &
        -.000000000027098646e0, &
         .000000000004648722e0, &
        -.000000000000797527e0, &
         .000000000000136827e0, &
        -.000000000000023475e0, &
         .000000000000004027e0, &
        -.000000000000000691e0, &
         .000000000000000118e0, &
        -.000000000000000020e0  ]
!
! series for apsi       on the interval  0.          to  2.50000d-01
!                                        with weighted error   5.54e-17
!                                         log weighted error  16.26
!                               significant figures required  14.42
!                                    decimal places required  16.86
!
    real, save :: apsics(16) = [ &
        -.0204749044678185e0,    &
        -.0101801271534859e0,    &
         .0000559718725387e0,    &
        -.0000012917176570e0,    &
         .0000000572858606e0,    &
        -.0000000038213539e0,    &
         .0000000003397434e0,    &
        -.0000000000374838e0,    &
         .0000000000048990e0,    &
        -.0000000000007344e0,    &
         .0000000000001233e0,    &
        -.0000000000000228e0,    &
         .0000000000000045e0,    &
        -.0000000000000009e0,    &
         .0000000000000002e0,    &
        -.0000000000000000e0     ]

    integer, save   :: ntpsi  = 0
    integer, save   :: ntapsi = 0
    real, save      :: xbig   = 0.0
    real, save      :: dxrel  = 0.0

    real            :: y, aux
    integer         :: i, n

    !
    ! Initialisation
    !
    if ( ntpsi == 0 ) then
        ntpsi  = inits (psics, 23, 0.1*r1mach(3))
        ntapsi = inits (apsics, 16, 0.1*r1mach(3))

        xbig   = 1.0/sqrt(r1mach(3))
        dxrel   = sqrt (r1mach(4))
    endif

    !
    ! Argument ranges:
    !            abs(x) <  2.0
    !            abs(x) >= 2.0
    !
    y = abs(x)
    if ( y < 2.0 ) then
        n = x
        if ( x < 0.0 ) n = n - 1
        y = x - float(n)
        n = n - 1
        psi = csevl (2.*y-1., psics, ntpsi)

        if ( n /= 0 ) then
            n = -n
            if ( x == 0.0 .or. ( x < 0.0 .and. x+float(n-2) == 0.0 ) ) then
                psi = ieee_value( x, ieee_positive_inf )
            else
                do i = 1,n
                    psi = psi - 1.0/(x+float(i-1))
                enddo
            endif
        endif
    else
        aux = 0.
        if ( y < xbig ) aux = csevl (8./y**2-1., apsics, ntapsi)

         if ( x < 0.0 ) then
             psi = log(abs(x)) - 0.5/x + aux - pi/tan(pi*x)
         else
             psi = log(x) - 0.5/x + aux
         endif
    endif
end function psi

end module fullerton_gamma
