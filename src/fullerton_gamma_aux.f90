! fullerton_gamma_aux.f90 --
!     Module with auxiliary functions for the gamma function and relatives
!
module fullerton_gamma_aux
    use ieee_arithmetic
    use fullerton_aux

    implicit none

    private :: gamma ! To ensure the Fullerton version of the gamma function is used

contains

! r9lgmc --
!     Original:
!     august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     compute the log gamma correction factor for x .ge. 10.0 so that
!     log (gamma(x)) = log(sqrt(2*pi)) + (x-.5)*log(x) - x + r9lgmc(x)
!
real function r9lgmc (x)
    real, intent(in) :: x

!
! series for algm       on the interval  0.          to  1.00000d-02
!                                        with weighted error   3.40e-16
!                                         log weighted error  15.47
!                               significant figures required  14.39
!                                    decimal places required  15.86
!
    real, save :: algmcs(6) = [ &
         .166638948045186e0,    &
        -.0000138494817606e0,   &
         .0000000098108256e0,   &
        -.0000000000180912e0,   &
         .0000000000000622e0,   &
        -.0000000000000003e0    ]

    integer, save :: nalgm = 0
    real, save    :: xbig  = 0.0
    real, save    :: xmax  = 0.0

    !
    ! Initialisation
    !
    if ( nalgm == 0 ) then
        nalgm = inits (algmcs, 6, r1mach(3))
        xbig = 1.0/sqrt(r1mach(3))
        xmax = exp (min(log(r1mach(2)/12.0), -log(12.0*r1mach(1))) )
    endif

    !
    ! Argument ranges:
    !     x < 10.0
    !     10.0 <= x < xmax
    !     x >= xmax
    !
    if ( x < 10.0 ) then
        r9lgmc = ieee_value( x, ieee_quiet_nan )

    elseif ( x < xmax ) then
        r9lgmc = 1.0/(12.0*x)
        if (x < xbig ) then
            r9lgmc = csevl (2.0*(10./x)**2-1., algmcs, nalgm)/x
        endif

    else
        r9lgmc = 0.0
    endif
end function r9lgmc

! r9gaml --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     calculate the minimum and maximum legal bounds for x in gamma(x).
!     xmin and xmax are not the only bounds, but they are the only non-
!     trivial ones to calculate.
!
!                 output arguments --
!     xmin   minimum legal value of x in gamma(x).  any smaller value of
!            x might result in underflow.
!     xmax   maximum legal value of x in gamma(x).  any larger value will
!            cause overflow.
!
subroutine r9gaml (xmin, xmax)
    real, intent(out) :: xmin, xmax

    real     :: alnsml, alnbig, xold, xln
    integer  :: i

    !
    ! Calculation
    !
    alnsml = log(r1mach(1))
    xmin   = -alnsml

    do i = 1,10
        xold = xmin
        xln = log(xmin)
        xmin = xmin - xmin*((xmin+0.5)*xln - xmin - 0.2258 + alnsml) / (xmin*xln + 0.5)
        if ( abs(xmin-xold) < 0.005 ) then
            exit
        endif
    enddo

    xmin = -xmin + 0.01

    alnbig = log(r1mach(2))
    xmax   = alnbig

    do i = 1,10
        xold = xmax
        xln = log(xmax)
        xmax = xmax - xmax*((xmax-0.5)*xln - xmax + 0.9189 - alnbig) / (xmax*xln - 0.5)
        if ( abs(xmax-xold) < 0.005 ) then
            exit
        endif
    enddo

    xmax = xmax - 0.01
    xmin = max( xmin, -xmax+1.0 )
end subroutine r9gaml

! r9lgit --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     compute the log of tricomi-s incomplete gamma function with perron-s
!     continued fraction for large x and for a .ge. x.
!
real function r9lgit (a, x, algap1)
    real, intent(in) :: a, x, algap1

    real, save :: eps   = 0.0
    real, save :: sqeps = 0.0

    real       :: ax, a1x, r, p, s, fk, t, hstar
    integer    :: k

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps   = 0.5*r1mach(3)
        sqeps = sqrt(r1mach(4))
    endif

    !
    ! Calculation
    !
    if ( x <= 0.0 .or. a < x) then
        r9lgit = ieee_value( x, ieee_quiet_nan )
    else
        ax  = a + x
        a1x = ax + 1.0
        r   = 0.0
        p   = 1.0
        s   = p

        do k = 1,200
            fk = k
            t = (a+fk)*x*(1.0+r)
            r = t/((ax+fk)*(a1x+fk)-t)
            p = r*p
            s = s + p
            if ( abs(p) < eps*s ) then
                exit
            endif
        enddo
    endif

    hstar  = 1.0 - x*s/a1x
    r9lgit = -x - algap1 - log(hstar)

end function r9lgit

! r9glic --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     compute the log complementary incomplete gamma function for large x
!     and for a .le. x.
!
real function r9lgic (a, x, alx)
    real, intent(in) :: a, x, alx

    real, save :: eps = 0.0

    real       :: xma, xpa, r, p, s, fk, t
    integer    :: k

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps = 0.5*r1mach(3)
    endif

    !
    ! Calculation
    !
    xpa = x + 1.0 - a
    xma = x - 1.0 - a

    r   = 0.0
    p   = 1.0
    s   = p

    do k = 1,200
        fk = k
        t  = fk*(a-fk)*(1.0+r)
        r  = -t/((xma+2.0*fk)*(xpa+2.0*fk)+t)
        p  = r*p
        s  = s + p
        if ( abs(p) < eps*s) then
            exit
        endif
    enddo

    r9lgic = a*alx - x + log(s/xpa)

end function r9lgic

! r9gmit --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     compute tricomi-s incomplete gamma function for small x.
!
real function r9gmit (a, x, algap1, sgngam, alx)
    real, intent(in) :: a, x, algap1, sgngam, alx

    real, save :: eps = 0.0
    real, save :: bot = 0.0

    real       :: ae, aeps, t, te, s, fk, algs, alg2, sgng2
    integer    :: k, m, ma

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps = 0.5*r1mach(3)
        bot = log(r1mach(1))
    endif

    !
    ! Calculation
    !
    if ( x <= 0.0 ) then
        r9gmit = ieee_value( x, ieee_quiet_nan )

    else
        ma = a + 0.5
        if ( a < 0.0 ) then
            ma = a - 0.5
        endif
        aeps = a - float(ma)

        ae = a
        if ( a < -0.5) then
            ae = aeps
        endif

        t  = 1.0
        te = ae
        s  = t
        do k = 1,200
            fk = k
            te = -x*te/fk
            t  = te/(ae+fk)
            s  = s + t
            if ( abs(t) < eps*abs(s) ) then
                exit
            endif
        enddo

        if ( a >= -0.5 ) then
            algs   = -algap1 + log(s)
            r9gmit = exp(algs)

        else
            algs = -alngam(1.0+aeps) + log(s)
            s    = 1.0
            m    = -ma - 1
            if ( m /= 0 ) then
                t = 1.0
                do k = 1,m
                    t = x*t/(aeps-float(m+1-k))
                    s = s + t
                    if ( abs(t) < eps*abs(s) ) then
                        exit
                    endif
                enddo

                r9gmit = 0.0
                algs   = -float(ma)*log(x) + algs

                if ( s == 0.0 .or. aeps == 0.0 ) then
                    r9gmit = exp(algs)

                else
                    sgng2 = sgngam*sign(1.0,s)
                    alg2  = -x - algap1 + log(abs(s))

                    if ( alg2 > bot ) then
                        r9gmit = sgng2*exp(alg2)
                    endif
                    if ( algs > bot ) then
                        r9gmit = r9gmit + exp(algs)
                    endif
                endif
            endif
        endif
    endif
end function r9gmit

! r9gmic --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     compute the complementary incomplete gamma function for a near
!     a negative integer and for small x.
!
real function r9gmic (a, x, alx)
    real, intent(in) :: a, x, alx

    real, save      :: eps   = 0.0
    real, save      :: bot   = 0.0

    real            :: fk, fm, t, te, s, fkp1, alng, sgng
    integer         :: m, ma, k, mm1

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps = 0.5*r1mach(3)
        bot = log(r1mach(1))
    endif

    !
    ! Calculation
    !
    if ( a > 0.0 .or. x <= 0.0 ) then
        r9gmic = ieee_value( x, ieee_quiet_nan )

    else
        ma  = a - 0.5
        fm = -ma
        m  = -ma

        te = 1.0
        t  = 1.0
        s  = t
        do k = 1,200
            fkp1 = k + 1
            te   = -x*te/(fm+fkp1)
            t    = te/fkp1
            s = s + t
            if ( abs(t) < eps*s ) then
                exit
            endif
        enddo

        r9gmic = -alx - euler + x*s/(fm+1.0)

        if ( m /= 0 ) then

            if ( m == 1 ) then
                r9gmic = -r9gmic - 1.0 + 1.0/x
            else
                te  = fm
                t   = 1.0
                s   = t
                mm1 = m - 1

                do k = 1,mm1
                    fk = k
                    te = -x*te/fk
                    t  = te/(fm-fk)
                    s  = s + t
                    if ( abs(t) < eps*abs(s) ) then
                        exit
                    endif
                enddo

                do k = 1,m
                    r9gmic = r9gmic + 1.0/float(k)
                enddo

                sgng = 1.0
                if ( mod(m,2) == 1 ) then
                    sgng = -1.0
                endif

                alng = log(r9gmic) - alngam(fm+1.0)

                r9gmic = 0.0
                if ( alng > bot ) then
                    r9gmic = sgng*exp(alng)
                endif
                if ( s /= 0.0 ) then
                    r9gmic = r9gmic + sign (exp(-fm*alx+log(abs(s)/fm)), s)
                endif
            endif
        endif
    endif
end function r9gmic

! alngam --
!     Original:
!     august 1980 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function alngam(x)
    real, intent(in) :: x

    real, save      :: xmax   = 0.0
    real, save      :: dxrel  = 0.0

    real            :: y, sinpiy

    !
    ! Initialisation
    !
    if ( xmax == 0.0 ) then
        xmax  = r1mach(2)/log(r1mach(2))
        dxrel = sqrt (r1mach(4))
    endif

    !
    ! Argument ranges:
    !            abs(x) <= 10.0
    !            abs(x) >  10.0
    !
    y = abs(x)

    if ( y <= 10.0 ) then
        !
        ! log (abs (gamma(x))) for  abs(x) <= 10.0
        !
        alngam = log (abs (gamma(x)))

    else
        !
        ! log (abs (gamma(x))) for abs(x) > 10.0
        !
        if ( y > xmax ) then
            alngam = ieee_value( x, ieee_positive_inf )
        else
            if ( x > 0.0 ) then
                alngam = sq2pil + (x-0.5)*log(x) - x + r9lgmc(y)

            else
                sinpiy = abs (sin(pi*y))

                if ( sinpiy ==0.0 ) then
                    alngam = ieee_value( x, ieee_positive_inf )
                else
                    alngam = sqpi2l + (x-0.5)*log(y) - x - log(sinpiy) - r9lgmc(y)
                endif
            endif
        endif
    endif
end function alngam

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

    integer, save :: ngcs = 0
    real, save    :: xmin = 0.0
    real, save    :: xmax = 0.0
    real, save    :: xsml = 0.0
    real, save    :: dxrel= 0.0

    real          :: y, sinpiy
    integer       :: i, n

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

end module fullerton_gamma_aux
