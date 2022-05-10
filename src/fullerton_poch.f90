! fullerton_poch.f90 --
!     Module for the Pochhammer symbol and related functions
!
module fullerton_poch
    use ieee_arithmetic
    use fullerton_aux
    use fullerton_misc
    use fullerton_gamma_aux
    use fullerton_gamma

    implicit none

    interface factorial
        module procedure fac
    end interface

    private
    public :: poch, poch1, factorial

contains

! poch --
!     Original:
!     august 1980 edition.  w. fullerton, c3, los alamos scientific lab.
!     error handling when poch (a, x) is less than half precision is
!     probably incorrect.  grossly wrong arguments are not handled right.
!
!     evaluate a generalization of pochhammer-s symbol
!     (a)-sub-x = gamma(a+x)/gamma(a).  for x a non-negative integer,
!     poch(a,x) is just pochhammer-s symbol.
!
real function poch (a, x)
    real, intent(in) :: a, x

    real, save       :: eps   = 0.0
    real, save       :: sqeps = 0.0

    real             :: ax, absa, absax, alnga, sgnga, alngax, sgngax, b
    real             :: cospia, cospix, den, err, errpch, sinpia, sinpix
    integer          :: i, n

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps = r1mach(4)
        sqeps = sqrt(eps)
    endif

    !
    ! Argument ranges:
    !              a+x  <=  0
    !              a+x  >   0
    !              a+x  = integer
    !              a+x  /= integer
    !
    ax = a + x

    if ( ax <= 0.0 .and. aint(ax) == ax ) then

        !
        ! Final check: a positive or not an integer
        !
        if (a > 0.0 .or. aint(a) /= a ) then
            poch = ieee_value( x, ieee_quiet_nan )
            return
        endif

        !
        ! we know here that both a+x and a are non-positive integers.
        !
        poch = 1.0

        if ( x == 0.0 ) then
            return
        else
            n = x
            if ( min(a+x,a) >= -20.0 ) then
                poch = (-1.0)**n * fac(-int(a))/fac(-int(a)-n)
            else
                poch = (-1.0)**n * exp ((a-0.5)*lnrel(x/(a-1.0)) + x*alog(-a+1.0-x) - x + r9lgmc(-a+1.0) - r9lgmc(-a-x+1.0) )
            endif
        endif
    else
        !
        ! here we know a+x is not zero or a negative integer.
        !
        poch = 0.0
        if ( a <= 0.0 .and. aint(a) == a ) then
            return
        endif

        n = abs(x)
        if ( float(n) == x .and. n <= 20 ) then

            !
            ! x is a small non-positive integer, presumably a common case.
            !
            poch = 1.0
            do i = 1,n
                poch = poch * (a+float(i-1))
            enddo
        else

            !
            ! x is not an integer or large
            !
            absax = abs(a+x)
            absa = abs(a)

            if ( max(absax,absa) <= 20.0 ) then
                poch = fgamma(a+x)*reciprocal_gamma(a)

                ! Note:
                ! error handling above is probably bad when a almost = -n and when x is
                ! small.  maybe use reflection formula below in modified form.
            else

                if ( abs(x) <= 0.5*absa ) then

                    !
                    ! here abs(x) is small and both abs(a+x) and abs(a) are large.  thus,
                    ! a+x and a must have the same sign.  for negative a, we use
                    ! gamma(a+x)/gamma(a) = gamma(-a+1)/gamma(-a-x+1) * sin(pi*a)/sin(pi*(a+x))
                    !
                    b = a
                    if ( b < 0.0 ) then
                        b = -a - x + 1.0
                    endif

                    poch = exp ((b-0.5)*lnrel(x/b) + x*log(b+x) - x + r9lgmc(b+x) - r9lgmc(b) )

                    if ( a < 0.0 .and. poch == 0.0 ) then
                        cospix = cos (pi*x)
                        sinpix = sin (pi*x)
                        cospia = cos (pi*a)
                        sinpia = sin (pi*a)

                        errpch = abs(x)*(1.0+alog(b))
                        den    = cospix + cospia*sinpix/sinpia
                        err    = (abs(x)*(abs(sinpix)+abs(cospia*cospix/sinpia)) + abs(a*sinpix)/sinpia**2)*pi
                        err    = errpch + err/abs(den)

                        poch = poch/den
                    else

                        call sign_log_gamma (a+x, alngax, sgngax)
                        call sign_log_gamma (a, alnga, sgnga)
                        poch = sgngax * sgnga * exp(alngax-alnga)
                    endif
                endif
            endif
        endif
    endif
end function poch

! poch1 --
!     Original:
!     august 1980 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate a generalization of pochhammer-s symbol for special
!     situations that require especially accurate values when x is small in
!            poch1(a,x) = (poch(a,x)-1)/x
!                       = (gamma(a+x)/gamma(a) - 1.0)/x .
!     this specification is particularly suited for stably computing
!     expressions such as
!            (gamma(a+x)/gamma(a) - gamma(b+x)/gamma(b))/x
!                 = poch1(a,x) - poch1(b,x)
!     note that poch1(a,0.0) = psi(a)
!
!     when abs(x) is so small that substantial cancellation will occur if
!     the straightforward formula is used, we  use an expansion due
!     to fields and discussed by y. l. luke, the special functions and their
!     approximations, vol. 1, academic press, 1969, page 34.
!
!     the ratio poch(a,x) = gamma(a+x)/gamma(a) is written by luke as
!            (a+(x-1)/2)**x * polynomial in (a+(x-1)/2)**(-2) .
!     in order to maintain significance in poch1, we write for positive a
!            (a+(x-1)/2)**x = exp(x*alog(a+(x-1)/2)) = exp(q)
!                           = 1.0 + q*exprel(q) .
!     likewise the polynomial is written
!            poly = 1.0 + x*poly1(a,x) .
!     thus,
!            poch1(a,x) = (poch(a,x) - 1) / x
!                       = exprel(q)*(q/x + q*poly1(a,x)) + poly1(a,x)
!
real function poch1 (a, x)
    real, intent(in) :: a, x

    !
    ! bern(i) is the 2*i bernoulli number divided by factorial(2*i).
    !
    real, save :: bern(9) = [    &
         .83333333333333333e-01, &
        -.13888888888888889e-02, &
         .33068783068783069e-04, &
        -.82671957671957672e-06, &
         .20876756987868099e-07, &
        -.52841901386874932e-09, &
         .13382536530684679e-10, &
        -.33896802963225829e-12, &
         .85860620562778446e-14  ]

    real, save :: sqtbig = 0.0
    real, save :: alneps = 0.0

    real       :: gbern(10)
    real       :: absa, absx, alnvar, b, binv, bp, gbk, poly1
    real       :: q, rho, sinpx2, sinpxx, term, trig, var, var2
    integer    :: i, ii, j, incr, k, nterms, ndx

    !
    ! Initialisation
    !
    if ( sqtbig == 0.0 ) then
        sqtbig = 1.0/sqrt(24.0*r1mach(1))
        alneps = log(r1mach(3))
    endif

    !
    ! Argument ranges:
    !              x = 0
    !         abs(x) > abs(a)
    !
    if ( x == 0.0 ) then
        poch1 = digamma(a)
    else

        absx = abs(x)
        absa = abs(a)

        if ( absx >  absa .or. absx*log(max(absa,2.0)) > 0.1 ) then
            poch1 = ( poch (a, x) - 1.0 ) / x
        else
            bp = a
            if ( a < -0.5 ) then
                bp = 1.0 - a - x
            endif

            incr = 0
            if ( bp < 10.0 ) then
                incr = 11.0 - bp
            endif

            b   = bp + float(incr)

            var    = b + 0.5*(x-1.0)
            alnvar = log(var)
            q      = x*alnvar

            poly1 = 0.0
            if ( var < sqtbig ) then
                var2     = (1.0/var)**2

                rho      = 0.5*(x+1.0)
                gbern(1) = 1.0
                gbern(2) = -rho/12.0
                term     = var2
                poly1    = gbern(2)*term

                nterms    = min( 9, int(-0.5*alneps/alnvar + 1.0) )

                do k = 2,nterms
                    gbk = 0.0
                    do j = 1,k
                       ndx = k - j + 1
                       gbk = gbk + bern(ndx)*gbern(j)
                    enddo
                    gbern(k+1) = -rho*gbk/float(k)

                    term  = term * (float(2*k-2)-x)*(float(2*k-1)-x)*var2
                    poly1 = poly1 + gbern(k+1)*term
                enddo
            endif

            poly1 = (x-1.0)*poly1
            poch1 = exprel(q)*(alnvar + q*poly1) + poly1

            !
            ! we have poch1(b,x).  but bp is small, so we use backwards recursion
            ! to obtain poch1(bp,x).
            !
            do ii = 1,incr
                i     = incr - ii
                binv  = 1.0/(bp+float(i))
                poch1 = (poch1-binv)/(1.0+x*binv)
            enddo

            if ( bp /= a ) then
                !
                ! we have poch1(bp,x), but a is lt -0.5.  we therefore use a reflection
                ! formula to obtain poch1(a,x).
                !
                sinpxx = sin(pi*x)/x
                sinpx2 = sin(0.5*pi*x)
                trig   = sinpxx*cot(pi*b) - 2.0*sinpx2*(sinpx2/x)

                poch1 = trig + (1.0 + x*trig) * poch1
            endif
        endif
    endif
end function poch1

! cot --
!     Original:
!     may 1980 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function cot (x)
    real, intent(in) :: x

    !
    ! series for cot        on the interval  0.          to  6.25000d-02
    !                                        with weighted error   3.76e-17
    !                                         log weighted error  16.42
    !                               significant figures required  15.51
    !                                    decimal places required  16.88
    !
    real, save :: cotcs(8) = [  &
         .24025916098295630e0,  &
        -.016533031601500228e0, &
        -.000042998391931724e0, &
        -.000000159283223327e0, &
        -.000000000619109313e0, &
        -.000000000002430197e0, &
        -.000000000000009560e0, &
        -.000000000000000037e0  ]

    integer, save   :: nterms = 0
    real, save      :: xmax   = 0.0
    real, save      :: xmin   = 0.0
    real, save      :: xsml   = 0.0
    real, save      :: sqeps  = 0.0

    real            :: ainty, ainty2, prodbg, y, yrem
    integer         :: ifn

    !
    ! Initialisation
    !
    if ( nterms == 0 ) then
        nterms = inits (cotcs, 8, 0.1*r1mach(3))
        xmax   = 1.0/r1mach(4)
        xsml   = sqrt (3.0*r1mach(3))
        xmin   = exp ( max(log(r1mach(1)), -log(r1mach(2))) + 0.01)
        sqeps  = sqrt (r1mach(4))
    endif

    !
    ! Argument ranges:
    !             x around n*pi
    !             y <= 0.25
    !             y <= 0.5
    !             y >  0.5
    !            (y reduced argument)
    !
    y = abs(x)

    if ( abs(x) < xmin ) then
        cot = ieee_value( x, ieee_quiet_nan )
    else
        !
        ! carefully compute y * (2/pi) = (aint(y) + rem(y)) * (.625 + pi2rec)
        ! = aint(.625*y) + rem(.625*y) + y*pi2rec  =  aint(.625*y) + z
        ! = aint(.625*y) + aint(z) + rem(z)
        !
        ainty  = aint (y)
        yrem   = y - ainty
        prodbg = 0.625*ainty
        ainty  = aint (prodbg)
        y      = (prodbg-ainty) + 0.625*yrem + y*pi2rec
        ainty2 = aint (y)
        ainty = ainty + ainty2
        y     = y - ainty2

        ifn   = mod (ainty, 2.0)
        if ( ifn == 1 ) then
            y = 1.0 - y
        endif

        if ( abs(x) > 0.5 .and. y < abs(x)*sqeps ) then
            cot = ieee_value( x, ieee_quiet_nan )
            return
        endif

        if ( y <= 0.25 ) then
            if ( y == 0.0 ) then
                cot = ieee_value( x, ieee_quiet_nan )
                return
            endif

            cot = 1.0 / y
            if ( y > xsml ) then
               cot = (0.5 + csevl (32.0*y**2-1., cotcs, nterms)) /y
            endif

        else if ( y <= 0.5 ) then
            cot = (0.5 + csevl (8.0*y**2-1., cotcs, nterms)) / (0.5*y)
            cot = (cot**2 - 1.0) * 0.5 / cot

        else
            cot = (0.5 + csevl (2.0*y**2-1., cotcs, nterms)) / (0.25*y)
            cot = (cot**2 - 1.0) * 0.5 / cot
            cot = (cot**2 - 1.0) * 0.5 / cot
        endif

        if ( x /= 0.0 ) then
            cot = sign (cot, x)
        endif
        if ( ifn == 1 ) then
            cot = -cot
        endif
    endif
end function cot

! fac --
!     Original:
!     june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function fac (n)
    integer, intent(in) :: n

    real, save :: facn(26) = [  &
        1.0e0,                  &
        1.0e0,                  &
        2.0e0,                  &
        6.0e0,                  &
        24.0e0,                 &
        120.0e0,                &
        720.0e0,                &
        5040.0e0,               &
        40320.0e0,              &
        362880.0e0,             &
        3628800.0e0,            &
        39916800.0e0,           &
        479001600.0e0,          &
        6227020800.0e0,         &
        87178291200.0e0,        &
        1307674368000.0e0,      &
        20922789888000.0e0,     &
        355687428096000.0e0,    &
        6402373705728000.0e0,   &
         .12164510040883200e18, &
         .24329020081766400e19, &
         .51090942171709440e20, &
         .11240007277776077e22, &
         .25852016738884977e23, &
         .62044840173323944e24, &
         .15511210043330986e26  ]

    integer, save   :: nmax   = 0

    real            :: x, xmin, xmax

    !
    ! Initialisation
    !
    if ( nmax == 0 ) then
        call r9gaml (xmin, xmax)
        nmax = xmax - 1.
    endif

    if ( n < 0 ) then
        fac = ieee_value( fac, ieee_quiet_nan )
    else
        if ( n <= 25 ) then
            fac = facn(n+1)

        else
            if ( n > nmax ) then
                fac = ieee_value( fac, ieee_positive_inf )

            else
                x = n + 1
                fac = exp ( (x-0.5)*log(x) - x + sq2pil + r9lgmc(x) )
            endif
        endif
    endif
end function fac

end module fullerton_poch
