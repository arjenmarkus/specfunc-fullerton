! fullerton_beta.f90 --
!     Module for the Beta function and several related functions
!
module fullerton_beta
    use ieee_arithmetic
    use fullerton_aux
    use fullerton_misc
    use fullerton_gamma_aux
    use fullerton_gamma

    implicit none

    interface beta_inc
        module procedure betai
    end interface

    interface log_beta
        module procedure albeta
    end interface

    private
    public :: beta, beta_inc, log_beta

contains

! beta --
!     Original:
!     june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function beta (a, b)
    real, intent(in) :: a, b

    real, save :: xmin   = 0.0
    real, save :: xmax   = 0.0
    real, save :: alnsml = 0.0

    !
    ! Initialisation
    !
    if ( xmin == 0.0 ) then
        call r9gaml (xmin, xmax)
        alnsml = log(r1mach(1))
    endif

    !
    ! Calculation
    !
    if ( a <= 0.0 .or. b <= 0.0 ) then
        beta = ieee_value( a, ieee_quiet_nan )
    else
        if ( a+b < xmax ) then
            beta = fgamma(a) * fgamma(b) / fgamma(a+b)

        else
            beta = albeta (a, b)
            beta = exp (beta)
        endif
    endif
end function beta

! betai --
!     Original:
!     august 1980 version.  w. fullerton, c3, los alamos scientific lab.
!     based on bosten and battiste, remark on algorithm 179, comm. acm,
!     v 17, p 153, (1974).
!
!                 input arguments --
!     x      upper limit of integration.  x must be in (0,1) inclusive.
!     p      first beta distribution parameter.  p must be gt 0.0.
!     q      second beta distribution parameter.  q must be gt 0.0.
!     betai  the incomplete beta function ratio is the probability that a
!            random variable from a beta distribution having parameters
!            p and q will be less than or equal to x.
!
real function betai (x, pin, qin)
    real, intent(in) :: x, pin, qin

    real, save :: eps    = 0.0
    real, save :: alneps = 0.0
    real, save :: sml    = 0.0
    real, save :: alnsml = 0.0

    real       :: y, p, q, p1, ps, finsum, term, xb, c
    integer    :: i, n, ib
    logical    :: yflip

    !
    !  Initialisation
    !
    if ( eps == 0.0 ) then
        eps    = r1mach(3)
        alneps = log(eps)
        sml    = r1mach(1)
        alnsml = log(sml)
    endif

    !
    ! Check arguments
    !
    if ( x < 0.0 .or. x > 1.0 ) then
        betai = ieee_value( x, ieee_quiet_nan )
        return
    endif
    if ( pin <= 0.0 .or. qin <= 0.0 ) then
        betai = ieee_value( x, ieee_quiet_nan )
        return
    endif

    !
    ! Calculation
    !
    y = x
    p = pin
    q = qin


    yflip = .true.
    if ( q <= p .and. x < 0.8 ) yflip = .false.
    if ( x < 0.2 ) yflip =.false.

    if ( yflip ) then
        y = 1.0 - y
        p = qin
        q = pin
    endif

    if ( (p+q)*y/(p+1.0) < eps ) then
        betai = 0.0
        xb = p*log(max(y,sml)) - log(p) - albeta(p,q)
        if ( xb > alnsml .and. y /= 0.0 ) then
            betai = exp (xb)
        endif
        if ( y /= x .or. p /= pin ) then
            betai = 1.0 - betai
        endif


    else
        !
        ! evaluate the infinite sum first.
        ! term will equal y**p/beta(ps,p) * (1.-ps)i * y**i / fac(i)
        !
        ps = q - aint(q)
        if ( ps == 0.0 ) then
            ps = 1.0
        endif

        xb    = p*alog(y) -  albeta(ps, p) - alog(p)
        betai = 0.0


        if ( xb >= alnsml ) then
            betai = exp (xb)
            term  = betai*p

            if ( ps /= 1.0 ) then
                n = max (alneps/log(y), 4.0)
                do i = 1,n
                    term = term*(float(i)-ps)*y/float(i)
                    betai = betai + term/(p+float(i))
                enddo
            endif
        endif

        !
        ! now evaluate the finite sum, maybe.
        !

        if ( q > 1.0 ) then
            xb     = p*alog(y) + q*log(1.0-y) - albeta(p,q) - log(q)
            ib     = max (xb/alnsml, 0.0)
            term   = exp (xb - float(ib)*alnsml)
            c      = 1.0/(1.0-y)
            p1     = q*c/(p+q-1.)

            finsum = 0.0

            n      = q
            if ( q == float(n) ) then
                n = n - 1
            endif

            do i = 1,n
                if ( p1 <= 1.0 .and. term/eps <= finsum ) then
                    exit
                endif
                term = (q-float(i-1))*c*term/(p+q-float(i))

                if ( term > 1.0 ) then
                    ib   = ib - 1
                    term = term*sml
                endif

                if ( ib == 0 ) then
                    finsum = finsum + term
                endif
            enddo

            betai = betai + finsum
        endif

        if ( y /= x .or. p /= pin ) then
            betai = 1.0 - betai
        endif

        betai = max (min (betai, 1.0), 0.0)
    endif

end function betai

! albeta --
!     Original:
!     july 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function albeta (a, b)
    real, intent(in) :: a, b

    real             :: p, q, corr

    !
    ! Calculation
    !
    p = min (a, b)
    q = max (a, b)

    if ( p <= 0.0 ) then
        albeta = ieee_value( p, ieee_quiet_nan )
    endif

    if ( p < 10.0 ) then
        if ( q < 10.0 ) then
            !
            ! p and q are small.
            !
            albeta = log(fgamma(p) * (fgamma(q)/fgamma(p+q)) )
        else
            !
            ! p is small, but q is big.
            !
            corr = r9lgmc(q) - r9lgmc(p+q)
            albeta = alngam(p) + corr + p - p*log(p+q) + (q-0.5)*lnrel(-p/(p+q))
        endif

    else
        !
        ! p and q are big.
        !
        corr = r9lgmc(p) + r9lgmc(q) - r9lgmc(p+q)
        albeta = -0.5*log(q) + sq2pil + corr + (p-0.5)*log(p/(p+q)) + q*lnrel(-p/(p+q))
    endif
end function albeta

end module fullerton_beta
