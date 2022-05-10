! fullerton_misc.f90 --
!     Module containing miscellaneous functions
!
module fullerton_misc
    use ieee_arithmetic
    use fullerton_aux

    implicit none

    private
    public :: cbrt, lnrel, exprel, integral_dawson, integral_spence, dawson_prime

    interface lnrel
        module procedure alnrel
    end interface

    interface dawson_prime
        module procedure dawsp
    end interface

    interface integral_dawson
        module procedure daws
    end interface

    interface integral_spence
        module procedure spenc
    end interface

contains

! cbrt --
!     Original:
!     june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function cbrt (x)
    real, intent(in) :: x

    real, save :: cbrt2(5) = [ &
       0.62996052494743658e0,  &
       0.79370052598409974e0,  &
       1.0e0,                  &
       1.25992104989487316e0,  &
       1.58740105196819947e0   ]

    integer, save :: niter = 0

    integer       :: irem, ixpnt, n, iter
    real          :: y, cbrtsq

    !
    ! Initialisation
    !
    if ( niter == 0 ) then
        niter = 1.443*log(-.106*log(0.1*r1mach(3))) + 1.0
    endif

    !
    ! Calculation
    !
    cbrt = 0.0

    if ( x /= 0.0 ) then
        call r9upak (abs(x), y, n)
        ixpnt = n/3
        irem  = n - 3*ixpnt + 3
        !
        ! the approximation below is a generalized chebyshev series converted
        ! to polynomial form.  the approx is nearly best in the sense of
        ! relative error with 4.085 digits accuracy.
        !
        cbrt = .439581e0 + y*(.928549e0 + y*(-.512653e0 + y*.144586e0))

        do iter=1,niter
           cbrtsq = cbrt*cbrt
           cbrt   = cbrt + (y-cbrt*cbrtsq)/(3.0*cbrtsq)
        enddo

        cbrt = r9pak (cbrt2(irem)*sign(cbrt,x), ixpnt)
    endif
end function cbrt

! alnrel --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
real function alnrel (x)
    real, intent(in) :: x

!
! series for alnr       on the interval -3.75000d-01 to  3.75000d-01
!                                        with weighted error   1.93e-17
!                                         log weighted error  16.72
!                               significant figures required  16.44
!                                    decimal places required  17.40
!
    real, save :: alnrcs(23) = [ &
        1.0378693562743770e0,    &
        -.13364301504908918e0,   &
         .019408249135520563e0,  &
        -.003010755112753577e0,  &
         .000486946147971548e0,  &
        -.000081054881893175e0,  &
         .000013778847799559e0,  &
        -.000002380221089435e0,  &
         .000000416404162138e0,  &
        -.000000073595828378e0,  &
         .000000013117611876e0,  &
        -.000000002354670931e0,  &
         .000000000425227732e0,  &
        -.000000000077190894e0,  &
         .000000000014075746e0,  &
        -.000000000002576907e0,  &
         .000000000000473424e0,  &
        -.000000000000087249e0,  &
         .000000000000016124e0,  &
        -.000000000000002987e0,  &
         .000000000000000554e0,  &
        -.000000000000000103e0,  &
         .000000000000000019e0   ]

    integer, save :: nlnrel = 0
    real, save    :: xmin   = 0.0

    !
    ! Initialisation
    !
    if ( nlnrel == 0 ) then
        nlnrel = inits (alnrcs, 23, 0.1*r1mach(3))
        xmin = -1.0 + sqrt(r1mach(4))
    endif

    !
    ! Argument ranges:
    !       x <= -1.0
    !    -1.0 <  x < 0.375
    !       x >  0.375
    !
    if ( x <= -1.0 ) then
        alnrel = ieee_value( x, ieee_quiet_nan )

    elseif ( abs(x) <= 0.375 ) then
        alnrel = x*(1. - x*csevl (x/.375, alnrcs, nlnrel))

    else
        alnrel = log (1.0+x)
   endif
end function alnrel

! exprel --
!     Original:
!     august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate  exprel(x) = (exp(x) - 1.0) / x.   for small abs(x) the
!     taylor series is used.  if x is negative the reflection formula
!             exprel(x) = exp(x) * exprel(abs(x))
!     may be used  this reflection formula will be of use when the
!     evaluation for small abs(x) is done by chebyshev series rather than
!     taylor series.
!
real function exprel (x)
    real, intent(in) :: x

    integer, save :: nterms = 0
    real, save    :: xbnd   = 0.0

    real          :: alneps, xn, xln, absx
    integer       :: i, ir, irold

    !
    ! Initialisation
    !
    if ( nterms == 0 ) then
        alneps = log(r1mach(3))
        xn     = 3.72 - 0.3*alneps
        xln    = log((xn+1.0)/1.36)
        nterms = xn - (xn*xln+alneps)/(xln+1.36) + 1.5
        xbnd   = r1mach(3)
    endif

    !
    ! Argument ranges:
    !     abs(x) < 0.5
    !     3.0 <= abs(x) <= 8.0
    !     abs(x) > 8.0
    !
    absx = abs(x)

    if ( absx < 0.5 ) then
        if ( absx < xbnd ) then
            exprel = 1.0
        else
            exprel = 0.0
            do i = 1,nterms
                exprel = 1.0 + exprel*x/float(nterms+2-i)
            enddo
        endif

    else
        exprel = (exp(x)-1.0)/x
    endif
end function exprel

! daws --
!     Original:
!     may 1980 edition.  w. fullerton c3, los alamos scientific lab.
!
real function daws (x)
    real, intent(in) :: x

!
! series for daw        on the interval  0.          to  1.00000d+00
!                                        with weighted error   3.83e-17
!                                         log weighted error  16.42
!                               significant figures required  15.78
!                                    decimal places required  16.97
!
    real, save :: dawcs(13) = [ &
        -.006351734375145949e0, &
        -.22940714796773869e0,  &
         .022130500939084764e0, &
        -.001549265453892985e0, &
         .000084973277156849e0, &
        -.000003828266270972e0, &
         .000000146285480625e0, &
        -.000000004851982381e0, &
         .000000000142146357e0, &
        -.000000000003728836e0, &
         .000000000000088549e0, &
        -.000000000000001920e0, &
         .000000000000000038e0  ]
!
! series for daw2       on the interval  0.          to  1.60000d+01
!                                        with weighted error   5.17e-17
!                                         log weighted error  16.29
!                               significant figures required  15.90
!                                    decimal places required  17.02
!
    real, save :: daw2cs(29) = [ &
        -.056886544105215527e0,  &
        -.31811346996168131e0,   &
         .20873845413642237e0,   &
        -.12475409913779131e0,   &
         .067869305186676777e0,  &
        -.033659144895270940e0,  &
         .015260781271987972e0,  &
        -.006348370962596214e0,  &
         .002432674092074852e0,  &
        -.000862195414910650e0,  &
         .000283765733363216e0,  &
        -.000087057549874170e0,  &
         .000024986849985481e0,  &
        -.000006731928676416e0,  &
         .000001707857878557e0,  &
        -.000000409175512264e0,  &
         .000000092828292216e0,  &
        -.000000019991403610e0,  &
         .000000004096349064e0,  &
        -.000000000800324095e0,  &
         .000000000149385031e0,  &
        -.000000000026687999e0,  &
         .000000000004571221e0,  &
        -.000000000000751873e0,  &
         .000000000000118931e0,  &
        -.000000000000018116e0,  &
         .000000000000002661e0,  &
        -.000000000000000377e0,  &
         .000000000000000051e0   ]
!
! series for dawa       on the interval  0.          to  6.25000d-02
!                                        with weighted error   2.24e-17
!                                         log weighted error  16.65
!                               significant figures required  14.73
!                                    decimal places required  17.36
!
    real, save :: dawacs(26) = [ &
         .01690485637765704e0,   &
         .00868325227840695e0,   &
         .00024248640424177e0,   &
         .00001261182399572e0,   &
         .00000106645331463e0,   &
         .00000013581597947e0,   &
         .00000002171042356e0,   &
         .00000000286701050e0,   &
        -.00000000019013363e0,   &
        -.00000000030977804e0,   &
        -.00000000010294148e0,   &
        -.00000000000626035e0,   &
         .00000000000856313e0,   &
         .00000000000303304e0,   &
        -.00000000000025236e0,   &
        -.00000000000042106e0,   &
        -.00000000000004431e0,   &
         .00000000000004911e0,   &
         .00000000000001235e0,   &
        -.00000000000000578e0,   &
        -.00000000000000228e0,   &
         .00000000000000076e0,   &
         .00000000000000038e0,   &
        -.00000000000000011e0,   &
        -.00000000000000006e0,   &
         .00000000000000002e0    ]

    integer, save :: ntdaw  = 0
    integer, save :: ntdaw2 = 0
    integer, save :: ntdawa = 0
    real, save    :: xsml   = 0.0
    real, save    :: xbig   = 0.0
    real, save    :: xmax   = 0.0

    real          :: eps, y

    !
    ! Initialisation
    !
    if ( ntdaw == 0 ) then
        eps    = r1mach(3)
        ntdaw  = inits (dawcs,  13, 0.1*eps)
        ntdaw2 = inits (daw2cs, 29, 0.1*eps)
        ntdawa = inits (dawacs, 26, 0.1*eps)

        xsml   = sqrt (1.5*eps)
        xbig   = sqrt (0.5/eps)
        xmax   = exp (min (-log(2.*r1mach(1)), log(r1mach(2))) - 0.01)
    endif

    !
    ! Argument ranges:
    !            abs(x) <= 1.0
    !     1.0 <  abs(x) <= 4.0
    !     4.0 <  abs(x) <= xmax
    !            abs(x) >  xmax
    !
    y = abs(x)
    if ( y <= 1.0 ) then
        daws = x
        if ( y > xsml ) then
            daws = x * (0.75 + csevl (2.0*y*y-1.0, dawcs, ntdaw))
        endif

    elseif ( y <= 4.0 ) then
      daws = x * (0.25 + csevl (0.125*y*y-1.0, daw2cs, ntdaw2))

    elseif ( y <= xmax ) then
        daws = 0.5/x
        if ( y < xbig ) then
           daws = (0.5 + csevl (32.0/y**2-1.0, dawacs, ntdawa)) / x
        endif

    else
      daws = 0.0
    endif
end function daws

! dawsp --
!     Calculate the first derivative of the Dawson integral
!     See:  https://mathworld.wolfram.com/DawsonsIntegral.html
!
real function dawsp (x)
    real, intent(in) :: x

    dawsp = 1.0 - 2.0 * x * daws(x)
end function dawsp

! spenc --
!     Original:
!     feb 1978 edition.   w. fullerton, c3, los alamos scientific lab.
!
!     evaluate a form of spence-s function defined by
!            integral from 0 to x of  -log(abs(1-y))/y  dy.
!     for abs(x).le.1, the uniformly convergent expansion
!            spenc = sum k=1,infinity  x**k / k**2     is valid.
!
!     spence-s function can be used to evaluate much more general integral
!     forms.  for example,
!           integral from 0 to z of  log(a*x+b)/(c*x+d)  dx  =
!                log(abs(b-a*d/c))*log(abs(a*(c*z+d)/(a*d-b*c)))/c
!                  - spenc (a*(c*z+d)/(a*d-b*c)) / c.
!
!     ref -- k. mitchell, philosophical magazine, 40, p.351 (1949).
!            stegun and abromowitz, ams 55, p.1004.
!
real function spenc (x)
    real, intent(in) :: x

!
! series for spen       on the interval  0.          to  5.00000d-01
!                                        with weighted error   6.82e-17
!                                         log weighted error  16.17
!                               significant figures required  15.22
!                                    decimal places required  16.81
!
    real, save :: spencs(19) = [ &
         .1527365598892406e0,    &
         .08169658058051014e0,   &
         .00581415714077873e0,   &
         .00053716198145415e0,   &
         .00005724704675185e0,   &
         .00000667454612164e0,   &
         .00000082764673397e0,   &
         .00000010733156730e0,   &
         .00000001440077294e0,   &
         .00000000198444202e0,   &
         .00000000027940058e0,   &
         .00000000004003991e0,   &
         .00000000000582346e0,   &
         .00000000000085767e0,   &
         .00000000000012768e0,   &
         .00000000000001918e0,   &
         .00000000000000290e0,   &
         .00000000000000044e0,   &
         .00000000000000006e0    ]
!
! pi26 = pi**2/6.0
    real, save :: pi26 = 1.644934066848226e0

    integer, save :: nspenc = 0
    real, save    :: xbig   = 0.0

    real          :: aln

    !
    ! Initialisation
    !
    if ( nspenc == 0 ) then
        nspenc = inits (spencs, 19, 0.1*r1mach(3))
        xbig = 1.0/r1mach(3)
    endif

    !
    ! Argument ranges: multiple
    !
    if ( x <= -1.0 ) then
        aln = log(1.0-x)
        spenc = -pi26 - 0.5*aln*(2.0*log(-x)-aln)
        if ( x > -xbig ) then
           spenc = spenc + (1.0 + csevl (4.0/(1.0-x)-1.0, spencs, nspenc)) / (1.0-x)
        endif

    elseif ( -1.0 < x .and. x < 0.0 ) then
        spenc = -0.5*log(1.0-x)**2 - x*(1.0 + csevl (4.0*x/(x-1.0)-1.0, spencs, nspenc)) / (x-1.0)

    elseif (  0.0 <= x .and. x <= 0.5 ) then
        spenc = x*(1.0 + csevl (4.0*x-1.0, spencs, nspenc))

    elseif (  0.5 < x .and. x <= 1.0 ) then
        spenc = pi26
        if ( x /= 1.0 ) then
            spenc = pi26 - log(x)*log(1.0-x) - (1.0-x)*(1.0 + csevl (4.0*(1.0-x)-1.0, spencs, nspenc))
        endif

    elseif (  1.0 < x .and. x <= 2.0 ) then
        spenc = pi26 - 0.5*log(x)*log((x-1.0)**2/x) + (x-1.)*(1.0 + csevl (4.0*(x-1.)/x-1.0, spencs, nspenc))/x

    else
        spenc = 2.0*pi26 - 0.5*log(x)**2

        if ( x < xbig ) then
            spenc = spenc - (1.0 + csevl (4.0/x-1.0, spencs, nspenc))/x
        endif
    endif
end function spenc

end module fullerton_misc
