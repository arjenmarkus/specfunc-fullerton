! fullerton_complex.f90 --
!     Module containing the complex versions of several sepcial functions
!
module fullerton_complex
     use ieee_arithmetic
     use fullerton_aux
     use fullerton_misc, only: lnrel, cbrt

     implicit none

     private
     public :: exprl, lnrel, cbrt, atan2, log10, fgamma, log_gamma, reciprocal_gamma, digamma, beta, log_beta

     interface exprl
         module procedure cexprl
     end interface

     interface lnrel
         module procedure clnrel
     end interface

     interface cbrt
         module procedure ccbrt
     end interface

     interface atan2
         module procedure catan2
     end interface

     interface log10
         module procedure clog10
     end interface

     interface fgamma
         module procedure cgamma ! In line with the alternative for the gamma function
     end interface

     interface log_gamma
         module procedure clngam
     end interface

     interface reciprocal_gamma
         module procedure cgamr
     end interface

     interface digamma
         module procedure cpsi
     end interface

     interface beta
         module procedure cbeta
     end interface

     interface log_beta
         module procedure clbeta
     end interface

contains

! cexprl --
!     Original:
!     august 1980 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate  (cexp(z)-1)/z .  for small cabs(z), we use the taylor
!     series.  we could instead use the expression
!            cexprl(z) = (exp(x)*exp(i*y)-1)/z
!                      = (x*exprel(x) * (1 - 2*sin(y/2)**2) - 2*sin(y/2)**2
!                                        + i*sin(y)*(1+x*exprel(x))) / z
!
complex function cexprl (z)
    complex, intent(in) :: z

    integer, save :: nterms = 0
    real, save    :: rbnd   = 0.0
    real, save    :: sqeps  = 0.0

    integer       :: i
    real          :: alneps, xln, xn, r

    !
    ! Initialisation
    !
    if ( nterms == 0 ) then
        alneps = log(r1mach(3))
        xn     = 3.72 - 0.3*alneps
        xln    = log((xn+1.0)/1.36)
        nterms = xn - (xn*xln+alneps)/(xln+1.36) + 1.5
        rbnd   = r1mach(3)
        sqeps  = 0.0
    endif

    !
    ! Calculation
    !
    r = abs(z)

    if (r > 0.5 ) then
        cexprl = (exp(z) - 1.0) / z

    else
        cexprl = (1.0, 0.0)

        if ( r > rbnd ) then
            cexprl = (0.0, 0.0)
            do i=1,nterms
                cexprl = 1.0 + cexprl*z/float(nterms+2-i)
            enddo
        endif
    endif
end function cexprl


! ccbrt --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
complex function ccbrt (z)
    complex, intent(in) :: z

    real :: r, theta

    theta = arg(z) / 3.0
    r     = cbrt (abs(z))

    ccbrt = cmplx (r*cos(theta), r*sin(theta))

end function ccbrt


! clnrel --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     clnrel(z) = clog(1+z) with relative error accuracy near z = 0.
!     let   rho = cabs(z)  and
!           r**2 = cabs(1+z)**2 = (1+x)**2 + y**2 = 1 + 2*x + rho**2 .
!     now if rho is small we may evaluate clnrel(z) accurately by
!           clog(1+z) = cmplx (alog(r), carg(1+z))
!                     = cmplx (0.5*alog(r**2), carg(1+z))
!                     = cmplx (0.5*alnrel(2*x+rho**2), carg(1+z))
!
complex function clnrel (z)
    complex, intent(in) :: z

    real, save :: sqeps = 0.0

    real       :: rho, x
    !
    ! Initialisation
    !
    if ( sqeps == 0.0 ) then
        sqeps = sqrt (r1mach(4))
    endif

    !
    ! Calculation
    !
    rho = abs(z)

    if ( rho <= 0.375 ) then
        clnrel = log (1.0+z)

    else
        x = real(z)
        clnrel = cmplx (0.5*lnrel(2.*x+rho**2), arg(1.0+z))
    endif
end function clnrel


! clog10 --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
complex function clog10 (z)
    complex, intent(in) :: z

    real, parameter :: aloge = 0.43429448190325182765e0

    clog10 = aloge * log(z)

end function clog10

! arg --
!     Original:
! april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
real function arg (z)
    complex, intent(in) :: z

    arg = 0.0
    if ( real(z) /= 0.0 .or. imag(z) /= 0.0 ) then
        arg = atan2 (imag(z), real(z))
    endif
end function arg


! catan2 --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
complex function catan2 (csn, ccs)
    complex, intent(in) :: csn, ccs

    real, parameter :: pi = 3.14159265358979323846e0

    if ( abs(ccs) /= 0.0 ) then
        catan2 = atan (csn/ccs)
        if ( real(ccs) < 0.0 ) then
            catan2 = catan2 + pi
        endif
        if ( real(catan2) > pi ) then
            catan2 = catan2 - 2.0*pi
        endif
    else
        if ( cabs(csn) /= 0.0 ) then
            catan2 = cmplx ( sign(0.5*pi,real(csn)), 0.0 )
        else
            catan2 = ieee_value( 0.0, ieee_quiet_nan )
        endif
   endif
end function catan2


! ccot --
!     Original:
!     march 1979 edition.  w. fullerton, c3, los alamos scientific lab.
!
complex function ccot (z)
    complex, intent(in) :: z

    real, save :: eps    = 0.0
    real, save :: xmax   = 0.0
    real, save :: ylarge = 0.0
    real, save :: ybig   = 0.0
    real, save :: rmin   = 0.0
    real, save :: ymin   = 1.5

    real       :: x, y, x2, y2, sn2x, den

    !
    ! Initialisation
    !
    if ( eps == 0.0 ) then
        eps    = r1mach(4)
        xmax   = 1.0/eps
        ylarge = -0.5*alog(0.5*r1mach(3))
        ybig   = -0.5*alog(0.5*sqrt(r1mach(3)))
        rmin   = exp (amax1(log(r1mach(1)), -log(r1mach(2))) + 0.01)
    endif

    !
    ! Calculation
    !
    x = real (z)
    y = imag (z)

    if ( abs(y) <= ylarge ) then
        x2 = 2.0*x
        y2 = 2.0*y

        if ( abs(x2) <= xmax ) then
            sn2x = sin(x2)
            den  = cosh(y2) - cos(x2)

            ccot = cmplx (sn2x/den, -sinh(y2)/den)

        else
            ccot = cmplx (0.0, 1.0/tanh(y2))
        endif
    else
        ccot = cmplx (0.0, sign (1.0, y))
    endif
end function ccot


! cbeta -
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
complex function cbeta (a, b)
    complex, intent(in) :: a, b

    real, save :: xmax = 0.0
    real       :: xmin

    !
    ! Initialisation
    !
    if ( xmax == 0.0 ) then
        call r9gaml (xmin, xmax)
    endif

    !
    ! Calculation
    !
    if ( real(a) <= 0.0 .or. real(b) <= 0.0 ) then
        cbeta = cmplx( ieee_value( 0.0, ieee_quiet_nan ), 0.0 )

    else
        if ( real(a)+real(b) < xmax ) then
            cbeta = cgamma(a) * ( cgamma(b) / cgamma(a+b) )

        else
            cbeta = cexp (clbeta(a, b))
        endif
    endif
end function cbeta


! clbeta --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!     a preliminary version that is portable, but not accurate enough.
!
complex function clbeta (a, b)
    complex, intent(in) :: a, b

    if ( real(a) <= 0.0 .or. real(b) <= 0.0 ) then
        clbeta = cmplx( ieee_value(0.0, ieee_quiet_nan), 0.0 )

    else
        clbeta = clngam(a) + clngam(b) - clngam(a+b)
    endif
end function clbeta

! cgamma --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!     a preliminary version that is portable, but not accurate enough.
!
complex function cgamma (z)
    complex, intent(in) :: z

    cgamma = exp (clngam(z))

end function cgamma


! cgamr --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!     this version is an inaccurate preliminary one.  eventually this
!     routine should be a fundamental routine with no dependence on cgamma.
!
complex function cgamr (z)
    complex, intent(in) :: z

    real :: x

    cgamr = (0.0, 0.0)
    x = real (z)

    if (x < 0.0 .and. aint(x) == x .and. imag(z) == 0.0 ) then
        return
    endif

    cgamr = clngam(z)
    cgamr = exp (-cgamr)
end function cgamr


! clngam --
!     Original:
!     august 1980 edition.  w. fullerton c3, los alamos scientific lab.
!     eventually clngam should make use of c8lgmc for all z except for
!     z in the vicinity of 1 and 2.
!
complex function clngam (zin)
    complex, intent(in) :: zin

    real, parameter :: pi     = 3.14159265358979324e0
    real, parameter :: sq2pil = 0.91893853320467274e0

    real, save :: bound = 0.0
    real, save :: dxrel = 0.0
    real, save :: rmax  = 0.0

    integer    :: i, n
    complex    :: z, corr
    real       :: x, y, argsum, cabsz

    !
    ! Initialisation
    !
    if ( bound == 0.0 ) then
        n     = -0.30*log(r1mach(3))
        ! bound = n*(0.1*eps)**(-1/(2*n-1))/(pi*exp(1))
        bound = 0.1171*float(n)*(0.1*r1mach(3))**(-1./(2.*float(n)-1.))
        dxrel = sqrt (r1mach(4))
        rmax  = r1mach(2)/log(r1mach(2))
    endif

    !
    ! Calculation
    !
    z = zin
    x = real(zin)
    y = imag(zin)

    corr  = (0.0, 0.0)
    cabsz = abs(z)
    if ( ( x >= 0.0 .and. cabsz > bound ) .or. &
         ( x <  0.0 .and. abs(y) > bound )     ) then
        !
        ! use stirling-s approximation for large z.
        !
        clngam = sq2pil + (z-0.5)*clog(z) - z + corr + c9lgmc(z)

    else
        if ( cabsz >= bound ) then
            !
            ! use the reflection formula for real(z) negative, cabs(z) large, and
            ! abs(aimag(y)) small.
            !
            if ( y > 0.0 ) then
                z = conjg (z)
            endif

            corr = exp (-cmplx(0.0,2.0*pi)*z)

            clngam = sq2pil + 1.0 - cmplx(0.0,pi)*(z-0.5) - clnrel(-corr) + (z-0.5)*clog(1.0-z) - z - c9lgmc(1.0-z)

            if ( y > 0.0 ) then
                clngam = conjg (clngam)
            endif

        else
            !
            ! use the recursion relation for cabs(z) small.
            !
            n      = sqrt (bound**2 - y**2) - x + 1.0
            argsum = 0.0
            corr   = (1.0, 0.0)

            do i=1,n
                argsum = argsum + arg(z)
                corr   = z*corr
                z      = 1.0 + z
            enddo

            corr = -cmplx (log(abs(corr)), argsum)

            !
            ! Now use Stirling's approximation again, with the calculated correction
            !
            clngam = sq2pil + (z-0.5)*clog(z) - z + corr + c9lgmc(z)
        endif
   endif
end function clngam


! cpsi --
!     Original:
!     may 1978 edition.  w. fullerton, c3, los alamos scientific lab.
!
complex function cpsi (zin)
    complex, intent(in) :: zin

    real, save :: bern(13) = [ &
        .83333333333333333e-1, &
       -.83333333333333333e-2, &
        .39682539682539683e-2, &
       -.41666666666666667e-2, &
        .75757575757575758e-2, &
       -.21092796092796093e-1, &
        .83333333333333333e-1, &
       -.44325980392156863e0,  &
        .30539543302701197e1,  &
       -.26456212121212121e2,  &
        .28146014492753623e3,  &
       -.34548853937728938e4,  &
        .54827583333333333e5   ]

    real, parameter :: pi = 3.141592653589793e0

    integer, save   :: nterm
    real, save      :: bound = 0.0
    real, save      :: dxrel = 0.0
    real, save      :: rmin  = 0.0
    real, save      :: rbig  = 0.0

    integer         :: i, n, ndx
    complex         :: z, z2inv, corr
    real            :: x, y, cabsz

    !
    ! Initialisation
    !
    if ( nterm == 0 ) then
        nterm = -0.30*log(r1mach(3))
        ! maybe bound = n*(0.1*eps)**(-1/(2*n-1)) / (pi*exp(1))
        bound = 0.1171*float(nterm) * (0.1*r1mach(3))**(-1.0/(2.0*float(nterm)-1.0))
        dxrel = sqrt(r1mach(4))
        rmin  = exp (amax1 (alog(r1mach(1)), -log(r1mach(2))) + 0.011 )
        rbig = 1.0/r1mach(3)
    endif

    !
    ! Calculation
    !
    z = zin
    x = real(z)
    y = imag(z)

    if ( y < 0.0 ) then
        z = conjg(z)
    endif

    corr = (0.0, 0.0)
    cabsz = abs(z)

    if ( ( x >= 0.0 .and. cabsz  > bound ) .or. &
         ( x <  0.0 .and. abs(y) > bound ) )  then

        !
        ! In this case, only the last part of the calculation is of interest
        !
    else
        if ( cabsz > bound ) then
            !
            ! use the reflection formula for real(z) negative, cabs(z) large, and
            ! abs(aimag(y)) small.
            !
            corr = -pi*ccot(pi*z)
            z = 1.0 - z
         else
            !
            ! use the recursion relation for cabs(z) small.
            !
            n = sqrt(bound**2-y**2) - x + 1.0
            do i=1,n
                corr = corr - 1.0/z
                z    = z + 1.0
            enddo
        endif
    endif

    !
    ! now evaluate the asymptotic series for suitably large z.
    !
    if ( cabsz > rbig ) then
        cpsi = log(z) + corr
    else
        cpsi = (0.0, 0.0)
        z2inv = 1.0/z**2
        do i=1,nterm
            ndx = nterm + 1 - i
            cpsi = bern(ndx) + z2inv*cpsi
        enddo
        cpsi = log(z) - 0.5/z - cpsi*z2inv + corr
    endif

    if ( y < 0.0 ) then
        cpsi = conjg(cpsi)
    endif
end function cpsi


! c9ln2r --
!     Original:
!     april 1978 edition.  w. fullerton c3, los alamos scientific lab.
!
!     evaluate  clog(1+z)  from 2-nd order with relative error accuracy so
!     that     clog(1+z) = z - z**2/2 + z**3*c9ln2r(z).
!
!     now  clog(1+z) = 0.5*alog(1+2*x+cabs(z)**2) + i*carg(1+z),
!     where x = real(z)  and  y = aimag(z).
!     we find
!         z**3 * c9ln2r(z) = -x*cabs(z)**2 - 0.25*cabs(z)**4
!            + (2*x+cabs(z)**2)**3 * r9ln2r(2*x+cabs(z)**2)
!            + i * (carg(1+z) + (x-1)*y)
!     the imaginary part must be evaluated carefully as
!         (atan(y/(1+x)) - y/(1+x)) + y/(1+x) - (1-x)*y
!           = (y/(1+x))**3 * r9atn1(y/(1+x)) + x**2*y/(1+x)
!
!     now we divide through by z**3 carefully.  write
!         1/z**3 = (x-i*y)/cabs(z)**3 * (1/cabs(z)**3)
!     then   c9ln2r(z) = ((x-i*y)/cabs(z))**3 * (-x/cabs(z) - cabs(z)/4
!            + 0.5*((2*x+cabs(z)**2)/cabs(z))**3 * r9ln2r(2*x+cabs(z)**2)
!            + i*y/(cabs(z)*(1+x)) * ((x/cabs(z))**2 +
!              + (y/(cabs(z)*(1+x)))**2 * r9atn1(y/(1+x)) ) )
!
!     if we let  xz = x/cabs(z)  and  yz = y/cabs(z)  we may write
!         c9ln2r(z) = (xz-i*yz)**3 * (-xz - cabs(z)/4
!            + 0.5*(2*xz+cabs(z))**3 * r9ln2r(2*x+cabs(z)**2)
!            + i*yz/(1+x) * (xz**2 + (yz/(1+x))**2*r9atn1(y/(1+x)) ))
!
complex function c9ln2r (z)
    complex, intent(In) :: z

    real :: x, y, xz, yz, cabsz, arg, rpart, y1x, aipart

    !
    ! Calculation
    !
    x = real (z)
    y = imag (z)

    cabsz = abs(z)

    if ( cabsz <= 0.8125 ) then
        c9ln2r = cmplx (1.0/3.0, 0.0)

        if ( cabsz /= 0.0) then
            xz    = x/cabsz
            yz    = y/cabsz

            arg   = 2.0*xz + cabsz
            rpart = 0.5*arg**3*r9ln2r(cabsz*arg) - xz - 0.25*cabsz
            y1x = yz/(1.0+x)
            aipart = y1x * (xz**2 + y1x**2*r9atn1(cabsz*y1x) )

            c9ln2r = cmplx(xz,-yz)**3 * cmplx(rpart,aipart)
        endif

    else
        c9ln2r = (log(1.0+z) - z*(1.0-0.5*z)) / z**3
    endif
end function c9ln2r


! c9lgmc --
!     Original:
!    april 1978 edition.  w. fullerton c3, los alamos scientific lab.
!
!    compute the log gamma correction term for large cabs(z) when real(z)
!    .ge. 0.0 and for large abs(aimag(y)) when real(z) .lt. 0.0.  we find
!    c9lgmc so that
!      clog(cgamma(z)) = 0.5*alog(2.*pi) + (z-0.5)*clog(z) - z + c9lgmc(z).
!
complex function c9lgmc (zin)
    complex, intent(in) :: zin

    real, save :: bern(11) = [    &
         .083333333333333333e0,   &
        -.0027777777777777778e0,  &
         .00079365079365079365e0, &
        -.00059523809523809524e0, &
         .00084175084175084175e0, &
        -.0019175269175269175e0,  &
         .0064102564102564103e0,  &
        -.029550653594771242e0,   &
         .17964437236883057e0,    &
       -1.3924322169059011e0,     &
       13.402864044168392e0       ]

    integer, save :: nterm = 0
    real, save    :: bound = 0.0
    real, save    :: xbig  = 0.0
    real, save    :: xmax  = 0.0

    integer       :: i, ndx
    complex       :: z, z2inv
    real          :: x, y, cabsz

    !
    ! Initialisation
    !
    if ( nterm == 0 ) then
        nterm = -0.30*log(r1mach(3))
        bound = 0.1170*float(nterm)* (0.1*r1mach(3))**(-1./(2.*float(nterm)-1.))
        xbig  = 1.0/sqrt(r1mach(3))
        xmax  = exp (min(log(r1mach(2)/12.0), -log(12.*r1mach(1))) )
    endif

    !
    ! Calculation
    !
    z     = zin
    x     = real (z)
    y     = aimag(z)
    cabsz = cabs(z)

    if ( ( x < 0.0 .and. abs(y) < bound ) .or. cabsz < bound ) then
        c9lgmc = cmplx( ieee_value( 0.0, ieee_quiet_nan ), 0.0 )

    else
        if ( cabsz < xmax ) then
            if ( cabsz > xbig ) then
                c9lgmc = 1.0/(12.0*z)

            else
                z2inv = 1.0/z**2
                c9lgmc = (0.0, 0.0)

                do i=1,nterm
                    ndx = nterm + 1 - i
                    c9lgmc = bern(ndx) + c9lgmc*z2inv
                enddo

                c9lgmc = c9lgmc/z
            endif
        else
            c9lgmc = (0.0, 0.0)
        endif
    endif
end function c9lgmc


! c0lgmc --
!     Original:
!     august 1980 edition.  w. fullerton c3, los alamos scientific lab.
!
!     evaluate  (z+0.5)*clog((z+1.0)/z) - 1.0  with relative error accuracy.
!     let q = 1.0/z so that
!         (z+0.5)*clog(1+1/z) - 1 = (z+0.5)*(clog(1+q) - q + q*q/2) - q*q/4
!            = (z+0.5)*q**3*c9ln2r(q) - q**2/4,
!     where  c9ln2r  is (clog(1+q) - q + 0.5*q**2) / q**3.
!
complex function c0lgmc (z)
    complex, intent(in) :: z

    real    :: cabsz
    complex :: q

    cabsz = abs(z)
    q     = 1.0/z

    if ( cabsz <= 1.23 ) then
        c0lgmc = (z+0.5)*log(1.0+q) - 1.0
    else
        c0lgmc = ((1.+.5*q)*c9ln2r(q) - .25) * q**2
    endif
end function c0lgmc


! r9atn1 --
!     Original:
!     april 1978 edition.  w. fullerton c3, los alamos scientific lab.
!
!     evaluate  atan(x)  from first order, that is, evaluate
!     (atan(x)-x)/x**3  with relative error accuracy so that
!            atan(x) = x + x**3*r9atn1(x).
!
real function r9atn1 (x)
    real, intent(in) :: x

!
! series for atn1       on the interval  0.          to  1.00000d+00
!                                        with weighted error   2.21e-17
!                                         log weighted error  16.66
!                               significant figures required  15.44
!                                    decimal places required  17.32
!
    real, save :: atn1cs(21) = [ &
        -.03283997535355202e0,   &
         .05833432343172412e0,   &
        -.00740036969671964e0,   &
         .00100978419933728e0,   &
        -.00014397871635652e0,   &
         .00002114512648992e0,   &
        -.00000317232107425e0,   &
         .00000048366203654e0,   &
        -.00000007467746546e0,   &
         .00000001164800896e0,   &
        -.00000000183208837e0,   &
         .00000000029019082e0,   &
        -.00000000004623885e0,   &
         .00000000000740552e0,   &
        -.00000000000119135e0,   &
         .00000000000019240e0,   &
        -.00000000000003118e0,   &
         .00000000000000506e0,   &
        -.00000000000000082e0,   &
         .00000000000000013e0,   &
        -.00000000000000002e0    ]

    integer, save :: ntatn1 = 0
    real, save    :: xsml   = 0.0
    real, save    :: xbig   = 0.0
    real, save    :: xmax   = 0.0

    real          :: eps, y

    !
    ! Initialisation
    !
    if ( ntatn1 == 0 ) then
        eps    = r1mach(3)
        ntatn1 = inits (atn1cs, 21, 0.1*eps)

        xsml   = sqrt (0.1*eps)
        xbig   = 1.571/sqrt(eps)
        xmax   = 1.571/eps
    endif

    !
    ! Calculation
    !
    y = abs(x)
    if ( y <= 1.0 ) then

        if ( y <= xsml ) then
            r9atn1 = -1.0/3.0

        else
            r9atn1 = -0.25 + csevl (2.0*y*y-1., atn1cs, ntatn1)
        endif

    else
        r9atn1 = (atan(x) - x) / x**3
    endif
end function r9atn1


! r9ln2r --
!     Original:
!     april 1978 edition.  w. fullerton c3, los alamos scientific lab.
!
!     evaluate  alog(1+x)  from 2-nd order with relative error accuracy so
!     that    alog(1+x) = x - x**2/2 + x**3*r9ln2r(x)
!
real function r9ln2r (x)
    real, intent(in) :: x

!
! series for ln21       on the interval -6.25000d-01 to  0.
!                                        with weighted error   2.49e-17
!                                         log weighted error  16.60
!                               significant figures required  15.87
!                                    decimal places required  17.31
!
    real, save :: ln21cs(26) = [ &
         .18111962513478810e0,   &
        -.15627123192872463e0,   &
         .028676305361557275e0,  &
        -.005558699655948139e0,  &
         .001117897665229983e0,  &
        -.000230805089823279e0,  &
         .000048598853341100e0,  &
        -.000010390127388903e0,  &
         .000002248456370739e0,  &
        -.000000491405927392e0,  &
         .000000108282565070e0,  &
        -.000000024025872763e0,  &
         .000000005362460047e0,  &
        -.000000001202995136e0,  &
         .000000000271078892e0,  &
        -.000000000061323562e0,  &
         .000000000013920858e0,  &
        -.000000000003169930e0,  &
         .000000000000723837e0,  &
        -.000000000000165700e0,  &
         .000000000000038018e0,  &
        -.000000000000008741e0,  &
         .000000000000002013e0,  &
        -.000000000000000464e0,  &
         .000000000000000107e0,  &
        -.000000000000000024e0   ]
!
! series for ln22       on the interval  0.          to  8.12500d-01
!                                        with weighted error   1.42e-17
!                                         log weighted error  16.85
!                               significant figures required  15.95
!                                    decimal places required  17.50
!
    real, save :: ln22cs(20) = [ &
        -.22242532535020461e0,   &
        -.061047100108078624e0,  &
         .007427235009750394e0,  &
        -.000933501826163697e0,  &
         .000120049907687260e0,  &
        -.000015704722952820e0,  &
         .000002081874781051e0,  &
        -.000000278919557764e0,  &
         .000000037693558237e0,  &
        -.000000005130902896e0,  &
         .000000000702714117e0,  &
        -.000000000096748595e0,  &
         .000000000013381046e0,  &
        -.000000000001858102e0,  &
         .000000000000258929e0,  &
        -.000000000000036195e0,  &
         .000000000000005074e0,  &
        -.000000000000000713e0,  &
         .000000000000000100e0,  &
        -.000000000000000014e0   ]

    integer, save :: ntln21 = 0
    integer, save :: ntln22 = 0
    real, save    :: xmin   = 0.0
    real, save    :: xbig   = 0.0
    real, save    :: xmax   = 0.0

    real          :: eps, sqeps

    !
    ! Initialisation
    !
    if ( ntln21 == 0 ) then
        eps    = r1mach(3)
        ntln21 = inits (ln21cs, 26, 0.1*eps)
        ntln22 = inits (ln22cs, 20, 0.1*eps)

        xmin   = -1.0 + sqrt(r1mach(4))
        sqeps  = sqrt(eps)
        xmax   = 6.0/sqeps
        xmax   = xmax - (eps*xmax**2 - 2.0*log(xmax)) / (2.0*eps*xmax)
        xbig   = 4.0/sqrt(sqeps)
        xbig   = xbig - (sqeps*xbig**2 - 2.0*log(xbig)) / (2.*sqeps*xbig)
    endif

    !
    ! Calculation
    !
    if (x >= -0.625 .and. x <= 0.8125 ) then
        if ( x < 0.0 ) then
            r9ln2r = 0.375 + csevl (16.*x/5.+1.0, ln21cs, ntln21)

        else
            r9ln2r = 0.375 + csevl (32.*x/13.-1.0, ln22cs, ntln22)
        endif

    else
        r9ln2r = (log(1.0+x) - x*(1.0-0.5*x) ) / x**3
    endif
end function r9ln2r

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

end module fullerton_complex
