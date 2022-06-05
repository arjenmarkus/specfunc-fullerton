! fullerton_bessel.f90 --
!     Module with functions for modified Bessel functions -I0, I1, K0 and K1
!
!     TODO:
!     Provide functions for In and Kn via recursive relations
!
module fullerton_bessel
    use fullerton_gamma
    use fullerton_aux
    use ieee_arithmetic

    implicit none

    private
    public :: bessel_i0, bessel_i1, bessel_in, bessel_k0, bessel_k1, bessel_kn, bessel_knu, &
              bessel_j0_prime, bessel_j1_prime, bessel_jn_prime,                            &
              bessel_y0_prime, bessel_y1_prime, bessel_yn_prime,                            &
              bessel_i0_prime, bessel_i1_prime, bessel_in_prime,                            &
              bessel_k0_prime, bessel_k1_prime, bessel_kn_prime

    !
    ! Rename the functions (combine single and double precision versions)
    !
    interface bessel_i0
        module procedure besi0
    end interface
    interface bessel_i1
        module procedure besi1
    end interface
    interface bessel_k0
        module procedure besk0
    end interface
    interface bessel_k1
        module procedure besk1
    end interface
    interface bessel_knu
        module procedure besks
    end interface

contains

! besi0 --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     series for bi0        on the interval  0.          to  9.00000d+00
!                                            with weighted error   2.46e-18
!                                             log weighted error  17.61
!                                   significant figures required  17.90
!                                        decimal places required  18.15
!
real function besi0 (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

    real, save :: bi0cs(12) = [    &
        -.07660547252839144951e0, &
        1.927337953993808270e0,   &
         .2282644586920301339e0,  &
         .01304891466707290428e0, &
         .00043442709008164874e0, &
         .00000942265768600193e0, &
         .00000014340062895106e0, &
         .00000000161384906966e0, &
         .00000000001396650044e0, &
         .00000000000009579451e0, &
         .00000000000000053339e0, &
         .00000000000000000245e0  ]

    integer, save :: nti0 = 0
    real, save    :: xsml = 0.0
    real, save    :: xmax = 0.0

    real          :: y

    if ( nti0 == 0 ) then
        nti0 = inits (bi0cs, 12, 0.1*r1mach(3))
        xsml = sqrt (4.0*r1mach(3))
        xmax = log (r1mach(2))
    endif

    !
    ! Absolute value |x| of the argument
    !
    y = abs(x)

    !
    ! Argument |x| <= 3.0 or > 3.0
    !
    if ( y <= 3.0 ) then
        besi0 = 1.0
        if (y > xsml) besi0 = 2.75 + csevl (y*y/4.5-1.0, bi0cs, nti0)
    else
        if (y > xmax) then
            besi0 = ieee_value( x, ieee_positive_inf )
        else
            besi0 = exp(y) * besi0e(x)
        endif
    endif
end function besi0

! besi0e --
!     Original:
!     july 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!
real function besi0e (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

! series for bi0        on the interval  0.          to  9.00000d+00
!                                        with weighted error   2.46e-18
!                                         log weighted error  17.61
!                               significant figures required  17.90
!                                    decimal places required  18.15
!
    real, save :: bi0cs(12) = [   &
        -.07660547252839144951e0, &
        1.927337953993808270e0,   &
         .2282644586920301339e0,  &
         .01304891466707290428e0, &
         .00043442709008164874e0, &
         .00000942265768600193e0, &
         .00000014340062895106e0, &
         .00000000161384906966e0, &
         .00000000001396650044e0, &
         .00000000000009579451e0, &
         .00000000000000053339e0, &
         .00000000000000000245e0  ]

! series for ai0        on the interval  1.25000d-01 to  3.33333d-01
!                                        with weighted error   7.87e-17
!                                         log weighted error  16.10
!                               significant figures required  14.69
!                                    decimal places required  16.76
!
    real, save :: ai0cs(21) = [ &
         .07575994494023796e0,  &
         .00759138081082334e0,  &
         .00041531313389237e0,  &
         .00001070076463439e0,  &
        -.00000790117997921e0,  &
        -.00000078261435014e0,  &
         .00000027838499429e0,  &
         .00000000825247260e0,  &
        -.00000001204463945e0,  &
         .00000000155964859e0,  &
         .00000000022925563e0,  &
        -.00000000011916228e0,  &
         .00000000001757854e0,  &
         .00000000000112822e0,  &
        -.00000000000114684e0,  &
         .00000000000027155e0,  &
        -.00000000000002415e0,  &
        -.00000000000000608e0,  &
         .00000000000000314e0,  &
        -.00000000000000071e0,  &
         .00000000000000007e0   ]

! series for ai02       on the interval  0.          to  1.25000d-01
!                                        with weighted error   3.79e-17
!                                         log weighted error  16.42
!                               significant figures required  14.86
!                                    decimal places required  17.09
!
    real, save :: ai02cs(22) = [ &
         .05449041101410882e0,   &
         .00336911647825569e0,   &
         .00006889758346918e0,   &
         .00000289137052082e0,   &
         .00000020489185893e0,   &
         .00000002266668991e0,   &
         .00000000339623203e0,   &
         .00000000049406022e0,   &
         .00000000001188914e0,   &
        -.00000000003149915e0,   &
        -.00000000001321580e0,   &
        -.00000000000179419e0,   &
         .00000000000071801e0,   &
         .00000000000038529e0,   &
         .00000000000001539e0,   &
        -.00000000000004151e0,   &
        -.00000000000000954e0,   &
         .00000000000000382e0,   &
         .00000000000000176e0,   &
        -.00000000000000034e0,   &
        -.00000000000000027e0,   &
         .00000000000000003e0    ]

    integer, save :: nti0   = 0
    integer, save :: ntai0  = 0
    integer, save :: ntai02 = 0
    real, save    :: xsml   = 0.0

    real          :: y

    !
    ! Initialisation
    !
    if ( nti0 == 0 ) then
        nti0   = inits (bi0cs, 12, 0.1*r1mach(3))
        ntai0  = inits (ai0cs, 21, 0.1*r1mach(3))
        ntai02 = inits (ai02cs, 22, 0.1*r1mach(3))
        xsml   = sqrt (4.0*r1mach(3))
    endif

    !
    ! Absolute value |x| of the argument
    !
    y = abs(x)

    !
    ! Argument |x| <= 3.0 or > 3.0
    !
    if ( y <= 3.0 ) then
        besi0e = 1.0
        if (y > xsml) besi0e = exp(-y) * ( 2.75 + csevl (y*y/4.5-1.0, bi0cs, nti0) )

    elseif ( y <= 8.0 ) then
        besi0e = (.375 + csevl ((48./y-11.)/5., ai0cs, ntai0) ) / sqrt(y)

    else
        besi0e = (.375 + csevl (16./y-1., ai02cs, ntai02)) / sqrt(y)
    endif

end function besi0e

! besi1 --
!     Original:
!     oct 1983 version.  w. fullerton, c3, los alamos scientific lab.
!
real function besi1 (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bi1        on the interval  0.          to  9.00000d+00
!                                        with weighted error   2.40e-17
!                                         log weighted error  16.62
!                               significant figures required  16.23
!                                    decimal places required  17.14
!
    real, save :: bi1cs(11) = [ &
        -.001971713261099859e0, &
         .40734887667546481e0,  &
         .034838994299959456e0, &
         .001545394556300123e0, &
         .000041888521098377e0, &
         .000000764902676483e0, &
         .000000010042493924e0, &
         .000000000099322077e0, &
         .000000000000766380e0, &
         .000000000000004741e0, &
         .000000000000000024e0  ]

    integer, save :: nti1 = 0
    real, save    :: xmin = 0.0
    real, save    :: xsml = 0.0
    real, save    :: xmax = 0.0

    real          :: y

    !
    ! Initialisation
    !
    if ( nti1 == 0 ) then
        nti1 = inits (bi1cs, 11, 0.1*r1mach(3))

        xmin = 2.0*r1mach(1)
        xsml = sqrt (8.0*r1mach(3))
        xmax = log (r1mach(2))
    endif

    !
    ! Argument ranges:
    !     abs(x) < 3.0
    !     3.0 <= abs(x) <= 8.0
    !     abs(x) > 8.0
    !
    y = abs(x)

    if ( y <= 3.0 ) then
        besi1 = 0.0

        if ( y > xmin) then
            if ( y <= xsml ) then
                besi1 = 0.5*x
            else
                besi1 = x * (.875 + csevl(y*y/4.5-1., bi1cs, nti1))
            endif
        endif

    elseif ( y < xmax ) then
        besi1 = exp(y) * besi1e(x)

    else
        besi1 = ieee_value( x, ieee_positive_inf )
    endif
end function besi1

! besi1e --
!     Original:
!     oct 1983 version.  w. fullerton, c3, los alamos scientific lab.

real function besi1e (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bi1        on the interval  0.          to  9.00000d+00
!                                        with weighted error   2.40e-17
!                                         log weighted error  16.62
!                               significant figures required  16.23
!                                    decimal places required  17.14
!
    real, save :: bi1cs(11) = [ &
        -.001971713261099859e0, &
         .40734887667546481e0,  &
         .034838994299959456e0, &
         .001545394556300123e0, &
         .000041888521098377e0, &
         .000000764902676483e0, &
         .000000010042493924e0, &
         .000000000099322077e0, &
         .000000000000766380e0, &
         .000000000000004741e0, &
         .000000000000000024e0  ]
!
! series for ai1        on the interval  1.25000d-01 to  3.33333d-01
!                                        with weighted error   6.98e-17
!                                         log weighted error  16.16
!                               significant figures required  14.53
!                                    decimal places required  16.82
!
    real, save :: ai1cs(21) = [ &
        -.02846744181881479e0,  &
        -.01922953231443221e0,  &
        -.00061151858579437e0,  &
        -.00002069971253350e0,  &
         .00000858561914581e0,  &
         .00000104949824671e0,  &
        -.00000029183389184e0,  &
        -.00000001559378146e0,  &
         .00000001318012367e0,  &
        -.00000000144842341e0,  &
        -.00000000029085122e0,  &
         .00000000012663889e0,  &
        -.00000000001664947e0,  &
        -.00000000000166665e0,  &
         .00000000000124260e0,  &
        -.00000000000027315e0,  &
         .00000000000002023e0,  &
         .00000000000000730e0,  &
        -.00000000000000333e0,  &
         .00000000000000071e0,  &
        -.00000000000000006e0   ]
!
! series for ai12       on the interval  0.          to  1.25000d-01
!                                        with weighted error   3.55e-17
!                                         log weighted error  16.45
!                               significant figures required  14.69
!                                    decimal places required  17.12
!
    real, save :: ai12cs(22) = [ &
         .02857623501828014e0,   &
        -.00976109749136147e0,   &
        -.00011058893876263e0,   &
        -.00000388256480887e0,   &
        -.00000025122362377e0,   &
        -.00000002631468847e0,   &
        -.00000000383538039e0,   &
        -.00000000055897433e0,   &
        -.00000000001897495e0,   &
         .00000000003252602e0,   &
         .00000000001412580e0,   &
         .00000000000203564e0,   &
        -.00000000000071985e0,   &
        -.00000000000040836e0,   &
        -.00000000000002101e0,   &
         .00000000000004273e0,   &
         .00000000000001041e0,   &
        -.00000000000000382e0,   &
        -.00000000000000186e0,   &
         .00000000000000033e0,   &
         .00000000000000028e0,   &
        -.00000000000000003e0    ]

    integer, save :: nti1   = 0
    integer, save :: ntai1  = 0
    integer, save :: ntai12 = 0
    real, save    :: xmin   = 0.0
    real, save    :: xsml   = 0.0

    real          :: y

    !
    ! Initialisation
    !
    if ( nti1 == 0 ) then
        nti1   = inits (bi1cs, 11, 0.1*r1mach(3))
        ntai1  = inits (ai1cs, 21, 0.1*r1mach(3))
        ntai12 = inits (ai12cs, 22, 0.1*r1mach(3))

        xmin   = 2.0*r1mach(1)
        xsml   = sqrt (8.0*r1mach(3))
    endif

    !
    ! Argument ranges:
    !     abs(x) < 3.0
    !     3.0 <= abs(x) <= 8.0
    !     abs(x) > 8.0
    !
    y = abs(x)
    if ( y <= 3.0 ) then
        besi1e = 0.0

        if ( y > xmin ) then
            besi1e = 0.5*x
            if ( y > xsml ) then
                besi1e = x * (.875 + csevl(y*y/4.5-1., bi1cs,nti1))
            endif
            besi1e = exp(-y) * besi1e
        endif
    else
        if ( y <= 8.0 ) then
            besi1e = (.375 + csevl ((48./y-11.)/5., ai1cs, ntai1)) / sqrt(y)
        else
            besi1e = (.375 + csevl (16./y-1.0, ai12cs, ntai12)) / sqrt(y)
        endif

        besi1e = sign (besi1e, x)
   endif
end function besi1e

! besk0 --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
!     series for bk0        on the interval  0.          to  4.00000d+00
!                                            with weighted error   3.57e-19
!                                             log weighted error  18.45
!                                   significant figures required  17.99
!                                        decimal places required  18.97
!
real function besk0 (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

    real       :: y

    real, save :: bk0cs(11) = [      &
        -.03532739323390276872e0, &
         .3442898999246284869e0,  &
         .03597993651536150163e0, &
         .00126461541144692592e0, &
         .00002286212103119451e0, &
         .00000025347910790261e0, &
         .00000000190451637722e0, &
         .00000000001034969525e0, &
         .00000000000004259816e0, &
         .00000000000000013744e0, &
         .00000000000000000035e0  ]

    integer, save :: ntk0 = 0
    real, save    :: xsml = 0.0
    real, save    :: xmax = 0.0

    !
    ! Initialisation
    !
    if ( ntk0 == 0 ) then
        ntk0 = inits (bk0cs, 11, 0.1*r1mach(3))
        xsml = sqrt (4.0*r1mach(3))
        xmax = -log(r1mach(1))
        xmax = xmax - 0.5*xmax*log(xmax)/(xmax+0.5) - 0.01
    endif

    !
    ! Argument x must be positive
    !
    if ( x <= 0.0 ) then
        besk0 = ieee_value( x, ieee_quiet_nan )
        return
    endif

    !
    ! Argument x positive: x <= 2.0 and x > 2.0
    !
    if ( x <= 2.0 ) then
        y = 0.
        if ( x > xsml ) y = x*x
        besk0 = -log(0.5*x)*besi0(x) - .25 + csevl (.5*y-1., bk0cs, ntk0)
    else
        besk0 = 0.
        if ( x > xmax ) then
            return ! Note: Underflow - but is that bad?
        endif

        besk0 = exp(-x) * besk0e(x)
    endif
end function besk0

! besk0e --
!     Original:
!     july 1980 version.  w. fullerton, c3, los alamos scientific lab.
!
real function besk0e (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!   series for bk0        on the interval  0.          to  4.00000d+00
!                                          with weighted error   3.57e-19
!                                           log weighted error  18.45
!                                 significant figures required  17.99
!                                      decimal places required  18.97
    real, save :: bk0cs(11) = [   &
        -.03532739323390276872e0, &
         .3442898999246284869e0,  &
         .03597993651536150163e0, &
         .00126461541144692592e0, &
         .00002286212103119451e0, &
         .00000025347910790261e0, &
         .00000000190451637722e0, &
         .00000000001034969525e0, &
         .00000000000004259816e0, &
         .00000000000000013744e0, &
         .00000000000000000035e0  ]

! series for ak0        on the interval  1.25000d-01 to  5.00000d-01
!                                        with weighted error   5.34e-17
!                                         log weighted error  16.27
!                               significant figures required  14.92
!                                    decimal places required  16.89
!
    real, save :: ak0cs(17) = [   &
        -.07643947903327941e0,    &
        -.02235652605699819e0,    &
         .00077341811546938e0,    &
        -.00004281006688886e0,    &
         .00000308170017386e0,    &
        -.00000026393672220e0,    &
         .00000002563713036e0,    &
        -.00000000274270554e0,    &
         .00000000031694296e0,    &
        -.00000000003902353e0,    &
         .00000000000506804e0,    &
        -.00000000000068895e0,    &
         .00000000000009744e0,    &
        -.00000000000001427e0,    &
         .00000000000000215e0,    &
        -.00000000000000033e0,    &
         .00000000000000005e0     ]

! series for ak02       on the interval  0.          to  1.25000d-01
!                                        with weighted error   2.34e-17
!                                         log weighted error  16.63
!                               significant figures required  14.67
!                                    decimal places required  17.20
!
    real, save :: ak02cs(14) = [  &
        -.01201869826307592e0,    &
        -.00917485269102569e0,    &
         .00014445509317750e0,    &
        -.00000401361417543e0,    &
         .00000015678318108e0,    &
        -.00000000777011043e0,    &
         .00000000046111825e0,    &
        -.00000000003158592e0,    &
         .00000000000243501e0,    &
        -.00000000000020743e0,    &
         .00000000000001925e0,    &
        -.00000000000000192e0,    &
         .00000000000000020e0,    &
        -.00000000000000002e0     ]

    integer, save :: ntk0   = 0
    integer, save :: ntak0  = 0
    integer, save :: ntak02 = 0
    real, save    :: xsml   = 0.0

    real          :: y

    !
    ! Initialisation
    !
    if ( ntk0 ==0 ) then
        ntk0  = inits (bk0cs, 11, 0.1*r1mach(3))
        ntak0  = inits (ak0cs, 17, 0.1*r1mach(3))
        ntak02 = inits (ak02cs, 14, 0.1*r1mach(3))
        xsml   = sqrt (4.0*r1mach(3))
    endif

    !
    ! Argument x must be positive
    !
    if ( x <= 0.0 ) then
        besk0e = ieee_value( x, ieee_quiet_nan )
        return
    endif

    !
    ! Argument x positive: x <= 2.0 or x > 2.0 and x x < 8.0 or x > 8.0
    !
    if ( x <= 2.0 ) then
        y = 0.
        if ( x > xsml ) y = x*x
        besk0e = exp(x) * (-log(0.5*x)*besi0(x) - .25 + csevl (.5*y-1., bk0cs, ntk0) )

    elseif ( x <= 8.0 ) then
        besk0e = (1.25 + csevl ((16./x-5.)/3., ak0cs, ntak0)) / sqrt(x)

    else
        besk0e = (1.25 + csevl (16./x-1., ak02cs, ntak02)) / sqrt(x)

    endif
end function besk0e

! besk1 --
!     Original:
!     july 1980 version.  w. fullerton, c3, los alamos scientific lab.
!
real function besk1 (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bk1        on the interval  0.          to  4.00000d+00
!                                        with weighted error   7.02e-18
!                                         log weighted error  17.15
!                               significant figures required  16.73
!                                    decimal places required  17.67
!
    real, save :: bk1cs(11) = [  &
         .0253002273389477705e0, &
        -.353155960776544876e0,  &
        -.122611180822657148e0,  &
        -.0069757238596398643e0, &
        -.0001730288957513052e0, &
        -.0000024334061415659e0, &
        -.0000000221338763073e0, &
        -.0000000001411488392e0, &
        -.0000000000006666901e0, &
        -.0000000000000024274e0, &
        -.0000000000000000070e0  ]

    integer, save :: ntk1 = 0
    real, save    :: xmin = 0.0
    real, save    :: xsml = 0.0
    real, save    :: xmax = 0.0

    real          :: y
    !
    ! Initialisation
    !
    if ( ntk1 == 0 ) then
        ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))

        xmin = exp (max(log(r1mach(1)), -log(r1mach(2))) + .01)
        xsml = sqrt (4.0*r1mach(3))
        xmax = -log(r1mach(1))
        xmax = xmax - 0.5*xmax*alog(xmax)/(xmax+0.5) - 0.01
    endif

    !
    ! Argument ranges:
    !       x <= 0.0
    !     0.0 <  x <= 2.0
    !     2.0 <  x < xmax
    !       x > xmax
    !
    if ( x <= 0.0 ) then
        besk1 = ieee_value( x, ieee_quiet_nan )

    elseif ( x <= 2.0 ) then
        y = x*x
        besk1 = log(0.5*x)*besi1(x) + (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x

    else
        besk1 = 0.0

        if ( x <= xmax ) then
            besk1 = exp(-x) * besk1e(x)
        endif
    endif
end function besk1

! besk1e --
!     Original:
!     july 1980 version.  w. fullerton, c3, los alamos scientific lab.
!
real function besk1e (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bk1        on the interval  0.          to  4.00000d+00
!                                        with weighted error   7.02e-18
!                                         log weighted error  17.15
!                               significant figures required  16.73
!                                    decimal places required  17.67
!
    real, save :: bk1cs(11) = [  &
         .0253002273389477705e0, &
        -.353155960776544876e0,  &
        -.122611180822657148e0,  &
        -.0069757238596398643e0, &
        -.0001730288957513052e0, &
        -.0000024334061415659e0, &
        -.0000000221338763073e0, &
        -.0000000001411488392e0, &
        -.0000000000006666901e0, &
        -.0000000000000024274e0, &
        -.0000000000000000070e0  ]
!
! series for ak1        on the interval  1.25000d-01 to  5.00000d-01
!                                        with weighted error   6.06e-17
!                                         log weighted error  16.22
!                               significant figures required  15.41
!                                    decimal places required  16.83
!
    real, save :: ak1cs(17) = [ &
         .2744313406973883e0,   &
         .07571989953199368e0,  &
        -.00144105155647540e0,  &
         .00006650116955125e0,  &
        -.00000436998470952e0,  &
         .00000035402774997e0,  &
        -.00000003311163779e0,  &
         .00000000344597758e0,  &
        -.00000000038989323e0,  &
         .00000000004720819e0,  &
        -.00000000000604783e0,  &
         .00000000000081284e0,  &
        -.00000000000011386e0,  &
         .00000000000001654e0,  &
        -.00000000000000248e0,  &
         .00000000000000038e0,  &
        -.00000000000000006e0   ]
!
! series for ak12       on the interval  0.          to  1.25000d-01
!                                        with weighted error   2.58e-17
!                                         log weighted error  16.59
!                               significant figures required  15.22
!                                    decimal places required  17.16
!
    real, save :: ak12cs(14) = [ &
         .06379308343739001e0,   &
         .02832887813049721e0,   &
        -.00024753706739052e0,   &
         .00000577197245160e0,   &
        -.00000020689392195e0,   &
         .00000000973998344e0,   &
        -.00000000055853361e0,   &
         .00000000003732996e0,   &
        -.00000000000282505e0,   &
         .00000000000023720e0,   &
        -.00000000000002176e0,   &
         .00000000000000215e0,   &
        -.00000000000000022e0,   &
         .00000000000000002e0    ]

    integer, save :: ntk1   = 0
    integer, save :: ntak1  = 0
    integer, save :: ntak12 = 0
    real, save    :: xmin   = 0.0
    real, save    :: xsml   = 0.0

    real          :: y

    !
    ! Initialisation
    !
    if ( ntk1 == 0 ) then
        ntk1 = inits (bk1cs, 11, 0.1*r1mach(3))
        ntak1 = inits (ak1cs, 17, 0.1*r1mach(3))
        ntak12 = inits (ak12cs, 14, 0.1*r1mach(3))

        xmin = exp (max(log(r1mach(1)), -log(r1mach(2))) + .01)
        xsml = sqrt (4.0*r1mach(3))
    endif

    !
    ! Argument ranges:
    !       x <= 0.0
    !     0.0 <  x <= 2.0
    !     2.0 <  x <= 8.0
    !       x >  8.0
    !
    if ( x <= 0.0 ) then
        besk1e = ieee_value( x, ieee_quiet_nan )

    elseif( x <= 2.0 ) then
        if ( x < xmin ) then
            besk1e = ieee_value( x, ieee_positive_inf )

        else
            y = 0.
            if ( x > xsml ) y = x*x

            besk1e = exp(x) * (alog(0.5*x)*besi1(x) + (0.75 + csevl (.5*y-1., bk1cs, ntk1))/x )
        endif

    elseif ( x <= 8.0 ) then
        besk1e = (1.25 + csevl ((16./x-5.)/3., ak1cs, ntak1)) / sqrt(x)

    else
        besk1e = (1.25 + csevl (16./x-1., ak12cs, ntak12)) / sqrt(x)
    endif
end function besk1e

! bessel_kn --
!     Evaluate the Bessel function Kn(x) using a recursion relation
!
real function bessel_kn( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    real                :: knm2, knm1, kn
    integer             :: i

    if ( n < 0 .or. x < 0.0 ) then
        bessel_kn = ieee_value( x, ieee_quiet_nan )

    elseif ( x == 0.0 ) then
        bessel_kn = ieee_value( x, ieee_positive_inf )

    elseif ( n == 0 ) then
        bessel_kn = besk0(x)

    elseif ( n == 1 ) then
        bessel_kn = besk1(x)

    else
        knm2 = besk0(x)
        knm1 = besk1(x)

        do i = 2,n
            kn   = knm2 - 2.0 * real(i) * knm1 / x

            knm2 = knm1
            knm1 = kn
        enddo

        bessel_kn = kn
    endif
end function bessel_kn

! bessel_in --
!     Evaluate the Bessel function In(x) using backward recursion
!
!     Note:
!     Following the receipe from R.F. Halbsgewachs Sandia Laboratories, 1969 (report SC-M-69-336)
!
real function bessel_in( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    real                :: inp2, inp1, in, sumin
    integer             :: i, niter

    if ( n < 0 .or. x < 0.0 ) then
        bessel_in = ieee_value( x, ieee_quiet_nan )

    elseif ( x == 0.0 ) then
        bessel_in = ieee_value( x, ieee_positive_inf )

    elseif ( n == 0 ) then
        bessel_in = besi0(x)

    elseif ( n == 1 ) then
        bessel_in = besi1(x)

    else
        inp2  = 0.0
        inp1  = sqrt(tiny(x))
        sumin = 2.0 * inp1

        niter = max( n, int(2.0*abs(x)) ) + 51

        do i = niter,1
            in   = inp2 + 2.0 * real(i) * inp1 / x

            inp2 = inp1
            inp1 = in

            if ( i > 1 ) then
                sumin = sumin + 2.0 * in
            endif

            if ( i == n ) then
                bessel_in = in
            endif
        enddo

        sumin     = sumin + in   ! Add I0(x) for the final sum

        bessel_in = exp(x) * bessel_in / sumin
    endif
end function bessel_in

! besks --
!     Original:
!         july 1980 edition.   w. fullerton, c3, los alamos scientific lab.
!
subroutine besks (xnu, x, bk, downward)
    real, intent(in)              :: xnu, x
    real, intent(out)             :: bk(:)
    logical, intent(in), optional :: downward

    real, save :: xmax = 0.0
    real       :: expxi
    logical    :: down
    integer    :: i

    if ( xmax == 0.0 ) then
        xmax = -log (r1mach(1))
        xmax = xmax + 0.5*log(3.14*0.5/xmax)
    endif

    !
    ! Scale the results
    !

    down = .false.
    if ( present(downward) ) then
        down = downward
    endif

    call beskes (xnu, x, bk, down)

    expxi = exp(-x)

    do i = 1,size(bk)
        bk(i) = expxi * bk(i)
    enddo

end subroutine besks

! beskes --
!     Original:
!         july 1980 edition.   w. fullerton, c3, los alamos scientific lab.
!
subroutine beskes (xnu, x, bke, downward)
    real, intent(in)              :: xnu, x
    real, intent(out)             :: bke(:)
    logical, intent(in), optional :: downward

    real, save :: alnbig = 0.0
    real       :: v, vend, vincr
    real       :: bknu1
    integer    :: i, n
    logical    :: down, switch

    if ( alnbig == 0.0 ) then
        alnbig = log (r1mach(2))
    endif

    v = abs(xnu)
    n = size(bke)

    if ( v >= 1.0 .or. x <= 0.0 ) then
        bke = ieee_value( x, ieee_quiet_nan )
        return
    endif

    if ( n == 0 ) then
        return ! Nothing to do!
    endif

    down = .false.
    if ( present(downward) ) then
        down = downward
    endif

    call r9knus (v, x, bke(1), bknu1, switch)

    if ( n > 1 ) then
        vincr = sign(1.0, float(n))

        if ( switch .and. down ) then
            bke = ieee_value( x, ieee_quiet_nan )
            return
        endif

        bke(2) = bknu1

        if (down) then
            call r9knus (abs(xnu+vincr), x, bke(2), bknu1, switch)
        endif

        if ( n > 2) then
            vend = abs(xnu+float(n)) - 1.0
            v    = xnu
            do i = 3,n
                v = v + vincr
                bke(i) = 2.0*v*bke(i-1)/x + bke(i-2)
            enddo
        endif
    endif
end subroutine beskes

! r9knus --
!     Original:
!         june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
!
!         compute bessel functions exp(x) * k-sub-xnu (x)  and
!         exp(x) * k-sub-xnu+1 (x) for 0.0 .le. xnu .lt. 1.0 .
!
subroutine r9knus (xnu, x, bknu, bknu1, switch)
    real, intent(in)     :: xnu, x
    real, intent(out)    :: bknu, bknu1
    logical, intent(out) :: switch

!
! series for c0k        on the interval  0.          to  2.50000d-01
!                                        with weighted error   1.60e-17
!                                         log weighted error  16.79
!                               significant figures required  15.99
!                                    decimal places required  17.40
!
    real, save :: c0kcs(16) = [ &
        .060183057242626108e0,  &
       -.15364871433017286e0,   &
       -.011751176008210492e0,  &
       -.000852487888919795e0,  &
       -.000061329838767496e0,  &
       -.000004405228124551e0,  &
       -.000000316312467283e0,  &
       -.000000022710719382e0,  &
       -.000000001630564460e0,  &
       -.000000000117069392e0,  &
       -.000000000008405206e0,  &
       -.000000000000603466e0,  &
       -.000000000000043326e0,  &
       -.000000000000003110e0,  &
       -.000000000000000223e0,  &
       -.000000000000000016e0   ]

!
! series for znu1       on the interval -7.00000d-01 to  0.
!                                        with weighted error   1.43e-17
!                                         log weighted error  16.85
!                               significant figures required  16.08
!                                    decimal places required  17.38
!
    real, save :: znu1cs(12) = [ &
        .20330675699419173e0,    &
        .14007793341321977e0,    &
        .007916796961001613e0,   &
        .000339801182532104e0,   &
        .000011741975688989e0,   &
        .000000339357570612e0,   &
        .000000008425941769e0,   &
        .000000000183336677e0,   &
        .000000000003549698e0,   &
        .000000000000061903e0,   &
        .000000000000000981e0,   &
        .000000000000000014e0    ]

   real, parameter :: euler = 0.57721566490153286e0
   real, parameter :: sqpi2 = 1.2533141373155003e0
   real, parameter :: aln2  = 0.69314718055994531e0

   integer, save   :: ntc0k  = 0
   integer, save   :: ntznu1 = 0

   real, save      :: xnusml = 0.0
   real, save      :: xsml   = 0.0
   real, save      :: alnsml = 0.0
   real, save      :: alnbig = 0.0
   real, save      :: alneps = 0.0

   real            :: alpha(15)
   real            :: beta(15)
   real            :: a(15)
   real            :: a0, an, alnz, b0, bn, bknu0, bknud, c0, expx
   real            :: p1, p2, p3, qq, result, sqrtx, v, vlnz
   real            :: x2n, x2tov, xi, xmu, z, ztov
   integer         :: i, ii, n, nterms, inu

   if ( ntc0k == 0 ) then
       ntc0k = inits (c0kcs, 16, 0.1*r1mach(3))
       ntznu1 = inits (znu1cs, 12, 0.1*r1mach(3))

       xnusml = sqrt (r1mach(3)/8.0)
       xsml   = 0.1*r1mach(3)
       alnsml = log (r1mach(1))
       alnbig = log (r1mach(2))
       alneps = log (0.1*r1mach(3))
   endif

   !
   ! Hard condition: xnu between 0 and 1
   !
   if ( xnu < 0.0 .or. xnu >= 1.0 ) then
       error stop  ! Impossible case - this would mean a problem in the calling routine
   endif

   switch = .false.

   if ( x <= 2.0 ) then
        !
        ! x is small.  compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
        ! then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-.5,+.5)
        ! then to (0., .5), because k of negative order (-nu) = k of positive
        !
        v = xnu
        if ( xnu > 0.5 ) then
            v = 1.0 - xnu
        endif

        !
        ! carefully find (x/2)**xnu and z**xnu where z = x*x/4.
        !
        alnz = 2.0 * (log(x) - aln2)

        if ( x <= xnu ) then
            if ( -0.5*xnu*alnz-aln2-alog(xnu) > alnbig ) then
                bknu  = ieee_value( x, ieee_positive_inf )
                bknu1 = ieee_value( x, ieee_positive_inf )
            endif
        endif

        vlnz  = v*alnz
        x2tov = exp (0.5*vlnz)
        ztov  = 0.0

        if ( vlnz > alnsml ) then
            ztov = x2tov**2
        endif

        a0 = 0.5*fgamma(1.0+v)
        b0 = 0.5*fgamma(1.0-v)
        c0 = -euler

        if ( ztov > 0.5 .and. v > xnusml ) then
            c0 = -0.75 + csevl ((8.0*v)*v-1., c0kcs, ntc0k)
        endif

        if ( ztov <= 0.5 ) then
            alpha(1) = (a0-ztov*b0)/v
        else
            alpha(1) = c0 - alnz*(0.75 + csevl (vlnz/0.35+1.0, znu1cs, ntznu1))*b0
        endif

        beta(1) = -0.5*(a0+ztov*b0)

        z = 0.0
        if ( x > xsml ) then
            z = 0.25*x**2
        endif

        nterms = max (2.0, 11.0+(8.*alnz-25.19-alneps)/(4.28-alnz))

        do i = 2,nterms
            xi = i - 1
            a0 = a0/(xi*(xi-v))
            b0 = b0/(xi*(xi+v))
            alpha(i) = (alpha(i-1)+2.0*xi*a0)/(xi*(xi+v))
            beta(i)  = (xi-0.5*v)*alpha(i) - ztov*b0
        enddo

        bknu  = alpha(nterms)
        bknud = beta(nterms)

        do ii = 2,nterms
            i = nterms + 1 - ii
            bknu  = alpha(i) + bknu*z
            bknud = beta(i) + bknud*z
        enddo

        expx = exp(x)
        bknu = expx*bknu/x2tov

        if ( -0.5*(xnu+1.)*alnz-2.0*aln2 > alnbig ) then
            switch = .true.
            return
        endif

        bknud = expx*bknud*2.0/(x2tov*x)

        if ( xnu <= 0.5 ) then
            bknu1 = v*bknu/x - bknud
        else

            bknu0 = bknu
            bknu = -v*bknu/x - bknud
            bknu1 = 2.0*xnu*bknu/x + bknu0
        endif
    else
        !
        ! x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
        ! rational expansion.
        !
        sqrtx = sqrt(x)

        if ( x <= 1.0/xsml ) then
            an = -1.56 + 4.0/x
            bn = -0.29 - 0.22/x

            nterms = min(15, int(max(3., an+bn*alneps)))

            do inu = 1,2
                xmu = 0.
                if (inu == 1 )then
                    if ( xnu > xnusml ) then
                        xmu = (4.0*xnu)*xnu
                    endif
                else
                    xmu = 4.0*(abs(xnu)+1.)**2
                endif

                a(1) =  1.0 - xmu
                a(2) =  9.0 - xmu
                a(3) = 25.0 - xmu

                if ( a(2) == 0.) then
                    result = sqpi2*(16.*x+xmu+7.)/(16.*x*sqrtx)
                else
                    alpha(1) = 1.0
                    alpha(2) = (16.*x+a(2))/a(2)
                    alpha(3) = ((768.*x+48.*a(3))*x + a(2)*a(3))/(a(2)*a(3))

                    beta(1) = 1.0
                    beta(2) = (16.*x+(xmu+7.))/a(2)
                    beta(3) = ((768.*x+48.*(xmu+23.))*x + ((xmu+62.)*xmu+129.)) / (a(2)*a(3))

                    if ( nterms >= 4 ) then
                        do i = 4,nterms
                            n = i - 1
                            x2n = 2*n - 1

                            a(i) = (x2n+2.)**2 - xmu
                            qq = 16.*x2n/a(i)
                            p1 = -x2n*(float(12*n*n-20*n)-a(1))/((x2n-2.)*a(i)) - qq*x
                            p2 = (float(12*n*n-28*n+8)-a(1))/a(i) - qq*x
                            p3 = -x2n*a(i-3)/((x2n-2.)*a(i))

                            alpha(i) = -p1*alpha(i-1) - p2*alpha(i-2) - p3*alpha(i-3)
                            beta(i) = -p1*beta(i-1) - p2*beta(i-2) - p3*beta(i-3)
                        enddo
                    endif

                    result = sqpi2*beta(nterms)/(sqrtx*alpha(nterms))
                endif

                select case (inu)
                    case( 1 )
                        bknu = result
                    case( 2 )
                        bknu1 = result
                    case default
                        error stop ! Impossible case!
                end select
            enddo
        else

            !
            ! Very large x
            !
            bknu = sqpi2/sqrtx
            bknu1 = bknu
        endif
    endif
end subroutine r9knus

! bessel_*_prime --
!     Evaluate the first derivative of the various Bessel functions, using the
!     well-known relations
!
!     Note:
!     The code relies on the standard intrinsic functions for the ordinary Bessel functions
!
real function bessel_j0_prime( x )
    real, intent(in) :: x

    bessel_j0_prime = -besj1(x)

end function bessel_j0_prime

real function bessel_j1_prime( x )
    real, intent(in) :: x

    bessel_j1_prime = 0.5 * (besj0(x) - besjn(2,x))

end function bessel_j1_prime

real function bessel_jn_prime( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    bessel_jn_prime = 0.5 * (besjn(n-1,x) - besjn(n+1,x))

end function bessel_jn_prime

real function bessel_y0_prime( x )
    real, intent(in) :: x

    bessel_y0_prime = -besy1(x)

end function bessel_y0_prime

real function bessel_y1_prime( x )
    real, intent(in) :: x

    bessel_y1_prime = 0.5 * (besy0(x) - besyn(2,x))

end function bessel_y1_prime

real function bessel_yn_prime( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    bessel_yn_prime = 0.5 * (besyn(n-1,x) - besyn(n+1,x))

end function bessel_yn_prime

real function bessel_i0_prime( x )
    real, intent(in) :: x

    bessel_i0_prime = bessel_i1(x)

end function bessel_i0_prime

real function bessel_i1_prime( x )
    real, intent(in) :: x

    bessel_i1_prime = 0.5 * (bessel_i0(x) + bessel_in(2,x))

end function bessel_i1_prime

real function bessel_in_prime( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    bessel_in_prime = 0.5 * (bessel_in(n-1,x) + bessel_in(n+1,x))

end function bessel_in_prime

real function bessel_k0_prime( x )
    real, intent(in) :: x

    bessel_k0_prime = -bessel_k1(x)

end function bessel_k0_prime

real function bessel_k1_prime( x )
    real, intent(in) :: x

    bessel_k1_prime = -0.5 * (bessel_k0(x) + bessel_kn(2,x))

end function bessel_k1_prime

real function bessel_kn_prime( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    bessel_kn_prime = 0.5 * (bessel_kn(n-1,x) + bessel_kn(n+1,x))

end function bessel_kn_prime

end module fullerton_bessel
