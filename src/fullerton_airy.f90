! fullerton_airy.f90 --
!     Module for the Airy function Ai and Bi
!
module fullerton_airy
    use fullerton_aux
    use ieee_arithmetic

    implicit none

    private
    public :: airy_ai, airy_bi, airy_aiprime, airy_biprime

    interface airy_ai
        module procedure ai
    end interface

    interface airy_bi
        module procedure bi
    end interface

    interface airy_aiprime
        module procedure aid
    end interface

    interface airy_biprime
        module procedure bid
    end interface

    ! atr = 16.0/(sqrt(8.)-1.)  and  btr = -(sqrt(8.)+1.)/(sqrt(8.)-1.)
    real, parameter :: atr = 8.7506905708484345e0
    real, parameter :: btr = -2.093836321356054e0

contains

! ai --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
real function ai (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for aif        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.09e-19
!                                         log weighted error  18.96
!                               significant figures required  17.76
!                                    decimal places required  19.44
!
    real, save :: aifcs(9) = [ &
        -.03797135849666999750e0, &
         .05919188853726363857e0, &
         .00098629280577279975e0, &
         .00000684884381907656e0, &
         .00000002594202596219e0, &
         .00000000006176612774e0, &
         .00000000000010092454e0, &
         .00000000000000012014e0, &
         .00000000000000000010e0  ]
!
! series for aig        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.51e-17
!                                         log weighted error  16.82
!                               significant figures required  15.19
!                                    decimal places required  17.27
!
    real, save :: aigcs(8) = [ &
         .01815236558116127e0, &
         .02157256316601076e0, &
         .00025678356987483e0, &
         .00000142652141197e0, &
         .00000000457211492e0, &
         .00000000000952517e0, &
         .00000000000001392e0, &
         .00000000000000001e0  ]

    integer, save :: naif  = 0
    integer, save :: naig  = 0
    real, save    :: x3sml = 0.0
    real, save    :: xmax  = 0.0

    real          :: xm, theta, z

    !
    ! Initialisation
    !
    if (naif == 0) then
        naif = inits (aifcs, 9, 0.1*r1mach(3))
        naig = inits (aigcs, 8, 0.1*r1mach(3))

        x3sml = r1mach(3)**0.3334
        xmax = (-1.5*log(r1mach(1)))**0.6667
        xmax = xmax - xmax*log(xmax)/(4.0*xmax*sqrt(xmax)+1.0) - 0.01
    endif

    !
    ! Argument ranges:
    !     x < -1.0
    !     -1.0 <= x <= xmax
    !     x > xmax
    !
    if ( x < -1.0 ) then
        call r9aimp (x, xm, theta)
        ai = xm * cos(theta)

    elseif ( x <= 1.0 ) then
         z = 0.0
         if (abs(x) > x3sml) z = x**3
         ai = 0.375 + (csevl (z, aifcs, naif) - x*(0.25 + csevl (z, aigcs, naig)) )

    elseif ( x <= xmax ) then
          ai = aie(x) * exp(-2.0*x*sqrt(x)/3.0)

    else
          ai = 0.0
    endif
end function ai

! aid --
!     Original:
!    july 1980 edition.  w. fullerton, bell labs.
!
!    evaluate the derivative of the airy function ai(x).
!
real function aid (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for aif on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   5.22e-18
!                                         log weighted error  17.28
!                               significant figures required  16.01
!                                    decimal places required  17.73
!
    real, save :: aifcs(8) = [ &
         0.10527461226531408809e0, &
         0.01183613628152997844e0, &
         0.00012328104173225664e0, &
         0.00000062261225638140e0, &
         0.00000000185298887844e0, &
         0.00000000000363328873e0, &
         0.00000000000000504622e0, &
         0.00000000000000000522e0  ]
!
! series for aig on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   3.14e-19
!                                         log weighted error  18.50
!                               significant figures required  17.44
!                                    decimal places required  18.98
!
    real, save :: aigcs(9) = [ &
         0.021233878150918666852e0, &
         0.086315930335214406752e0, &
         0.001797594720383231358e0, &
         0.000014265499875550693e0, &
         0.000000059437995283683e0, &
         0.000000000152403366479e0, &
         0.000000000000264587660e0, &
         0.000000000000000331562e0, &
         0.000000000000000000314e0  ]

    integer, save :: naif  = 0
    integer, save :: naig  = 0
    real, save    :: x2sml = 0.0
    real, save    :: x3sml = 0.0

    real          :: phi, x2, x3, xn

    !
    ! Initialisation
    !
    if ( naif == 0 ) then
        naif = inits (aifcs, 8, 0.1*r1mach(3))
        naig = inits (aigcs, 9, 0.1*r1mach(3))

        x3sml = r1mach(3)**0.3334
        x2sml = sqrt(r1mach(3))
    endif

    !
    ! Argument ranges:
    !     x < -1.0
    !     -1 <= x <= 1.0
    !     x > xmax
    !
    if ( x <  -1.0 ) then
        call r9admp (x, xn, phi)
        aid = xn * cos(phi)

    elseif ( x <= 1.0 ) then
        x3 = 0.0
        if ( abs(x) > x3sml ) x3 = x**3
        x2 = 0.0
        if ( abs(x) > x2sml ) x2 = x*x
        aid = (x2*(0.125 + csevl (x3, aifcs, naif)) - csevl (x3, aigcs, naig)) - 0.25

    else
        aid = aide(x) * exp(-2.0*x*sqrt(x)/3.0)
    endif
end function aid

! aie --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!     evaluate ai(x) for x .le. 0.0 and  ai(x)*exp(zeta)  where
!     zeta = 2/3 * x**(3/2)  for x .ge. 0.0
!
real function aie (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for aif        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.09e-19
!                                         log weighted error  18.96
!                               significant figures required  17.76
!                                    decimal places required  19.44
!
    real, save :: aifcs(9) = [    &
        -.03797135849666999750e0, &
         .05919188853726363857e0, &
         .00098629280577279975e0, &
         .00000684884381907656e0, &
         .00000002594202596219e0, &
         .00000000006176612774e0, &
         .00000000000010092454e0, &
         .00000000000000012014e0, &
         .00000000000000000010e0  ]
!
! series for aig        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.51e-17
!                                         log weighted error  16.82
!                               significant figures required  15.19
!                                    decimal places required  17.27
!
    real, save :: aigcs(8) = [  &
         .01815236558116127e0, &
         .02157256316601076e0, &
         .00025678356987483e0, &
         .00000142652141197e0, &
         .00000000457211492e0, &
         .00000000000952517e0, &
         .00000000000001392e0, &
         .00000000000000001e0  ]
!
! series for aip        on the interval  0.          to  1.00000d+00
!                                        with weighted error   5.10e-17
!                                         log weighted error  16.29
!                               significant figures required  14.41
!                                    decimal places required  17.06
!
    real, save :: aipcs(34) = [ &
        -.0187519297793868e0,   &
        -.0091443848250055e0,   &
         .0009010457337825e0,   &
        -.0001394184127221e0,   &
         .0000273815815785e0,   &
        -.0000062750421119e0,   &
         .0000016064844184e0,   &
        -.0000004476392158e0,   &
         .0000001334635874e0,   &
        -.0000000420735334e0,   &
         .0000000139021990e0,   &
        -.0000000047831848e0,   &
         .0000000017047897e0,   &
        -.0000000006268389e0,   &
         .0000000002369824e0,   &
        -.0000000000918641e0,   &
         .0000000000364278e0,   &
        -.0000000000147475e0,   &
         .0000000000060851e0,   &
        -.0000000000025552e0,   &
         .0000000000010906e0,   &
        -.0000000000004725e0,   &
         .0000000000002076e0,   &
        -.0000000000000924e0,   &
         .0000000000000417e0,   &
        -.0000000000000190e0,   &
         .0000000000000087e0,   &
        -.0000000000000040e0,   &
         .0000000000000019e0,   &
        -.0000000000000009e0,   &
         .0000000000000004e0,   &
        -.0000000000000002e0,   &
         .0000000000000001e0,   &
        -.0000000000000000e0    ]

    integer, save :: naif   = 0
    integer, save :: naig   = 0
    integer, save :: naip   = 0
    real, save    :: x3sml  = 0.0
    real, save    :: x32sml = 0.0
    real, save    :: xbig   = 0.0

    real          :: eta, sqrtx, theta, xm, z

    !
    ! Initialisation
    !
    if ( naif == 0 ) then
        eta = 0.1*r1mach(3)
        naif   = inits (aifcs , 9, eta)
        naig   = inits (aigcs , 8, eta)
        naip   = inits (aipcs , 34, eta)

        x3sml  = eta**0.3333
        x32sml = 1.3104*x3sml**2
        xbig   = r1mach(2)**0.6666
    endif

    !
    ! Argument ranges:
    !     x < -1.0
    !     -1.0 <= x <= 1.0
    !     x > 1.0
    !

    if ( x < -1.0 ) then
        call r9aimp (x, xm, theta)
        aie = xm * cos(theta)

    elseif ( x <= 1.0 ) then
        z = 0.0
        if ( abs(x) > x3sml ) z = x**3
        aie = 0.375 + (csevl (z, aifcs, naif) - x*(0.25 + csevl (z, aigcs, naig)) )
        if ( x > x32sml ) aie = aie * exp(2.0*x*sqrt(x)/3.0)

    else
        sqrtx = sqrt(x)
        z = -1.0
        if (x < xbig) z = 2.0/(x*sqrtx) - 1.0
        aie = (.28125 + csevl (z, aipcs, naip))/sqrt(sqrtx)
    endif
end function aie

! aide --
!     Original:
!     july 1980 edition.  w. fullerton, bell labs.
!
!     evaluate the derivative of the airy function for x .le. 0
!     and evaluate aid(x)*exp(zeta) for x .ge. 0 where
!     zeta = 2/3 * x**(3/2)
!
real function aide (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for aif on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   5.22e-18
!                                         log weighted error  17.28
!                               significant figures required  16.01
!                                    decimal places required  17.73
!
    real, save :: aifcs(8) = [      &
          0.10527461226531408809e0, &
          0.01183613628152997844e0, &
          0.00012328104173225664e0, &
          0.00000062261225638140e0, &
          0.00000000185298887844e0, &
          0.00000000000363328873e0, &
          0.00000000000000504622e0, &
          0.00000000000000000522e0  ]
!
! series for aig on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   3.14e-19
!                                         log weighted error  18.50
!                               significant figures required  17.44
!                                    decimal places required  18.98
!
    real, save :: aigcs(9) = [       &
          0.021233878150918666852e0, &
          0.086315930335214406752e0, &
          0.001797594720383231358e0, &
          0.000014265499875550693e0, &
          0.000000059437995283683e0, &
          0.000000000152403366479e0, &
          0.000000000000264587660e0, &
          0.000000000000000331562e0, &
          0.000000000000000000314e0  ]
!
! series for aip2 on the interval  0.00000e+00 to  1.25000e-01
!                                        with weighted error   2.15e-17
!                                         log weighted error  16.67
!                               significant figures required  14.27
!                                    decimal places required  17.26
!
    real, save :: aip2cs(15) = [   &
          0.0065457691989713757e0, &
          0.0023833724120774592e0, &
         -0.0000430700770220586e0, &
          0.0000015629125858629e0, &
         -0.0000000815417186163e0, &
          0.0000000054103738057e0, &
         -0.0000000004284130883e0, &
          0.0000000000389497963e0, &
         -0.0000000000039623161e0, &
          0.0000000000004428184e0, &
         -0.0000000000000536297e0, &
          0.0000000000000069650e0, &
         -0.0000000000000009620e0, &
          0.0000000000000001403e0, &
         -0.0000000000000000215e0  ]
!
! series for aip1 on the interval  1.25000e-01 to  1.00000e+00
!                                        with weighted error   2.60e-17
!                                         log weighted error  16.58
!                               significant figures required  14.91
!                                    decimal places required  17.28
!
    real, save :: aip1cs(25) = [   &
          0.0358865097808301538e0, &
          0.0114668575627764899e0, &
         -0.0007592073583861400e0, &
          0.0000869517610893841e0, &
         -0.0000128237294298592e0, &
          0.0000022062695681038e0, &
         -0.0000004222295185921e0, &
          0.0000000874686415726e0, &
         -0.0000000192773588418e0, &
          0.0000000044668460054e0, &
         -0.0000000010790108052e0, &
          0.0000000002700029447e0, &
         -0.0000000000696480108e0, &
          0.0000000000184489907e0, &
         -0.0000000000050027817e0, &
          0.0000000000013852243e0, &
         -0.0000000000003908218e0, &
          0.0000000000001121536e0, &
         -0.0000000000000326862e0, &
          0.0000000000000096619e0, &
         -0.0000000000000028935e0, &
          0.0000000000000008770e0, &
         -0.0000000000000002688e0, &
          0.0000000000000000832e0, &
         -0.0000000000000000260e0  ]

    integer, save :: naif   = 0
    integer, save :: naig   = 0
    integer, save :: naip1  = 0
    integer, save :: naip2  = 0
    real, save    :: x2sml  = 0.0
    real, save    :: x3sml  = 0.0
    real, save    :: x32sml = 0.0
    real, save    :: xbig   = 0.0

    real          :: eta, xn, x2, x3, sqrtx, phi, z

    !
    ! Initialisation
    !
    if ( naif == 0 ) then
        eta    = 0.1*r1mach(3)
        naif   = inits (aifcs, 8, eta)
        naig   = inits (aigcs, 9, eta)
        naip1  = inits (aip1cs, 25, eta)
        naip2  = inits (aip2cs, 15, eta)

        x2sml  = sqrt (eta)
        x3sml  = eta**0.3333
        x32sml = 1.3104*x3sml**2
        xbig  = r1mach(2)**0.6666
    endif

    !
    ! Argument ranges:
    !     x < -1.0
    !     -1.0 <= x <= 1.0
    !      1.0 <= x <= 4.0
    !     x > 4.0
    !
    if ( x < -1.0 ) then
        call r9admp (x, xn, phi)
        aide = xn * cos(phi)

    elseif( x <= 1.0 ) then
        x3 = 0.0
        if ( abs(x) > x3sml ) x3 = x**3
        x2 = 0.0
        if ( abs(x) > x2sml ) x2 = x*x
        aide = (x2*(0.125 + csevl (x3, aifcs, naif)) - csevl (x3, aigcs, naig)) - 0.25

        if ( x > x32sml ) aide = aide * exp (2.0*x*sqrt(x)/3.0)

    elseif ( x <= 4.0 ) then
        sqrtx = sqrt(x)
        z = (16.0/(x*sqrtx) - 9.0)/7.0
        aide = (-0.28125 - csevl (z, aip1cs, naip1)) * sqrt(sqrtx)

    else
        sqrtx = sqrt(x)
        z = -1.0
        if ( x < xbig ) z = 16.0/(x*sqrtx) - 1.0
        aide = (-0.28125 - csevl (z, aip2cs, naip2)) * sqrt(sqrtx)
    endif
end function aide

! bi --
!     Original:
! july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
real function bi (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bif        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.88e-19
!                                         log weighted error  18.72
!                               significant figures required  17.74
!                                    decimal places required  19.20
!
    real, save :: bifcs(9) = [ &
        -.01673021647198664948e0, &
         .1025233583424944561e0,  &
         .00170830925073815165e0, &
         .00001186254546774468e0, &
         .00000004493290701779e0, &
         .00000000010698207143e0, &
         .00000000000017480643e0, &
         .00000000000000020810e0, &
         .00000000000000000018e0  ]
!
! series for big        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   2.61e-17
!                                         log weighted error  16.58
!                               significant figures required  15.17
!                                    decimal places required  17.03
!
    real, save :: bigcs(8) = [ &
        .02246622324857452e0, &
        .03736477545301955e0, &
        .00044476218957212e0, &
        .00000247080756363e0, &
        .00000000791913533e0, &
        .00000000001649807e0, &
        .00000000000002411e0, &
        .00000000000000002e0  ]
!
! series for bif2       on the interval  1.00000d+00 to  8.00000d+00
!                                        with weighted error   1.11e-17
!                                         log weighted error  16.95
!                        approx significant figures required  16.5
!                                    decimal places required  17.45
!
    real, save :: bif2cs(10) = [ &
        0.09984572693816041e0,   &
         .478624977863005538e0,  &
         .0251552119604330118e0, &
         .0005820693885232645e0, &
         .0000074997659644377e0, &
         .0000000613460287034e0, &
         .0000000003462753885e0, &
         .0000000000014288910e0, &
         .0000000000000044962e0, &
         .0000000000000000111e0  ]
!
! series for big2       on the interval  1.00000d+00 to  8.00000d+00
!                                        with weighted error   1.19e-18
!                                         log weighted error  17.92
!                        approx significant figures required  17.2
!                                    decimal places required  18.42
!
    real, save :: big2cs(10) = [ &
        .033305662145514340e0,   &
        .161309215123197068e0,   &
        .0063190073096134286e0,  &
        .0001187904568162517e0,  &
        .0000013045345886200e0,  &
        .0000000093741259955e0,  &
        .0000000000474580188e0,  &
        .0000000000001783107e0,  &
        .0000000000000005167e0,  &
        .0000000000000000011e0   ]

    integer, save :: nbif  = 0
    integer, save :: nbig  = 0
    integer, save :: nbif2 = 0
    integer, save :: nbig2 = 0
    real, save    :: x3sml = 0.0
    real, save    :: xmax  = 0.0

    real          :: eta, theta, xm, z

    !
    ! Initialisation
    !
    if ( nbif == 0 ) then
        eta   = 0.1*r1mach(3)
        nbif  = inits (bifcs , 9, eta)
        nbig  = inits (bigcs , 8, eta)
        nbif2 = inits (bif2cs, 10, eta)
        nbig2 = inits (big2cs, 10, eta)
        x3sml = eta**0.3333
        xmax  = (1.5*log(r1mach(2)))**0.6666
    endif

    !
    ! Argument ranges:
    !             x <  -1.0
    !     -1.0 <= x <=  1.0
    !      1.0 <  x <=  2.0
    !      2.0 <  x <=  xmax
    !             x >   xmax
    !

    if ( x < -1.0 ) then
        call r9aimp (x, xm, theta)
        bi = xm * sin(theta)

    elseif ( x <= 1.0 ) then
        z = 0.0
        if (abs(x) > x3sml) z = x**3
        bi = 0.625 + csevl (z, bifcs, nbif) + x*(0.4375 + csevl (z, bigcs, nbig))

    elseif ( x <= 2.0 ) then

         z = (2.0*x**3 - 9.0) / 7.0
         bi = 1.125 + csevl (z, bif2cs, nbif2) + x*(0.625 + csevl (z, big2cs, nbig2))

    elseif ( x <= xmax ) then
         bi = bie(x) * exp(2.0*x*sqrt(x)/3.0)

    else
         bi = ieee_value( x, ieee_positive_inf )
    endif
end function bi

! bid.f90 --
!     Original:
!     july 1980 edition.  w. fullerton, bell labs.
!
!     evaluate the derivative of the airy function bi(x).
!

real function bid (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bif on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   9.05e-18
!                                         log weighted error  17.04
!                               significant figures required  15.83
!                                    decimal places required  17.49
!
    real, save :: bifcs(8) = [     &
          0.1153536790828570243e0, &
          0.0205007894049192875e0, &
          0.0002135290278902876e0, &
          0.0000010783960614677e0, &
          0.0000000032094708833e0, &
          0.0000000000062930407e0, &
          0.0000000000000087403e0, &
          0.0000000000000000090e0  ]
!
! series for big on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   5.44e-19
!                                         log weighted error  18.26
!                               significant figures required  17.46
!                                    decimal places required  18.74
!
    real, save :: bigcs(9) = [       &
         -0.097196440416443537390e0, &
          0.149503576843167066571e0, &
          0.003113525387121326042e0, &
          0.000024708570579821297e0, &
          0.000000102949627731379e0, &
          0.000000000263970373987e0, &
          0.000000000000458279271e0, &
          0.000000000000000574283e0, &
          0.000000000000000000544e0  ]
!
! series for bif2 on the interval  1.00000e+00 to  8.00000e+00
!                                        with weighted error   3.82e-19
!                                         log weighted error  18.42
!                               significant figures required  17.68
!                                    decimal places required  18.92
!
    real, save :: bif2cs(10) = [     &
          0.323493987603522033521e0, &
          0.086297871535563559139e0, &
          0.002994025552655397426e0, &
          0.000051430528364661637e0, &
          0.000000525840250036811e0, &
          0.000000003561751373958e0, &
          0.000000000017146864007e0, &
          0.000000000000061663520e0, &
          0.000000000000000171911e0, &
          0.000000000000000000382e0  ]
!
! series for big2 on the interval  1.00000e+00 to  8.00000e+00
!                                        with weighted error   3.35e-17
!                                         log weighted error  16.48
!                               significant figures required  16.52
!                                    decimal places required  16.98
!
    real, save :: big2cs(10) = [     &
          1.6062999463621294578e0, &
          0.7449088819876088652e0, &
          0.0470138738610277380e0, &
          0.0012284422062548239e0, &
          0.0000173222412256624e0, &
          0.0000001521901652368e0, &
          0.0000000009113560249e0, &
          0.0000000000039547918e0, &
          0.0000000000000130017e0, &
          0.0000000000000000335e0  ]

    integer, save :: nbif  = 0
    integer, save :: nbig  = 0
    integer, save :: nbif2 = 0
    integer, save :: nbig2 = 0

    real, save    :: x2sml = 0.0
    real, save    :: x3sml = 0.0
    real, save    :: xmax  = 0.0

    real          :: eta, x2, x3, xn, z, phi

    !
    ! Initialisation
    !
    if ( nbif == 0 ) then
        eta   = 0.1*r1mach(3)
        nbif  = inits (bifcs, 8, eta)
        nbig  = inits (bigcs, 9, eta)
        nbif2 = inits (bif2cs, 10, eta)
        nbig2 = inits (big2cs, 10, eta)

        x2sml = sqrt (eta)
        x3sml = eta**0.3333
        xmax = (1.5*log(r1mach(2)))**0.6666
    endif

    !
    ! Argument ranges:
    !             x <  -1.0
    !     -1.0 <= x <=  1.0
    !      1.0 <  x <=  2.0
    !      2.0 <  x <=  xmax
    !             x >   xmax
    !
    if ( x < -1.0 ) then
        call r9admp (x, xn, phi)
        bid = xn * sin (phi)

    elseif ( x <= 1.0 ) then
        x3 = 0.0
        if ( abs(x) > x3sml ) x3 = x**3
        x2 = 0.0
        if (abs(x)  > x2sml ) x2 = x*x
        bid = x2*(csevl (x3, bifcs, nbif) + 0.25) + csevl (x3, bigcs, nbig) + 0.5

    elseif ( x <= 2.0 ) then
        z = (2.0*x**3 - 9.0) / 7.0
        bid = x*x*(csevl (z, bif2cs, nbif2) + 0.25) + csevl (z, big2cs, nbig2) + 0.5

    else
        if (x <= xmax ) then
            bid = bide(x) * exp (2.0*x*sqrt(x)/3.0)
        else
            bid = ieee_value( x, ieee_positive_inf )
        endif
    endif
end function bid

! bie --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate bi(x) for x .le. 0  and  bi(x)*exp(-zeta)  where
!     zeta = 2/3 * x**(3/2)  for x .ge. 0.0

real function bie (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bif        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   1.88e-19
!                                         log weighted error  18.72
!                               significant figures required  17.74
!                                    decimal places required  19.20
!
    real, save :: bifcs(9) = [    &
        -.01673021647198664948e0, &
         .1025233583424944561e0,  &
         .00170830925073815165e0, &
         .00001186254546774468e0, &
         .00000004493290701779e0, &
         .00000000010698207143e0, &
         .00000000000017480643e0, &
         .00000000000000020810e0, &
         .00000000000000000018e0  ]
!
! series for big        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   2.61e-17
!                                         log weighted error  16.58
!                               significant figures required  15.17
!                                    decimal places required  17.03
!
    real, save :: bigcs(8) = [    &
         .02246622324857452e0,    &
         .03736477545301955e0,    &
         .00044476218957212e0,    &
         .00000247080756363e0,    &
         .00000000791913533e0,    &
         .00000000001649807e0,    &
         .00000000000002411e0,    &
         .00000000000000002e0     ]
!
! series for bif2       on the interval  1.00000d+00 to  8.00000d+00
!                                        with weighted error   1.11e-17
!                                         log weighted error  16.95
!                        approx significant figures required  16.5
!                                    decimal places required  17.45
!
    real, save :: bif2cs(10) = [  &
        0.09984572693816041e0,    &
         .478624977863005538e0,   &
         .0251552119604330118e0,  &
         .0005820693885232645e0,  &
         .0000074997659644377e0,  &
         .0000000613460287034e0,  &
         .0000000003462753885e0,  &
         .0000000000014288910e0,  &
         .0000000000000044962e0,  &
         .0000000000000000111e0   ]
!
! series for big2       on the interval  1.00000d+00 to  8.00000d+00
!                                        with weighted error   1.19e-18
!                                         log weighted error  17.92
!                        approx significant figures required  17.2
!                                    decimal places required  18.42
!
    real, save :: big2cs(10) = [  &
         .033305662145514340e0,   &
         .161309215123197068e0,   &
         .0063190073096134286e0,  &
         .0001187904568162517e0,  &
         .0000013045345886200e0,  &
         .0000000093741259955e0,  &
         .0000000000474580188e0,  &
         .0000000000001783107e0,  &
         .0000000000000005167e0,  &
         .0000000000000000011e0   ]
!
! series for bip        on the interval  1.25000d-01 to  3.53553d-01
!                                        with weighted error   1.91e-17
!                                         log weighted error  16.72
!                               significant figures required  15.35
!                                    decimal places required  17.41
!
    real, save :: bipcs(24) = [ &
        -.08322047477943447e0,  &
         .01146118927371174e0,  &
         .00042896440718911e0,  &
        -.00014906639379950e0,  &
        -.00001307659726787e0,  &
         .00000632759839610e0,  &
        -.00000042226696982e0,  &
        -.00000019147186298e0,  &
         .00000006453106284e0,  &
        -.00000000784485467e0,  &
        -.00000000096077216e0,  &
         .00000000070004713e0,  &
        -.00000000017731789e0,  &
         .00000000002272089e0,  &
         .00000000000165404e0,  &
        -.00000000000185171e0,  &
         .00000000000059576e0,  &
        -.00000000000012194e0,  &
         .00000000000001334e0,  &
         .00000000000000172e0,  &
        -.00000000000000145e0,  &
         .00000000000000049e0,  &
        -.00000000000000011e0,  &
         .00000000000000001e0   ]
!
! series for bip2       on the interval  0.          to  1.25000d-01
!                                        with weighted error   1.05e-18
!                                         log weighted error  17.98
!                               significant figures required  16.74
!                                    decimal places required  18.71
!
    real, save :: bip2cs(29) = [   &
        -.113596737585988679e0,    &
         .0041381473947881595e0,   &
         .0001353470622119332e0,   &
         .0000104273166530153e0,   &
         .0000013474954767849e0,   &
         .0000001696537405438e0,   &
        -.0000000100965008656e0,   &
        -.0000000167291194937e0,   &
        -.0000000045815364485e0,   &
         .0000000003736681366e0,   &
         .0000000005766930320e0,   &
         .0000000000621812650e0,   &
        -.0000000000632941202e0,   &
        -.0000000000149150479e0,   &
         .0000000000078896213e0,   &
         .0000000000024960513e0,   &
        -.0000000000012130075e0,   &
        -.0000000000003740493e0,   &
         .0000000000002237727e0,   &
         .0000000000000474902e0,   &
        -.0000000000000452616e0,   &
        -.0000000000000030172e0,   &
         .0000000000000091058e0,   &
        -.0000000000000009814e0,   &
        -.0000000000000016429e0,   &
         .0000000000000005533e0,   &
         .0000000000000002175e0,   &
        -.0000000000000001737e0,   &
        -.0000000000000000010e0    ]

    integer, save   :: nbif   = 0
    integer, save   :: nbig   = 0
    integer, save   :: nbif2  = 0
    integer, save   :: nbig2  = 0
    integer, save   :: nbip   = 0
    integer, save   :: nbip2  = 0

    real, save      :: x3sml  = 0.0
    real, save      :: x32sml = 0.0
    real, save      :: xbig   = 0.0

    real            :: eta, xm, theta, z, sqrtx

    !
    ! Initialisation
    !
    if ( nbif == 0 ) then
        eta   = 0.1*r1mach(3)
        nbif  = inits (bifcs, 9, eta)
        nbig  = inits (bigcs, 8, eta)
        nbif2 = inits (bif2cs, 10, eta)
        nbig2 = inits (big2cs, 10, eta)
        nbip  = inits (bipcs , 24, eta)
        nbip2 = inits (bip2cs, 29, eta)

        x3sml = eta**0.3333
        x32sml = 1.3104*x3sml**2
        xbig = r1mach(2)**0.6666
    endif

    !
    ! Argument ranges:
    !             x <  -1.0
    !     -1.0 <= x <=  1.0
    !      1.0 <  x <=  2.0
    !      1.0 <  x <=  4.0
    !      4.0 <  x <=  xmax
    !             x >   xmax
    !
      if ( x < -1.0 ) then
          call r9aimp (x, xm, theta)
         bie = xm * sin(theta)

      elseif (x <= 1.0 ) then
          z = 0.0
          if ( abs(x) > x3sml ) z = x**3
          bie = 0.625 + csevl (z, bifcs, nbif) + x*(0.4375 + csevl (z, bigcs, nbig))
          if ( x > x32sml ) bie = bie * exp(-2.0*x*sqrt(x)/3.0)

      elseif ( x <= 2.0 ) then
          z = (2.0*x**3 - 9.0) / 7.0
          bie = exp(-2.0*x*sqrt(x)/3.0) * (1.125 + csevl (z, bif2cs, nbif2) + x*(0.625 + csevl (z, big2cs, nbig2)) )

      elseif ( x <= 4.0 ) then
          sqrtx = sqrt(x)
          z = atr/(x*sqrtx) + btr
          bie = (0.625 + csevl (z, bipcs, nbip)) / sqrt(sqrtx)

      else
          sqrtx = sqrt(x)
          z = -1.0
          if ( x < xbig ) z = 16.0/(x*sqrtx) - 1.0
          bie = (0.625 + csevl (z, bip2cs, nbip2))/sqrt(sqrtx)
      endif
end function bie

! bide --
!     Original:
! july 1980 edition.  w. fullerton, bell labs.
!
! evaluate the dervative of the airy function bi(x) for x.lt.0
! and evaluate bid(x) * exp(-zeta) for x .ge. 0 where
! zeta = 2/3 * x**(3/2)

real function bide (x)
    use ieee_arithmetic

    implicit none

    real, intent(in) :: x

!
! series for bif on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   9.05e-18
!                                         log weighted error  17.04
!                               significant figures required  15.83
!                                    decimal places required  17.49
!
    real, save :: bifcs(8) = [     &
          0.1153536790828570243e0, &
          0.0205007894049192875e0, &
          0.0002135290278902876e0, &
          0.0000010783960614677e0, &
          0.0000000032094708833e0, &
          0.0000000000062930407e0, &
          0.0000000000000087403e0, &
          0.0000000000000000090e0  ]
!
! series for big on the interval -1.00000e+00 to  1.00000e+00
!                                        with weighted error   5.44e-19
!                                         log weighted error  18.26
!                               significant figures required  17.46
!                                    decimal places required  18.74
!
    real, save :: bigcs(9) = [       &
         -0.097196440416443537390e0, &
          0.149503576843167066571e0, &
          0.003113525387121326042e0, &
          0.000024708570579821297e0, &
          0.000000102949627731379e0, &
          0.000000000263970373987e0, &
          0.000000000000458279271e0, &
          0.000000000000000574283e0, &
          0.000000000000000000544e0  ]
!
! series for bif2 on the interval  1.00000e+00 to  8.00000e+00
!                                        with weighted error   3.82e-19
!                                         log weighted error  18.42
!                               significant figures required  17.68
!                                    decimal places required  18.92
!
    real, save :: bif2cs(10) = [     &
          0.323493987603522033521e0, &
          0.086297871535563559139e0, &
          0.002994025552655397426e0, &
          0.000051430528364661637e0, &
          0.000000525840250036811e0, &
          0.000000003561751373958e0, &
          0.000000000017146864007e0, &
          0.000000000000061663520e0, &
          0.000000000000000171911e0, &
          0.000000000000000000382e0  ]
!
! series for big2 on the interval  1.00000e+00 to  8.00000e+00
!                                        with weighted error   3.35e-17
!                                         log weighted error  16.48
!                               significant figures required  16.52
!                                    decimal places required  16.98
!
    real, save :: big2cs(10) = [    &
          1.6062999463621294578e0,  &
          0.7449088819876088652e0,  &
          0.0470138738610277380e0,  &
          0.0012284422062548239e0,  &
          0.0000173222412256624e0,  &
          0.0000001521901652368e0,  &
          0.0000000009113560249e0,  &
          0.0000000000039547918e0,  &
          0.0000000000000130017e0,  &
          0.0000000000000000335e0   ]
!
! series for bip2 on the interval  0.00000e+00 to  1.25000e-01
!                                        with weighted error   2.07e-18
!                                         log weighted error  17.69
!                               significant figures required  16.51
!                                    decimal places required  18.42
!
    real, save :: bip2cs(29) = [    &
         -0.13269705443526630495e0, &
         -0.00568443626045977481e0, &
         -0.00015643601119611610e0, &
         -0.00001136737203679562e0, &
         -0.00000143464350991284e0, &
         -0.00000018098531185164e0, &
          0.00000000926177343611e0, &
          0.00000001710005490721e0, &
          0.00000000476698163504e0, &
         -0.00000000035195022023e0, &
         -0.00000000058890614316e0, &
         -0.00000000006678499608e0, &
          0.00000000006395565102e0, &
          0.00000000001554529427e0, &
         -0.00000000000792397000e0, &
         -0.00000000000258326243e0, &
          0.00000000000121655048e0, &
          0.00000000000038707207e0, &
         -0.00000000000022487045e0, &
         -0.00000000000004953477e0, &
          0.00000000000004563782e0, &
          0.00000000000000332998e0, &
         -0.00000000000000921750e0, &
          0.00000000000000094157e0, &
          0.00000000000000167154e0, &
         -0.00000000000000055134e0, &
         -0.00000000000000022369e0, &
          0.00000000000000017487e0, &
          0.00000000000000000207e0  ]
!
! series for bip1 on the interval  1.25000e-01 to  3.53553e-01
!                                        with weighted error   1.86e-17
!                                         log weighted error  16.73
!                               significant figures required  15.67
!                                    decimal places required  17.42
!
    real, save :: bip1cs(24) = [    &
         -0.1729187351079553719e0,  &
         -0.0149358492984694364e0,  &
         -0.0005471104951678566e0,  &
          0.0001537966292958408e0,  &
          0.0000154353476192179e0,  &
         -0.0000065434113851906e0,  &
          0.0000003728082407879e0,  &
          0.0000002072078388189e0,  &
         -0.0000000658173336470e0,  &
          0.0000000074926746354e0,  &
          0.0000000011101336884e0,  &
         -0.0000000007265140553e0,  &
          0.0000000001782723560e0,  &
         -0.0000000000217346352e0,  &
         -0.0000000000020302035e0,  &
          0.0000000000019311827e0,  &
         -0.0000000000006044953e0,  &
          0.0000000000001209450e0,  &
         -0.0000000000000125109e0,  &
         -0.0000000000000019917e0,  &
          0.0000000000000015154e0,  &
         -0.0000000000000004977e0,  &
          0.0000000000000001155e0,  &
         -0.0000000000000000186e0   ]

    integer, save   :: nbif   = 0
    integer, save   :: nbig   = 0
    integer, save   :: nbif2  = 0
    integer, save   :: nbig2  = 0
    integer, save   :: nbip1  = 0
    integer, save   :: nbip2  = 0
    real, save      :: x2sml  = 0.0
    real, save      :: x3sml  = 0.0
    real, save      :: x32sml = 0.0
    real, save      :: xbig   = 0.0

    real            :: eta, phi, sqrtx, x2, x3, xn, z

    !
    ! Initialisation
    !
    if ( nbif == 0 ) then
        eta = 0.1*r1mach(3)
        nbif   = inits (bifcs, 8, eta)
        nbig   = inits (bigcs, 9, eta)
        nbif2  = inits (bif2cs, 10, eta)
        nbig2  = inits (big2cs, 10, eta)
        nbip1  = inits (bip1cs, 24, eta)
        nbip2  = inits (bip2cs, 29, eta)

        x2sml  = sqrt (eta)
        x3sml  = eta**0.3333
        x32sml = 1.3104*x3sml**2
        xbig   = r1mach(2)**0.6666
    endif

    !
    ! Argument ranges:
    !     x < -1.0
    !     -1.0 <= x <= 1.0
    !      1.0 <= x <= 2.0
    !      2.0 <= x <= 4.0
    !     x > 4.0
    !
    if ( x < -1.0 ) then
        call r9admp (x, xn, phi)
        bide = xn * sin (phi)

    elseif ( x <= 1.0 ) then
        x3 = 0.0
        if ( abs(x) > x3sml ) x3 = x**3
        x2 = 0.0
        if ( abs(x) > x2sml ) x2 = x*x
        bide = x2 * (csevl(x3, bifcs, nbif) + 0.25) + csevl (x3, bigcs, nbig) + 0.5
        if (x > x32sml) bide = bide * exp (-2.0*x*sqrt(x)/3.0)

    elseif ( x <= 2.0 ) then
        z = (2.0*x**3 - 9.0) / 7.0
        bide = exp (-2.0*x*sqrt(x)/3.0) * (x*x * (0.25 + csevl (z, bif2cs, nbif2)) + 0.5 + csevl (z, big2cs, nbig2))

    elseif ( x <= 4.0 ) then
        sqrtx = sqrt(x)
        z = atr/(x*sqrtx) + btr
        bide = (0.625 + csevl (z, bip1cs, nbip1)) * sqrt(sqrtx)

    else
       sqrtx = sqrt(x)
        z = -1.0
        if (x < xbig ) z = 16.0/(x*sqrtx) - 1.0
        bide = (0.625 + csevl (z, bip2cs, nbip2)) * sqrt(sqrtx)
    endif
end function bide

! r9admp --
!     Original:
!     july 1980 edition.  w. fullerton, bell labs.  revised nov 1983.
!
!     evaluate the derivative of the airy function modulus and
!     phase for x .le. -1.0.
!
subroutine r9admp (x, ampl, phi)
    use ieee_arithmetic

    implicit none

    real, intent(in)  :: x
    real, intent(out) :: ampl, phi

!
! series for an22 on the interval -1.00000e+00 to -1.25000e-01
!                                        with weighted error   3.30e-17
!                                         log weighted error  16.48
!                               significant figures required  14.95
!                                    decimal places required  17.24
!
    real, save :: an22cs(33) = [  &
         0.0537418629629794329e0, &
        -0.0126661435859883193e0, &
        -0.0011924334106593007e0, &
        -0.0002032327627275655e0, &
        -0.0000446468963075164e0, &
        -0.0000113359036053123e0, &
        -0.0000031641352378546e0, &
        -0.0000009446708886149e0, &
        -0.0000002966562236472e0, &
        -0.0000000969118892024e0, &
        -0.0000000326822538653e0, &
        -0.0000000113144618964e0, &
        -0.0000000040042691002e0, &
        -0.0000000014440333684e0, &
        -0.0000000005292853746e0, &
        -0.0000000001967763374e0, &
        -0.0000000000740800096e0, &
        -0.0000000000282016314e0, &
        -0.0000000000108440066e0, &
        -0.0000000000042074801e0, &
        -0.0000000000016459150e0, &
        -0.0000000000006486827e0, &
        -0.0000000000002574095e0, &
        -0.0000000000001027889e0, &
        -0.0000000000000412846e0, &
        -0.0000000000000166711e0, &
        -0.0000000000000067657e0, &
        -0.0000000000000027585e0, &
        -0.0000000000000011296e0, &
        -0.0000000000000004645e0, &
        -0.0000000000000001917e0, &
        -0.0000000000000000794e0, &
        -0.0000000000000000330e0  ]
!
! series for an21 on the interval -1.25000e-01 to -1.56250e-02
!                                        with weighted error   3.43e-17
!                                         log weighted error  16.47
!                               significant figures required  14.48
!                                    decimal places required  17.16
!
    real, save :: an21cs(24) = [  &
         0.0198313155263169394e0, &
        -0.0029376249067087533e0, &
        -0.0001136260695958196e0, &
        -0.0000100554451087156e0, &
        -0.0000013048787116563e0, &
        -0.0000002123881993151e0, &
        -0.0000000402270833384e0, &
        -0.0000000084996745953e0, &
        -0.0000000019514839426e0, &
        -0.0000000004783865344e0, &
        -0.0000000001236733992e0, &
        -0.0000000000334137486e0, &
        -0.0000000000093702824e0, &
        -0.0000000000027130128e0, &
        -0.0000000000008075954e0, &
        -0.0000000000002463214e0, &
        -0.0000000000000767656e0, &
        -0.0000000000000243883e0, &
        -0.0000000000000078831e0, &
        -0.0000000000000025882e0, &
        -0.0000000000000008619e0, &
        -0.0000000000000002908e0, &
        -0.0000000000000000993e0, &
        -0.0000000000000000343e0  ]
!
! series for an20 on the interval -1.56250e-02 to  0.00000e+00
!                                        with weighted error   4.41e-17
!                                         log weighted error  16.36
!                               significant figures required  14.16
!                                    decimal places required  16.96
!
    real, save :: an20cs(16) = [  &
         0.0126732217145738027e0, &
        -0.0005212847072615621e0, &
        -0.0000052672111140370e0, &
        -0.0000001628202185026e0, &
        -0.0000000090991442687e0, &
        -0.0000000007438647126e0, &
        -0.0000000000795494752e0, &
        -0.0000000000104050944e0, &
        -0.0000000000015932426e0, &
        -0.0000000000002770648e0, &
        -0.0000000000000535343e0, &
        -0.0000000000000113062e0, &
        -0.0000000000000025772e0, &
        -0.0000000000000006278e0, &
        -0.0000000000000001621e0, &
        -0.0000000000000000441e0  ]
!
! series for aph2 on the interval -1.00000e+00 to -1.25000e-01
!                                        with weighted error   2.94e-17
!                                         log weighted error  16.53
!                               significant figures required  15.58
!                                    decimal places required  17.28
!
    real, save :: aph2cs(32) = [ &
       -0.2057088719781465107e0, &
        0.0422196961357771922e0, &
        0.0020482560511207275e0, &
        0.0002607800735165006e0, &
        0.0000474824268004729e0, &
        0.0000105102756431612e0, &
        0.0000026353534014668e0, &
        0.0000007208824863499e0, &
        0.0000002103236664473e0, &
        0.0000000644975634555e0, &
        0.0000000205802377264e0, &
        0.0000000067836273921e0, &
        0.0000000022974015284e0, &
        0.0000000007961306765e0, &
        0.0000000002813860610e0, &
        0.0000000001011749057e0, &
        0.0000000000369306738e0, &
        0.0000000000136615066e0, &
        0.0000000000051142751e0, &
        0.0000000000019351689e0, &
        0.0000000000007393607e0, &
        0.0000000000002849792e0, &
        0.0000000000001107281e0, &
        0.0000000000000433412e0, &
        0.0000000000000170801e0, &
        0.0000000000000067733e0, &
        0.0000000000000027017e0, &
        0.0000000000000010835e0, &
        0.0000000000000004367e0, &
        0.0000000000000001769e0, &
        0.0000000000000000719e0, &
        0.0000000000000000294e0  ]
!
! series for aph1 on the interval -1.25000e-01 to -1.56250e-02
!                                        with weighted error   6.38e-17
!                                         log weighted error  16.20
!                               significant figures required  14.91
!                                    decimal places required  16.87
!
    real, save :: aph1cs(22) = [ &
       -0.1024172908077571694e0, &
        0.0071697275146591248e0, &
        0.0001209959363122329e0, &
        0.0000073361512841220e0, &
        0.0000007535382954272e0, &
        0.0000001041478171741e0, &
        0.0000000174358728519e0, &
        0.0000000033399795033e0, &
        0.0000000007073075174e0, &
        0.0000000001619187515e0, &
        0.0000000000394539982e0, &
        0.0000000000101192282e0, &
        0.0000000000027092778e0, &
        0.0000000000007523806e0, &
        0.0000000000002156369e0, &
        0.0000000000000635283e0, &
        0.0000000000000191757e0, &
        0.0000000000000059143e0, &
        0.0000000000000018597e0, &
        0.0000000000000005950e0, &
        0.0000000000000001934e0, &
        0.0000000000000000638e0  ]
!
! series for aph0 on the interval -1.56250e-02 to  0.00000e+00
!                                        with weighted error   2.29e-17
!                                         log weighted error  16.64
!                               significant figures required  15.27
!                                    decimal places required  17.23
!
    real, save :: aph0cs(15) = [ &
       -0.0855849241130933257e0, &
        0.0011214378867065261e0, &
        0.0000042721029353664e0, &
        0.0000000817607381483e0, &
        0.0000000033907645000e0, &
        0.0000000002253264423e0, &
        0.0000000000206284209e0, &
        0.0000000000023858763e0, &
        0.0000000000003301618e0, &
        0.0000000000000527010e0, &
        0.0000000000000094555e0, &
        0.0000000000000018709e0, &
        0.0000000000000004024e0, &
        0.0000000000000000930e0, &
        0.0000000000000000229e0  ]

    real, save    :: pi34 = 2.3561944901923449e0

!   real, save :: pi34 = 2.35619449019234492884698253745962716313d0, &

    integer, save :: nan20  = 0
    integer, save :: nan21  = 0
    integer, save :: nan22  = 0
    integer, save :: naph0  = 0
    integer, save :: naph1  = 0
    integer, save :: naph2  = 0

    real, save    :: xsml   = 0.0

    real          :: eta, z, sqrtx

    !
    ! Initialisation
    !
    if ( nan20 == 0 ) then
        eta   = 0.1*r1mach(3)
        nan20 = inits (an20cs, 16, eta)
        nan21 = inits (an21cs, 24, eta)
        nan22 = inits (an22cs, 33, eta)
        naph0 = inits (aph0cs, 15, eta)
        naph1 = inits (aph1cs, 22, eta)
        naph2 = inits (aph2cs, 32, eta)

        xsml  = -(128.0/r1mach(3))**0.3333
    endif

    !
    ! Argument ranges:
    !             x <  -4.0
    !     -4.0 <= x <  -2.0
    !     -2.0 <= x <= -1.0
    !             x >  -1.0
    !
    if ( x < -4.0 ) then
        z = 1.0
        if (x > xsml) z = 128.0/x**3 + 1.0
        ampl = 0.3125 + csevl (z, an20cs, nan20)
        phi  = -0.625 + csevl (z, aph0cs, naph0)

    elseif ( x < -2.0 ) then
        z    = (128.0/x**3 + 9.0) / 7.0
        ampl = 0.3125 + csevl (z, an21cs, nan21)
        phi  = -0.625 + csevl (z, aph1cs, naph1)

    elseif ( x <= -1.0 ) then
        z    = (16.0/x**3 + 9.0) / 7.0
        ampl = 0.3125 + csevl (z, an22cs, nan22)
        phi  = -0.625 + csevl (z, aph2cs, naph2)

    else
        ampl = ieee_value( x, ieee_quiet_nan )
        phi  = ieee_value( x, ieee_quiet_nan )
    endif

    !
    ! Finish the calculations
    !
    sqrtx = sqrt (-x)
    ampl  = sqrt (ampl*sqrtx)
    phi  = pi34 - x*sqrtx*phi
end subroutine r9admp

! r9aimp --
!     Original:
!     july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
!     evaluate the airy modulus and phase for x .le. -1.0
!
subroutine r9aimp (x, ampl, theta)
    use ieee_arithmetic

    implicit none

    real, intent(in)  :: x
    real, intent(out) :: ampl, theta

!
! series for am21       on the interval -1.25000d-01 to  0.
!                                        with weighted error   2.89e-17
!                                         log weighted error  16.54
!                               significant figures required  14.15
!                                    decimal places required  17.34
!
    real, save :: am21cs(40) = [ &
        .0065809191761485e0,     &
        .0023675984685722e0,     &
        .0001324741670371e0,     &
        .0000157600904043e0,     &
        .0000027529702663e0,     &
        .0000006102679017e0,     &
        .0000001595088468e0,     &
        .0000000471033947e0,     &
        .0000000152933871e0,     &
        .0000000053590722e0,     &
        .0000000020000910e0,     &
        .0000000007872292e0,     &
        .0000000003243103e0,     &
        .0000000001390106e0,     &
        .0000000000617011e0,     &
        .0000000000282491e0,     &
        .0000000000132979e0,     &
        .0000000000064188e0,     &
        .0000000000031697e0,     &
        .0000000000015981e0,     &
        .0000000000008213e0,     &
        .0000000000004296e0,     &
        .0000000000002284e0,     &
        .0000000000001232e0,     &
        .0000000000000675e0,     &
        .0000000000000374e0,     &
        .0000000000000210e0,     &
        .0000000000000119e0,     &
        .0000000000000068e0,     &
        .0000000000000039e0,     &
        .0000000000000023e0,     &
        .0000000000000013e0,     &
        .0000000000000008e0,     &
        .0000000000000005e0,     &
        .0000000000000003e0,     &
        .0000000000000001e0,     &
        .0000000000000001e0,     &
        .0000000000000000e0,     &
        .0000000000000000e0,     &
        .0000000000000000e0      ]
!
! series for ath1       on the interval -1.25000d-01 to  0.
!                                        with weighted error   2.53e-17
!                                         log weighted error  16.60
!                               significant figures required  15.15
!                                    decimal places required  17.38
!
    real, save :: ath1cs(36) = [ &
        -.07125837815669365e0,   &
        -.00590471979831451e0,   &
        -.00012114544069499e0,   &
        -.00000988608542270e0,   &
        -.00000138084097352e0,   &
        -.00000026142640172e0,   &
        -.00000006050432589e0,   &
        -.00000001618436223e0,   &
        -.00000000483464911e0,   &
        -.00000000157655272e0,   &
        -.00000000055231518e0,   &
        -.00000000020545441e0,   &
        -.00000000008043412e0,   &
        -.00000000003291252e0,   &
        -.00000000001399875e0,   &
        -.00000000000616151e0,   &
        -.00000000000279614e0,   &
        -.00000000000130428e0,   &
        -.00000000000062373e0,   &
        -.00000000000030512e0,   &
        -.00000000000015239e0,   &
        -.00000000000007758e0,   &
        -.00000000000004020e0,   &
        -.00000000000002117e0,   &
        -.00000000000001132e0,   &
        -.00000000000000614e0,   &
        -.00000000000000337e0,   &
        -.00000000000000188e0,   &
        -.00000000000000105e0,   &
        -.00000000000000060e0,   &
        -.00000000000000034e0,   &
        -.00000000000000020e0,   &
        -.00000000000000011e0,   &
        -.00000000000000007e0,   &
        -.00000000000000004e0,   &
        -.00000000000000002e0    ]
!
! series for am22       on the interval -1.00000d+00 to -1.25000d-01
!                                        with weighted error   2.99e-17
!                                         log weighted error  16.52
!                               significant figures required  14.57
!                                    decimal places required  17.28
!
    real, save :: am22cs(33) = [ &
        -.01562844480625341e0,   &
        .00778336445239681e0,    &
        .00086705777047718e0,    &
        .00015696627315611e0,    &
        .00003563962571432e0,    &
        .00000924598335425e0,    &
        .00000262110161850e0,    &
        .00000079188221651e0,    &
        .00000025104152792e0,    &
        .00000008265223206e0,    &
        .00000002805711662e0,    &
        .00000000976821090e0,    &
        .00000000347407923e0,    &
        .00000000125828132e0,    &
        .00000000046298826e0,    &
        .00000000017272825e0,    &
        .00000000006523192e0,    &
        .00000000002490471e0,    &
        .00000000000960156e0,    &
        .00000000000373448e0,    &
        .00000000000146417e0,    &
        .00000000000057826e0,    &
        .00000000000022991e0,    &
        .00000000000009197e0,    &
        .00000000000003700e0,    &
        .00000000000001496e0,    &
        .00000000000000608e0,    &
        .00000000000000248e0,    &
        .00000000000000101e0,    &
        .00000000000000041e0,    &
        .00000000000000017e0,    &
        .00000000000000007e0,    &
        .00000000000000002e0     ]
!
! series for ath2       on the interval -1.00000d+00 to -1.25000d-01
!                                        with weighted error   2.57e-17
!                                         log weighted error  16.59
!                               significant figures required  15.07
!                                    decimal places required  17.34
!
    real, save :: ath2cs(32) = [ &
        .00440527345871877e0,    &
        -.03042919452318455e0,   &
        -.00138565328377179e0,   &
        -.00018044439089549e0,   &
        -.00003380847108327e0,   &
        -.00000767818353522e0,   &
        -.00000196783944371e0,   &
        -.00000054837271158e0,   &
        -.00000016254615505e0,   &
        -.00000005053049981e0,   &
        -.00000001631580701e0,   &
        -.00000000543420411e0,   &
        -.00000000185739855e0,   &
        -.00000000064895120e0,   &
        -.00000000023105948e0,   &
        -.00000000008363282e0,   &
        -.00000000003071196e0,   &
        -.00000000001142367e0,   &
        -.00000000000429811e0,   &
        -.00000000000163389e0,   &
        -.00000000000062693e0,   &
        -.00000000000024260e0,   &
        -.00000000000009461e0,   &
        -.00000000000003716e0,   &
        -.00000000000001469e0,   &
        -.00000000000000584e0,   &
        -.00000000000000233e0,   &
        -.00000000000000093e0,   &
        -.00000000000000037e0,   &
        -.00000000000000015e0,   &
        -.00000000000000006e0,   &
        -.00000000000000002e0    ]

    real, save    :: pi4 = 0.78539816339744831e0

    integer, save :: nam21 = 0
    integer, save :: nath1 = 0
    integer, save :: nath2 = 0
    integer, save :: nam22 = 0

    real, save    :: xsml  = 0.0

    real          :: eta, z, sqrtx

    !
    ! Initialisation
    !
    if ( nam21 == 0 ) then
        eta   = 0.1*r1mach(3)
        nam21 = inits (am21cs, 40, eta)
        nath1 = inits (ath1cs, 36, eta)
        nam22 = inits (am22cs, 33, eta)
        nath2 = inits (ath2cs, 32, eta)

        xsml  = -(16.0/r1mach(3))**0.3333
    endif

    !
    ! Argument ranges:
    !             x <  -2.0
    !     -2.0 <= x <= -1.0
    !             x >  -1.0
    !
    if ( x < -2.0 ) then
        z = 1.0
        if (x > xsml) z = 16.0/x**3 + 1.0
        ampl  = 0.3125 + csevl(z, am21cs, nam21)
        theta = -0.625 + csevl (z, ath1cs, nath1)

    elseif ( x <= -1.0 ) then
        z     = (16.0/x**3 + 9.0)/7.0
        ampl  = 0.3125 + csevl (z, am22cs, nam22)
        theta = -0.625 + csevl (z, ath2cs, nath2)

    else
        ampl  = ieee_value( x, ieee_quiet_nan )
        theta = ieee_value( x, ieee_quiet_nan )
        return
    endif

    !
    ! Finish the calculations
    !
    sqrtx = sqrt(-x)
    ampl  = sqrt (ampl/sqrtx)
    theta = pi4 - x*sqrtx * theta
end subroutine r9aimp

end module fullerton_airy
