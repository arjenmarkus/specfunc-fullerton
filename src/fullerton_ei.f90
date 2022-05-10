! fullerton_ei.f90 --
!     Module containing functions for evaluating the exponential integral and
!     related functions
!
!     TODO:
!     Implement En(x) on the basis of recursive relations
!
module fullerton_ei
    use ieee_arithmetic
    use fullerton_aux

    implicit none

    private
    public :: integral_ei, integral_e1, integral_li, integral_si, integral_ci, integral_shi, integral_chi, integral_en

    interface integral_ei
        module procedure ei
    end interface

    interface integral_e1
        module procedure e1
    end interface

    interface integral_li
        module procedure li
    end interface

    interface integral_si
        module procedure si
    end interface

    interface integral_ci
        module procedure ci
    end interface

    interface integral_shi
        module procedure shi
    end interface

    interface integral_chi
        module procedure chi
    end interface

contains

! ei --
!     Original:
!     april 1977 version.  w. fullerton, c3, los alamos scientific lab.
!
real function ei (x)
    real, intent(in) :: x

    ei = -e1(-x)

end function ei

! e1 --
!    Original:
!    july 1977 edition.  w. fullerton, c3, los alamos scientific lab.
!
real function e1 (x)
    real, intent(in) :: x

!
! series for ae11       on the interval -1.00000d-01 to  0.
!                                        with weighted error   1.76e-17
!                                         log weighted error  16.75
!                               significant figures required  15.70
!                                    decimal places required  17.55
!
    real, save :: ae11cs(39) = [ &
         .12150323971606579e0,   &
        -.065088778513550150e0,  &
         .004897651357459670e0,  &
        -.000649237843027216e0,  &
         .000093840434587471e0,  &
         .000000420236380882e0,  &
        -.000008113374735904e0,  &
         .000002804247688663e0,  &
         .000000056487164441e0,  &
        -.000000344809174450e0,  &
         .000000058209273578e0,  &
         .000000038711426349e0,  &
        -.000000012453235014e0,  &
        -.000000005118504888e0,  &
         .000000002148771527e0,  &
         .000000000868459898e0,  &
        -.000000000343650105e0,  &
        -.000000000179796603e0,  &
         .000000000047442060e0,  &
         .000000000040423282e0,  &
        -.000000000003543928e0,  &
        -.000000000008853444e0,  &
        -.000000000000960151e0,  &
         .000000000001692921e0,  &
         .000000000000607990e0,  &
        -.000000000000224338e0,  &
        -.000000000000200327e0,  &
        -.000000000000006246e0,  &
         .000000000000045571e0,  &
         .000000000000016383e0,  &
        -.000000000000005561e0,  &
        -.000000000000006074e0,  &
        -.000000000000000862e0,  &
         .000000000000001223e0,  &
         .000000000000000716e0,  &
        -.000000000000000024e0,  &
        -.000000000000000201e0,  &
        -.000000000000000082e0,  &
         .000000000000000017e0   ]
!
! series for ae12       on the interval -2.50000d-01 to -1.00000d-01
!                                        with weighted error   5.83e-17
!                                         log weighted error  16.23
!                               significant figures required  15.76
!                                    decimal places required  16.93
!
    real, save :: ae12cs(25) = [ &
         .58241749513472674e0,   &
        -.15834885090578275e0,   &
        -.006764275590323141e0,  &
         .005125843950185725e0,  &
         .000435232492169391e0,  &
        -.000143613366305483e0,  &
        -.000041801320556301e0,  &
        -.000002713395758640e0,  &
         .000001151381913647e0,  &
         .000000420650022012e0,  &
         .000000066581901391e0,  &
         .000000000662143777e0,  &
        -.000000002844104870e0,  &
        -.000000000940724197e0,  &
        -.000000000177476602e0,  &
        -.000000000015830222e0,  &
         .000000000002905732e0,  &
         .000000000001769356e0,  &
         .000000000000492735e0,  &
         .000000000000093709e0,  &
         .000000000000010707e0,  &
        -.000000000000000537e0,  &
        -.000000000000000716e0,  &
        -.000000000000000244e0,  &
        -.000000000000000058e0   ]
!
! series for e11        on the interval -4.00000d+00 to -1.00000d+00
!                                        with weighted error   1.08e-18
!                                         log weighted error  17.97
!                               significant figures required  19.02
!                                    decimal places required  18.61
!
    real, save :: e11cs(19) = [   &
        -16.113461655571494026e0, &
         7.7940727787426802769e0, &
        -1.9554058188631419507e0, &
         .37337293866277945612e0, &
        -.05692503191092901938e0, &
         .00721107776966009185e0, &
        -.00078104901449841593e0, &
         .00007388093356262168e0, &
        -.00000620286187580820e0, &
         .00000046816002303176e0, &
        -.00000003209288853329e0, &
         .00000000201519974874e0, &
        -.00000000011673686816e0, &
         .00000000000627627066e0, &
        -.00000000000031481541e0, &
         .00000000000001479904e0, &
        -.00000000000000065457e0, &
         .00000000000000002733e0, &
        -.00000000000000000108e0  ]
!
! series for e12        on the interval -1.00000d+00 to  1.00000d+00
!                                        with weighted error   3.15e-18
!                                         log weighted error  17.50
!                        approx significant figures required  15.8
!                                    decimal places required  18.10
!
    real, save :: e12cs(16) = [   &
        -.037390214792202795e0,   &
         .042723986062209577e0,   &
        -.1303182079849700544e0,  &
         .01441912402469889073e0, &
        -.00134617078051068022e0, &
         .00010731029253063780e0, &
        -.00000742999951611943e0, &
         .00000045377325690753e0, &
        -.00000002476417211390e0, &
         .00000000122076581374e0, &
        -.00000000005485141480e0, &
         .00000000000226362142e0, &
        -.00000000000008635897e0, &
         .00000000000000306291e0, &
        -.00000000000000010148e0, &
         .00000000000000000315e0  ]
!
! series for ae13       on the interval  2.50000d-01 to  1.00000d+00
!                                        with weighted error   2.34e-17
!                                         log weighted error  16.63
!                               significant figures required  16.14
!                                    decimal places required  17.33
!
    real, save :: ae13cs(25) = [ &
        -.60577324664060346e0,   &
        -.11253524348366090e0,   &
         .013432266247902779e0,  &
        -.001926845187381145e0,  &
         .000309118337720603e0,  &
        -.000053564132129618e0,  &
         .000009827812880247e0,  &
        -.000001885368984916e0,  &
         .000000374943193568e0,  &
        -.000000076823455870e0,  &
         .000000016143270567e0,  &
        -.000000003466802211e0,  &
         .000000000758754209e0,  &
        -.000000000168864333e0,  &
         .000000000038145706e0,  &
        -.000000000008733026e0,  &
         .000000000002023672e0,  &
        -.000000000000474132e0,  &
         .000000000000112211e0,  &
        -.000000000000026804e0,  &
         .000000000000006457e0,  &
        -.000000000000001568e0,  &
         .000000000000000383e0,  &
        -.000000000000000094e0,  &
         .000000000000000023e0   ]
!
! series for ae14       on the interval  0.          to  2.50000d-01
!                                        with weighted error   5.41e-17
!                                         log weighted error  16.27
!                               significant figures required  15.38
!                                    decimal places required  16.97
!
    real, save :: ae14cs(26) = [ &
        -.1892918000753017e0,    &
        -.08648117855259871e0,   &
         .00722410154374659e0,   &
        -.00080975594575573e0,   &
         .00010999134432661e0,   &
        -.00001717332998937e0,   &
         .00000298562751447e0,   &
        -.00000056596491457e0,   &
         .00000011526808397e0,   &
        -.00000002495030440e0,   &
         .00000000569232420e0,   &
        -.00000000135995766e0,   &
         .00000000033846628e0,   &
        -.00000000008737853e0,   &
         .00000000002331588e0,   &
        -.00000000000641148e0,   &
         .00000000000181224e0,   &
        -.00000000000052538e0,   &
         .00000000000015592e0,   &
        -.00000000000004729e0,   &
         .00000000000001463e0,   &
        -.00000000000000461e0,   &
         .00000000000000148e0,   &
        -.00000000000000048e0,   &
         .00000000000000016e0,   &
        -.00000000000000005e0    ]

    integer, save :: ntae11 = 0
    integer, save :: ntae12 = 0
    integer, save :: nte11  = 0
    integer, save :: nte12  = 0
    integer, save :: ntae13 = 0
    integer, save :: ntae14 = 0
    real, save    :: xmax   = 0.0

    real          :: eta

    !
    ! Initialisation
    !
    if ( ntae11 == 0 ) then
        eta    = 0.1*r1mach(3)
        ntae11 = inits (ae11cs, 39, eta)
        ntae12 = inits (ae12cs, 25, eta)
        nte11  = inits (e11cs, 19, eta)
        nte12  = inits (e12cs, 16, eta)
        ntae13 = inits (ae13cs, 25, eta)
        ntae14 = inits (ae14cs, 26, eta)

        xmax = -alog (r1mach(1))
        xmax = xmax - alog(xmax)
    endif

    !
    ! Argument ranges:
    !     x <= -10.0
    !     -4.0 <= x <= -1.0
    !     -1.0 <= x <= 1.0
    !      1.0 <= x <= 4.0
    !      4.0 <= x <= xmax
    !     x > xmax
    !
    if ( x <= -10.0 ) then
        !
        !  e1(x) = -ei(-x) for x .le. -10.
        !
        e1 = exp(-x)/x * (1.+csevl (20./x+1., ae11cs, ntae11))

    elseif (x <= -4.0 ) then
        !
        ! e1(x) = -ei(-x) for -10. .lt. x .le. -4.
        !
        e1 = exp(-x)/x * (1.+csevl ((40./x+7.)/3., ae12cs, ntae12))

    elseif (x <= -1.0 ) then
        !
        ! e1(x) = -ei(-x) for -4. .lt. x .le. -1.
        !
        e1 = -alog(abs(x)) + csevl ((2.*x+5.)/3., e11cs, nte11)

    elseif ( x <= 1.0 ) then

        if (x == 0.0 ) then
            e1 = ieee_value( x, ieee_quiet_nan )
        else
            !
            ! e1(x) = -ei(-x) for -1. .lt. x .le. 1.,  x .ne. 0.
            !
            e1 = (-alog(abs(x)) - 0.6875 + x) + csevl (x, e12cs, nte12)
        endif

    elseif (x <= 4.0 ) then
        !
        ! e1(x) = -ei(-x) for 1. .lt. x .le. 4.
        !
        e1 = exp(-x)/x * (1.+csevl ((8./x-5.)/3., ae13cs, ntae13))

    elseif ( x <= xmax ) then
        !
        ! e1(x) = -ei(-x) for 4. .lt. x .le. xmax
        !
        e1 = exp(-x)/x * (1. + csevl (8./x-1., ae14cs, ntae14))

    else
        e1 = 0.
    endif
end function e1

! li --
!     Original:
!     august 1980 edition.   w. fullerton, c3, los alamos scientific lab.
!
real function li (x)
    real, intent(in) :: x

    real, save :: sqeps = 0.0

    !
    ! Initialisation
    !
    if ( sqeps == 0.0 ) then
        sqeps = sqrt(r1mach(3))
    endif

    !
    ! Argument ranges:
    !     x <= 0.0
    !     x >  0.0
    !
    if ( x <= 0.0 ) then
        li = ieee_value( x, ieee_quiet_nan )

    else
        li = ei (log(x) )
    endif
end function li

! si --
!     Original:
!     december 1980 edition,  w. fullerton, bell labs.
!
real function si (x)
    real, intent(in) :: x

!
! series for si   on the interval  0.00000e+00 to  1.60000e+01
!                                        with weighted error   1.22e-17
!                                         log weighted error  16.91
!                               significant figures required  16.37
!                                    decimal places required  17.45
!
    real, save :: sics(12) = [     &
         -0.1315646598184841929e0, &
         -0.2776578526973601892e0, &
          0.0354414054866659180e0, &
         -0.0025631631447933978e0, &
          0.0001162365390497009e0, &
         -0.0000035904327241606e0, &
          0.0000000802342123706e0, &
         -0.0000000013562997693e0, &
          0.0000000000179440722e0, &
         -0.0000000000001908387e0, &
          0.0000000000000016670e0, &
         -0.0000000000000000122e0  ]

    integer, save :: nsi  = 0
    real, save    :: xsml = 0.0

    real          :: absx, cosx, f, g

    !
    ! Initialisation
    !
    if ( nsi == 0 ) then
        nsi  = inits (sics, 12, 0.1*r1mach(3))
        xsml = sqrt (r1mach(3))
    endif

    !
    ! Argument ranges:
    !          abs(x) <= 4.0
    !          abs(x) >  4.0
    !
    absx = abs(x)
    if ( absx <= 4.0 ) then
        si = x
        if ( absx < xsml ) then
            return
        endif

        si = x*(0.75 + csevl ((x**2-8.0)*0.125, sics, nsi))
    else
        call r9sifg (absx, f, g)
        cosx = cos (absx)
        si   = pi2 - f*cosx - g*sin(x)
        if ( x < 0.0 ) si = -si
    endif
end function si

! ci --
!     Original:
!     december 1980 edition, w. fullerton, bell labs.
!
real function ci (x)
    real, intent(in) :: x

!
! series for ci   on the interval  0.00000e+00 to  1.60000e+01
!                                        with weighted error   1.94e-18
!                                         log weighted error  17.71
!                               significant figures required  17.74
!                                    decimal places required  18.27
!
    real, save :: cics(13) = [       &
          -0.34004281856055363156e0, &
          -1.03302166401177456807e0, &
           0.19388222659917082877e0, &
          -0.01918260436019865894e0, &
           0.00110789252584784967e0, &
          -0.00004157234558247209e0, &
           0.00000109278524300229e0, &
          -0.00000002123285954183e0, &
           0.00000000031733482164e0, &
          -0.00000000000376141548e0, &
           0.00000000000003622653e0, &
          -0.00000000000000028912e0, &
           0.00000000000000000194e0  ]

    integer, save :: nci  = 0
    real, save    :: xsml = 0.0

    real          :: f, g, y, sinx

    !
    ! Initialisation
    !
    if ( nci == 0 ) then
        nci  = inits (cics, 13, 0.1*r1mach(3))
        xsml = sqrt (r1mach(3))
    endif

    !
    ! Argument ranges:
    !     x <= 0.0
    !     0.0 <= x <= 4.0
    !     x > 4.0
    !
    if (x <= 0.0) then
       ci = ieee_value( x, ieee_quiet_nan )

    elseif ( x <= 4.0 ) then
       y = -1.0
       if (x > xsml) y = (x**2-8.0)*0.125
       ci = log(x) - 0.5 + csevl (y, cics, nci)

    else
        call r9sifg (x, f, g)
        sinx = sin (x)
        ci = f*sinx - g*cos(x)
    endif
end function ci

! shi --
!     Original:
!     december 1980 edition.  w. fullerton, bell labs.
!
!     evaluate the hyperbolic sine integral
!             shi = integral from 0 to x of  sinh(t)/t dt.
!
real function shi (x)
    real, intent(in) :: x

!
! series for shi  on the interval  0.00000e+00 to  1.40625e-01
!                                        with weighted error   4.67e-20
!                                         log weighted error  19.33
!                               significant figures required  17.07
!                                    decimal places required  19.75
!
    real, save :: shics(7) = [      &
        0.0078372685688900950695e0, &
        0.0039227664934234563973e0, &
        0.0000041346787887617267e0, &
        0.0000000024707480372883e0, &
        0.0000000000009379295591e0, &
        0.0000000000000002451817e0, &
        0.0000000000000000000467e0  ]

    integer, save :: nshi = 0
    real, save    :: xsml = 0.0

    real          :: absx

    !
    ! Initialisation
    !
    if ( nshi == 0 ) then
        nshi = inits (shics, 7, 0.1*r1mach(3))
        xsml = sqrt (r1mach(3))
    endif

    !
    ! Argument ranges:
    !     abs(x) <= xsml
    !     xsml < abs(x) <= 0.375
    !     x > 0.375
    !
    absx = abs(x)
    shi = x

    if ( absx > xsml .and. absx <= 0.375 ) then
        shi = x*(1.0 + csevl (128.*x**2/9.-1., shics, nshi))

    else
        shi = 0.5*(ei(x) + e1(x))
    endif
end function shi

! chi --
!     Original:
!     december 1980 edition.  w. fullerton, bell labs.
!
!     evaluate the hyperbolic cosine integral.  when x is negative, the
!     principal value is used.
!
real function chi (x)
    real, intent(in) :: x

    chi = 0.5 * (ei(x) - e1(x))
end function chi

! r9sifg --
!     Original:
!     december 1980 edition.  w. fullerton, bell labs
!
subroutine r9sifg (x, f, g)
    real, intent(in)  :: x
    real, intent(out) :: f, g

!
! series for f1   on the interval  2.00000e-02 to  6.25000e-02
!                                        with weighted error   2.82e-17
!                                         log weighted error  16.55
!                               significant figures required  15.36
!                                    decimal places required  17.20
!
    real, save :: f1cs(20) = [     &
         -0.1191081969051363610e0, &
         -0.0247823144996236248e0, &
          0.0011910281453357821e0, &
         -0.0000927027714388562e0, &
          0.0000093373141568271e0, &
         -0.0000011058287820557e0, &
          0.0000001464772071460e0, &
         -0.0000000210694496288e0, &
          0.0000000032293492367e0, &
         -0.0000000005206529618e0, &
          0.0000000000874878885e0, &
         -0.0000000000152176187e0, &
          0.0000000000027257192e0, &
         -0.0000000000005007053e0, &
          0.0000000000000940241e0, &
         -0.0000000000000180014e0, &
          0.0000000000000035063e0, &
         -0.0000000000000006935e0, &
          0.0000000000000001391e0, &
         -0.0000000000000000282e0  ]

!
! series for f2   on the interval  0.00000e+00 to  2.00000e-02
!                                        with weighted error   4.32e-17
!                                         log weighted error  16.36
!                               significant figures required  14.75
!                                    decimal places required  17.10
!
    real, save :: f2cs(29) = [     &
         -0.0348409253897013234e0, &
         -0.0166842205677959686e0, &
          0.0006752901241237738e0, &
         -0.0000535066622544701e0, &
          0.0000062693421779007e0, &
         -0.0000009526638801991e0, &
          0.0000001745629224251e0, &
         -0.0000000368795403065e0, &
          0.0000000087202677705e0, &
         -0.0000000022601970392e0, &
          0.0000000006324624977e0, &
         -0.0000000001888911889e0, &
          0.0000000000596774674e0, &
         -0.0000000000198044313e0, &
          0.0000000000068641396e0, &
         -0.0000000000024731020e0, &
          0.0000000000009226360e0, &
         -0.0000000000003552364e0, &
          0.0000000000001407606e0, &
         -0.0000000000000572623e0, &
          0.0000000000000238654e0, &
         -0.0000000000000101714e0, &
          0.0000000000000044259e0, &
         -0.0000000000000019634e0, &
          0.0000000000000008868e0, &
         -0.0000000000000004074e0, &
          0.0000000000000001901e0, &
         -0.0000000000000000900e0, &
          0.0000000000000000432e0  ]
!                                ,
! series for g1   on the interval  2.00000e-02 to  6.25000e-02
!                                        with weighted error   5.48e-17
!                                         log weighted error  16.26
!                               significant figures required  15.47
!                                    decimal places required  16.92
!
    real, save :: g1cs(21) = [     &
         -0.3040578798253495954e0, &
         -0.0566890984597120588e0, &
          0.0039046158173275644e0, &
         -0.0003746075959202261e0, &
          0.0000435431556559844e0, &
         -0.0000057417294453025e0, &
          0.0000008282552104503e0, &
         -0.0000001278245892595e0, &
          0.0000000207978352949e0, &
         -0.0000000035313205922e0, &
          0.0000000006210824236e0, &
         -0.0000000001125215474e0, &
          0.0000000000209088918e0, &
         -0.0000000000039715832e0, &
          0.0000000000007690431e0, &
         -0.0000000000001514697e0, &
          0.0000000000000302892e0, &
         -0.0000000000000061400e0, &
          0.0000000000000012601e0, &
         -0.0000000000000002615e0, &
          0.0000000000000000548e0  ]
!
! series for g2   on the interval  0.00000e+00 to  2.00000e-02
!                                        with weighted error   5.01e-17
!                                         log weighted error  16.30
!                               significant figures required  15.12
!                                    decimal places required  17.07
!
    real, save :: g2cs(34) = [     &
         -0.0967329367532432218e0, &
         -0.0452077907957459871e0, &
          0.0028190005352706523e0, &
         -0.0002899167740759160e0, &
          0.0000407444664601121e0, &
         -0.0000071056382192354e0, &
          0.0000014534723163019e0, &
         -0.0000003364116512503e0, &
          0.0000000859774367886e0, &
         -0.0000000238437656302e0, &
          0.0000000070831906340e0, &
         -0.0000000022318068154e0, &
          0.0000000007401087359e0, &
         -0.0000000002567171162e0, &
          0.0000000000926707021e0, &
         -0.0000000000346693311e0, &
          0.0000000000133950573e0, &
         -0.0000000000053290754e0, &
          0.0000000000021775312e0, &
         -0.0000000000009118621e0, &
          0.0000000000003905864e0, &
         -0.0000000000001708459e0, &
          0.0000000000000762015e0, &
         -0.0000000000000346151e0, &
          0.0000000000000159996e0, &
         -0.0000000000000075213e0, &
          0.0000000000000035970e0, &
         -0.0000000000000017530e0, &
          0.0000000000000008738e0, &
         -0.0000000000000004487e0, &
          0.0000000000000002397e0, &
         -0.0000000000000001347e0, &
          0.0000000000000000801e0, &
         -0.0000000000000000501e0  ]

    integer, save :: nf1   = 0
    integer, save :: nf2   = 0
    integer, save :: ng1   = 0
    integer, save :: ng2   = 0

    real, save    :: xbnd  = 0.0
    real, save    :: xbig  = 0.0
    real, save    :: xmaxf = 0.0
    real, save    :: xmaxg = 0.0

    real          :: tol
    !
    ! Initialisation
    !
    if ( nf1 == 0 ) then
        tol   = 0.1*r1mach(3)
        nf1   = inits (f1cs, 20, tol)
        nf2   = inits (f2cs, 29, tol)
        ng1   = inits (g1cs, 21, tol)
        ng2   = inits (g2cs, 34, tol)

        xbig  = sqrt (1.0/r1mach(3))
        xmaxf = exp (amin1(-alog(r1mach(1)), alog(r1mach(2))) - 0.01)
        xmaxg = 1.0/sqrt(r1mach(1))
        xbnd  = sqrt(50.0)
    endif

    !
    ! Argument ranges:
    !      x < 4.0
    !      4.0 <= x <= xbig
    !      x > xbig
    !
    if ( x < 4.0 ) then
        f = ieee_value( x, ieee_quiet_nan )
        g = ieee_value( x, ieee_quiet_nan )
    else
        if ( x <= xbnd ) then
            f = (1.0 + csevl ((1.0/x**2-0.04125)/0.02125, f1cs, nf1))/x
            g = (1.0 + csevl ((1.0/x**2-0.04125)/0.02125, g1cs, ng1))/x**2

        elseif ( x <= xbig ) then
            f = (1.0 + csevl (100./x**2-1., f2cs, nf2))/x
            g = (1.0 + csevl (100./x**2-1., g2cs, ng2))/x**2

        else
            f = 0.0
            if ( x < xmaxf ) f = 1.0/x
            g = 0.0
            if ( x < xmaxg ) g = 1.0/x**2
        endif
    endif
end subroutine r9sifg

! en --
!     Evaluate the function En(x) via a recursion relation
!
real function integral_en( n, x )
    integer, intent(in) :: n
    real, intent(in)    :: x

    real                :: en
    integer             :: i

    if ( n < 1 ) then
        integral_en = ieee_value( x, ieee_quiet_nan )

    else
        en = e1(x)

        do i = 1,n-1
            en = (exp(-x) - x * en) / real(i)
        enddo
        integral_en = en
    endif
end function integral_en

end module fullerton_ei
