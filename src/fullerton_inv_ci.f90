! fullerton_inv_ci.f90 --
!     Module for the inverse cosine and cosh functions
!
module fullerton_inv_ci
    use ieee_arithmetic
    use fullerton_aux
    use fullerton_ei

    implicit none

    private
    public :: inverse_ci, inverse_chi

    interface inverse_ci
        module procedure cin
    end interface

    interface inverse_chi
        module procedure cinh
    end interface

contains

! cinh --
!     Original:
!     december 1980 edition.  w. fullerton, bell labs.
!
!     evaluate the hyperbolic cin function.
!             cinh(x) = integral from 0 to x of (cosh(t) - 1)/t dt.
!
real function cinh (x)
    real, intent(in) :: x

!
! series for cinh on the interval  0.00000e+00 to  9.00000e+00
!                                        with weighted error   6.64e-17
!                                         log weighted error  16.18
!                               significant figures required  15.08
!                                    decimal places required  16.68
!
    real, save :: cinhcs(10) = [  &
         0.1093291636520734431e0, &
         0.0573928847550379676e0, &
         0.0028095756978830353e0, &
         0.0000828780840721357e0, &
         0.0000016278596173914e0, &
         0.0000000227809519256e0, &
         0.0000000002384484842e0, &
         0.0000000000019360830e0, &
         0.0000000000000125454e0, &
         0.0000000000000000664e0  ]

    integer, save   :: ncinh = 0
    real, save      :: xsml  = 0.0
    real, save      :: xmin  = 0.0

    real            :: absx, y

    !
    ! Initialisation
    !
    if ( ncinh == 0 ) then
        ncinh = inits (cinhcs, 10, 0.1*r1mach(3))
        xsml = sqrt (r1mach(3))
        xmin = 2.0*sqrt(r1mach(1))
    endif

    !
    ! Argument ranges:
    !     abs(x) <= 3.0
    !     abs(x) >  3.0
    !
    absx = abs(x)

    if ( absx <= 3.0 ) then
        cinh = 0.0

        if ( absx > xmin ) then
            y = -1.0
            if ( absx > xsml ) then
                y = x*x/9.0 - 1.0
            endif

            cinh = x*x* (0.25 + csevl (y, cinhcs, ncinh))
        endif

    else
        cinh = integral_chi(absx) - euler - log(absx)
    endif
end function cinh

! cin --
!     Original:
!     december 1980 edition.  w. fullerton, bell labs.
!
real function cin (x)
    real, intent(in) :: x

!
! series for cin  on the interval  0.00000e+00 to  1.60000e+01
!                                        with weighted error   7.33e-17
!                                         log weighted error  16.14
!                               significant figures required  15.42
!                                    decimal places required  16.66
!
    real, save :: cincs(11) = [    &
         0.3707450175090968874e0, &
        -0.0589357489636444683e0, &
         0.0053818964211356912e0, &
        -0.0002986005284196214e0, &
         0.0000109557257532162e0, &
        -0.0000002840545487735e0, &
         0.0000000054697399488e0, &
        -0.0000000000812418746e0, &
         0.0000000000009586859e0, &
        -0.0000000000000092027e0, &
         0.0000000000000000733e0  ]

   integer, save :: ncin = 0
   real, save    :: xmin = 0.0

   real          :: absx, sinx, f, g

   !
   ! Initialisation
   !
   if ( ncin == 0 ) then
       ncin = inits (cincs, 11, 0.1*r1mach(3))
       xmin = sqrt (r1mach(1))
   endif

   !
   ! Argument ranges:
   !    x <= 0.0
   !    0.0 < x <= 4.0
   !    x > 4.0
   !
   cin  = 0.0
   absx = abs(x)

   if ( abs(x) > xmin ) then
       if ( x <= 4.0 ) then
           cin = x*x*csevl ((x*x-8.0)*.125, cincs, ncin)
       else
           call r9sifg (x, f, g)
           sinx = sin (absx)
           cin = -f*sinx + g*cos(absx) + log(absx) + euler
       endif
   endif
end function cin

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

end module fullerton_inv_ci
