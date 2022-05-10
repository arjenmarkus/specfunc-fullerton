! specfunc_fullerton.f90 --
!     Library for evaluating special mathematical functions - overall module
!
!     Functions:
!     - Airy functions: airy_ai, airy_bi, airy_aiprime, airy_biprime
!     - Modified Bessel functions: bessel_i0, bessel_i1, bessel_in, bessel_k0, bessel_k1, bessel_kn
!     - Beta functions: beta, beta_inc, log_beta
!     - Exponential and logarithmic integrals: integral_ei, integral_e1, integral_en, integral_li, integral_li,
!           integral_si, integral_ci, integral_shi, integral_chi
!     - Gamma functions: fgamma, log_gamma, sign_log_gamma, inc_gamma, inc_gamma_compl, gamma_tricomi,
!           reciprocal_gamma, digamma
!     - Inverse cosine and cosine hyperbolic integrals: cin, cinh
!     - Miscellaneous functions: cbrt, lnrel, exprel, integral_dawson, dawson_prime, integral_spence
!     - Pochhammer symbols: poc, poch1, factorial
!
!     Note:
!     fgamma is an alternative for the standard gamma function, it is used internally to guarantee
!     identical results
!
!     Origin:
!     The code on netlib.org for the library by W. Fullerton was modernised to obtain this new
!     version.
!
!     Limitations:
!     For the moment only single-precision real functions are supported.
!     The complex and real double-precision versions are planned.
!
!     TODO:
!     Dawson and Spence functions
!
module specfunc_fullerton
    use fullerton_airy
    use fullerton_bessel
    use fullerton_beta
    use fullerton_ei
    use fullerton_gamma
    use fullerton_inv_ci
    use fullerton_misc
    use fullerton_poch
end module specfunc_fullerton
