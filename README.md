# specfunc_fullerton

Library for evaluating special mathematical functions, based on the library at [netlib.org](http://netlib.org)

It currently offers single-precision implementations of the following functions:

* Airy functions: airy_ai, airy_bi, airy_aiprime, airy_biprime
* Modified Bessel functions: bessel_i0, bessel_i1, bessel_in, bessel_k0, bessel_k1, bessel_kn
* Beta functions: beta, beta_inc, log_beta
* Exponential and logarithmic integrals: integral_ei, integral_e1, integral_en, integral_li, integral_li, integral_si, integral_ci, integral_shi, integral_chi
* Gamma functions: fgamma, log_gamma, sign_log_gamma, inc_gamma, inc_gamma_compl, gamma_tricomi, reciprocal_gamma, digamma
* Inverse cosine and cosine hyperbolic integrals: cin, cinh
* Miscellaneous functions: cbrt, lnrel, exprel, integral_dawson, dawson_prime, integral_spence
* Pochhammer symbols: poc, poch1, factorial

The overall module is `specfunc_fullerton`, which you can use to get access to the entire library. Various specific
modules are available, but mostly for convenience during the development.

The test programs produce output files showing the results of the functions for various argument ranges. Special attention
is paid to the "breakpoints" in the code where a different algorithm is used. The idea is that the functions should be
continuous, so that on either side of the breakpoint the values as calculated should not differ too much.


# Development

The library by W. Fullerton has been around on netlib for a long time and I decided to modernise the code:

* Use modern constructs, starting with free-form source code and `implicit none`, eliminating the `goto` statements and so on
* Remove the specific implementations of standard functions like `log` and `sin`. These were introduced to guarantee platform-independence.
* Eliminate the error and warning messages. Instead, either `NaN` is produced or the warning about the inaccuracy of the result is suppressed.
  This is a rather radical change in some respects and perhaps a flag should be set to report such a situation. Alternatively:
  an auxiliary routine to check that this might occur.
* Test the results after each transformation against the results with the original code.
* Turn the library into a set of modules, with one overall module, `specfunc_fullerton`.

Not all potentially useful functions are made public, an example being the scaled versions of the Bessel functions.
(Since the Bessel functions of the first kind are part of the Fortran standard, only the functions of the second kind have
been included.)

*Note:* For the moment, the library only contains the single-precision versions. The double-precision versions and the complex
versions that are also part of the original library will be added in the near future.

# Documentation

The `doc` directory contains a PDF file documenting the functions in some detail. As mathematical expressions are involved
to precisely define these functions, the documentation source is in Latex.
