# continued-fraction-math-library
Cfr is a FreeBASIC library for computing with continued fractions.  
  
  
The Cf length (thus the maximum precision) is fixed at compile time,  
also the number of mpz state registers (a b c d) / (e f g h).  
Both are set in include file modules\cfr_lib.bi  
  
  
### Contents of the Cfs packet
  
  
Makefiles are in the base directory  
(-W1ndows only, but easy to adapt):  
  
1_make_cfr_dll.bat  
2_make_cf_demos.bat  
3_run_all_demos.bat  
  
_make_one_demo.bat  
_run_one_demo.bat  
  
  
#### Cfs\workdir\  
   Will hold logfiles  
  
  
#### Cfs\library\  
  
cfr_arith.bas  
&ensp;R. W. Gosper's continued fraction arithmetic: a FreeBASIC library  
&emsp;using GMP functions for performing integer arithmetic  
  
cfr_math.bas  
   Generate algebraic and transcendental functions using cfr_arith  
  
libgmp-10.dll.a  
   GMP import library for dynamic linking  
  
libcfr_math.dll.a  
   cfr_math import library  
  
  
#### Cfs\modules\bin\  
  
libgmp-10.dll  
   GMP dynamic link library, 32-bit, build 6.1.1  
  
cfr_math.dll  
   cfr_math dynamic link library, 32-bit  
  
#### Cfs\modules\  
  
cfr_lib.bi  
   Include file for continued fraction arithmetic  
  
gmp_short.bi  
   Include file for the GMP functions called by cfr  
  
cf_first.bas  
   First steps: input conversion, convergents,  
   building Cf power and cube root functions  
  
cf_Gosper.bas  
   Knuth exercises and Gosper's Appendix 2,  
   Gosper's numbers in HAKMEM 239, items 97-101  
  
cf_constants.bas  
   Generate mathematical constants from continued fractions:  
      sin(1 degree)  
      i to the power i  
      Omega (w * exp(w) = 1)  
      Euler's constant  
      ln(2)  
      Gauss' constant  
      Catalan's constant  
      Apéry's constant zeta(3)  
      gamma(1/4)  
  
cf_ratarg.bas  
   Applications with (Gaussian) rational arguments:  
      sin((x + yi) / w)  
      exp(pi * (x + yi) / w)  
      principal value of Log(x + yi)  
      complete elliptic integrals with rational parameter m = x / y  
      Lambert's function: x / y = z, solve Cf_w * exp(Cf_w) = z  
      demo run for functions with rational arguments x / y in cfr_math  
  
cf_trans.bas  
   Demo run for functions with Cf arguments in cfr_math library  
  
cf_polroot  
   Continued fraction expansion of real polynomial roots  
  
cf_trees.bas  
   P. J. Potts' expression trees for transcendental functions.  
   Iterating the logistic map  
  
  
#### Copyright:  
          (C) 2020 Djoser.j.Spacher, All rights reserved  
  
#### License:  
          GNU General Public License, GPL  
  
          ______________________________________________  
  
Find HAKMEM paper (Gosper on continued fractions) [here](https://perl.plover.com/classes/cftalk/INFO/hakmem.html)  
  
Find 'Appendix 2' on continued fraction arithmetic [here](https://perl.plover.com/classes/cftalk/INFO/gosper.txt)  
  
2017-09-08 W1ndows port of the [GMP library](https://sourceforge.net/projects/mingw/files/MinGW/Base/gmp/gmp-6.1.2/)  
