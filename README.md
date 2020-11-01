# continued-fraction-math-library  
  
Cfr is a FreeBASIC library for computing with continued fractions.  
Unpack to the base directory of your FreeBasic installation.  
  
The Cf length (thus the maximum precision) is fixed at compile time,  
also the number of mpz state registers (a b c d) / (e f g h).  
Both are set in include file modules\cfr_lib.bi  
  
  
### Contents of the Cfs packet  
  
  
#### Cfs\library\  
  
cfr_arith.bas  
&emsp; R. W. Gosper's continued fraction arithmetic: a FreeBASIC library  
&emsp; using GMP functions for performing integer arithmetic  
  
cfr_math.bas  
&emsp; Generate algebraic and transcendental functions using cfr_arith  
  
libgmp-10.dll.a  
&emsp; GMP import library for dynamic linking  
  
libcfr_math.dll.a  
&emsp; cfr_math import library  
  
  
#### Cfs\modules\bin\  
  
libgmp-10.dll  
&emsp; GMP dynamic link library, 32-bit, build 6.1.1  
  
cfr_math.dll  
&emsp; cfr_math dynamic link library, 32-bit  
  
#### Cfs\modules\  
  
cfr_lib.bi  
&emsp; Include file for continued fraction arithmetic  
  
gmp_short.bi  
&emsp; Include file for the GMP functions called by cfr  
  
cf_first.bas  
&emsp; First steps: input conversion, convergents,  
&emsp; building Cf power and cube root functions  
  
cf_Gosper.bas  
&emsp; Knuth exercises and Gosper's Appendix 2,  
&emsp; Gosper's numbers in HAKMEM 239, items 97-101  
  
cf_constants.bas  
&emsp; Generate mathematical constants from continued fractions:  
&emsp; &emsp; sin(1 degree)  
&emsp; &emsp; i to the power i  
&emsp; &emsp; Omega (w * exp(w) = 1)  
&emsp; &emsp; Euler's constant  
&emsp; &emsp; ln(2)  
&emsp; &emsp; Gauss' constant  
&emsp; &emsp; Catalan's constant  
&emsp; &emsp; Apéry's constant zeta(3)  
&emsp; &emsp; gamma(1/4)  
  
cf_ratarg.bas  
&emsp; Applications with (Gaussian) rational arguments:  
&emsp; &emsp; sin((x + yi) / w)  
&emsp; &emsp; exp(pi * (x + yi) / w)  
&emsp; &emsp; principal value of Log(x + yi)  
&emsp; &emsp; complete elliptic integrals with rational parameter m = x / y  
&emsp; &emsp; Lambert's function: x / y = z, solve Cf_w * exp(Cf_w) = z  
&emsp; &emsp; demo run for functions with rational arguments x / y in cfr_math  
  
cf_trans.bas  
&emsp; Demo run for functions with Cf arguments in cfr_math library  
  
cf_polroot  
&emsp; Continued fraction expansion of real polynomial roots  
  
cf_trees.bas  
&emsp; P. J. Potts' expression trees for transcendental functions.  
&emsp; Iterating the logistic map  
  
  
#### Cfs\workdir\  
&emsp; Will hold logfiles  

Makefiles are in the **base directory**  
(-W1ndows only, but easy to adapt):  
  
1_make_cfr_dll.bat  
2_make_cf_demos.bat  
3_run_all_demos.bat  
  
_make_one_demo.bat  
_run_one_demo.bat  
  
&nbsp;  
  
#### Copyright:  
&emsp;&emsp;&emsp;&emsp;&emsp; (C) 2020 Djoser.j.Spacher, All rights reserved  
  
#### License:  
&emsp;&emsp;&emsp;&emsp;&emsp; GNU General Public License, GPL  
  
&emsp;&emsp;&emsp; ______________________________________________  
  
Find HAKMEM paper (Gosper on continued fractions) [here](https://perl.plover.com/classes/cftalk/INFO/hakmem.html)  
  
Find 'Appendix 2' on continued fraction arithmetic [here](https://perl.plover.com/classes/cftalk/INFO/gosper.txt)  
  
2017-09-08 W1ndows port of the [GMP library](https://sourceforge.net/projects/mingw/files/MinGW/Base/gmp/gmp-6.1.2/)  
