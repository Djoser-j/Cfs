Continued fraction math library  
 =============================
Cfs is a FreeBasic library for computing with continued fractions.  
Unpack to the base directory of your FreeBasic installation.  
  
The Cf length (thus the maximum precision) is fixed at compile time,  
also the number of mpz state registers (a b c d) / (e f g h).  
Both are set in include file \modules\cfr_lib.bi  
The same file laconically doubles as library documentation.  
  
  
### Contents of the Cfs packet  
  
  
Makefiles are in the **base directory**  
  
Windows  
  
fbc_path.bat  
make_cfr_dll.bat  
1_make_cf_demos.bat  
2_run_all_demos.bat  
   
Linux  
  
Requires GMP library libgmp10 and developers tools libgmp-dev to be installed.  
  
make_cfr_so.sh  
1_make_cf_demos.sh  
2_installso.sh  
3_run_all_demos.sh  
  
A few additional scripts are in **Cfs\library\more\\**  
  
  
#### Cfs\library\  
  
cfr_arith.bas  
  R. W. Gosper's continued fraction arithmetic realized in FreeBasic  
  
cfr_math.bas  
  Generate algebraic and transcendental functions using cfr_arith  
  
  
#### Cfs\modules\bin\  
  
libgmp-32.dll  
  GMP dynamic link library, 32-bit, build 6.2.1  
  
libgmp-64.dll  
  GMP dynamic link library, 64-bit, build 6.2.1  
  
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
 (Compare the code with the papers to get acquainted with the lib.)  
  
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
    complete elliptic integrals with rational parameter m  
    Lambert's function: solve Cf_w * exp(Cf_w) = x / y  
    demo run for functions with rational arguments  
  
cf_trans.bas  
  Demo run for functions with Cf arguments in cfr_math  
  
cf_polroot  
  Continued fraction expansion of real polynomial roots  
  
cf_trees.bas  
  P. J. Potts' expression trees for transcendental functions  
  Iterating the logistic map  
  
  
#### Cfs\workdir\  
  Will hold logfiles  
  
   
  
#### Copyright:  
        (C) 2023 Djoser.j.Spacher, All rights reserved  
  
#### License:  
        GNU General Public License, GPL  
  
      ______________________________________________  
  
[HAKMEM paper](http://www.inwap.com/pdp10/hbaker/hakmem/cf.html)
(Gosper on continued fractions)  
  
His ['Appendix 2'](https://perl.plover.com/classes/cftalk/INFO/gosper.txt)
on continued fraction arithmetic  
  
[Potts' PhD Thesis](https://www.doc.ic.ac.uk/~ae/papers/potts-phd.pdf)
'Exact Real Arithmetic...' (pdf)  
  
[The FreeBASIC compiler](https://sourceforge.net/projects/fbc/files/)  
  
[Windows 32-bit and 64-bit ports of the GMP library](https://drive.google.com/file/d/1PNbZB-Ia7-YM7aI3PvmMo6fS2BtMA3hd/view?usp=sharing)  
