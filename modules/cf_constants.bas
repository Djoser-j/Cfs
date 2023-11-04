' *****************************************************************************
'Subject: Generate mathematical constants from continued fractions.
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#include "cfr_lib.bi"
#inclib "cfr_math"

' *****************************************************************************
' OEIS A016730
sub log2 (byval k as long = 2)
dim a as cf, clog as log_t

   print : cvcf a, k
   clog.ini a: a.ini
   do : push a, clog.f
   loop until tm(a)
    outcf "ln(2)", a

   clrs
end sub

' *****************************************************************************
' OEIS A053002
'Gauss (1799)
sub Gauss ()
dim as cf a, b
dim sqrt2 as sqrt_2
dim Cagm as agm_t

print : ? "Gauss's constant" : ?

   cvcf a, 1
   b.ini: sqrt2.ini
   do : push b, sqrt2.f
   loop until tm(b)
   Cagm.f a, a, b: Cinv a '            a:= 1/agm(1,sqrt(2))
    outcf "1/agm(1,sqrt(2))", a

   clrs
end sub

' *****************************************************************************
'gamma(1/4) = sqrt(2pi^(3/2) / agm(1,sqrt(2)))
sub gamma025 ()
dim as cfa m, n
dim as cf a, b, c
dim cpi as pi_t
dim sqrt2 as sqrt_2
dim Csqrt as sqrtN_t
dim Cagm as agm_t
dim d as long
dim prc as cnv = cnvmax
m.ini: n.ini

print : ? "Gauss and the agM" : ?

   tset n, 0,2,0,0, 0,0,0,1
   c.ini: cpi.ini 1
   do
      push n.x, cpi.f
      push c, nextlfd(n) '             2pi
   loop until tm(c)
   Csqrt.f n.x, c
   least prc, n.x
    outcf "sqrt(2pi)", n.x '           OEIS A058293

   cvcf a, 1
   b.ini: sqrt2.ini
   do : push b, sqrt2.f
   loop until tm(b)
   Cagm.f m.x, a, b '                  1 / Gauss's constant
   least prc, m.x
    outcf "agm(1,sqrt(2))", m.x '      1 / A053002

   m.y = c
   tset m, 0,0,1,0, 0,2,0,0
   tset n, 2,0,0,0, 0,0,0,1
   b.ini: c.ini
   do
      d = nextbfd(m) '                 pi/ agm(1,sqrt(2))
      push b, d
      push n.y, d
      push c, nextbfd(n) '                               * 2sqrt(2pi)
   loop until tm(b) and tm(c)
   least b, prc
    outcf "lemniscate constant", b '   A062540

   Csqrt.f c, c
   least c, prc
    outcf "gamma(1/4)", c '            A068153

   clrs
end sub

' *****************************************************************************
' OEIS A019810 (decimal expansion)
sub sinpi (byval r as long = 1, byval s as long = 180)
dim c as cf, csinpi as sinpi_t

   print
   c.ini: csinpi.ini r, s
   do : push c, csinpi.f
   loop until tm(c)
    outcf "sin(1)", c
   ''rate csinpi.t

   clrs
end sub

' *****************************************************************************
'iterative method for sin(1)
'Ptolemaeus (140), al-Kashi (1425), Raphson (1690)
sub sin1ter ()
dim as cfa m, n, p, q
dim as cf a, b, c
dim as poly f, g
dim as polr_t cpolroot1, cpolroot2
dim as sqrtN_t Csqrt
dim d as long, t as integer
dim prc as cnv
m.ini: n.ini: p.ini: q.ini

print : ? "sin(1) iteration" : ?  '    angles in degrees

   tset n, 0,1,-1,4, 0,0,0,8 '         sin(3)^2
   readcoef f, "1,0,-5,0,5"
   cpolroot1.ini f
   cvcf a, "-2,1"
   readcoef g, "1,0,-9,0,9"
   cpolroot2.ini g
   cvcf b, "2,1"
   c.ini
   do
      push n.x, cpolroot1.f(a) '       -sqrt((5 - sqrt(5)) / 2)
      push n.y, cpolroot2.f(b) '       sqrt((9 + 3sqrt(5)) / 2)
      push c, nextbfd(n) '             (4 - 4cos(36 - 30)) / 8
   loop until tm(c)
    ''outcf "-4sin(36)sin(30)", n.x
    ''outcf " 4cos(36)cos(30)", n.y
   Csqrt.f p.y, c
   prc = p.y.s
    ''outcf "sin(3)", p.y

   cvcf a, "0": t = 1 '                trisection c = 3x - 4x^3
   while (t shr 1) < Cfx
      m.x = a: m.y = a: n.x = a '      feed back
      tset m, 4,0,0,0,  0,0,0,1 '      operations
      tset n, 2,0,0,0,  0,0,0,1
      tset p, 0,1,-1,0, 0,0,0,1
      tset q, 0,1,0,0,  0,0,3,-3
      a.ini '                          Newton iteration
      do
         d = nextbfd(m) '              4x^2
         push q.y, d: push n.y, d
         push p.x, nextbfd(n) '        8x^3
         push q.x, nextbfd(p) '       (8x^3 - c)
         push a, nextbfd(q) '                   / (12x^2 - 3)
      loop until tm(a) or (a.w > t)
      push a, Inf

      if 0 then
         print "t_sin(1) ="; t
         outcf "X", a '                approximate sin(1)
      end if
      t shl= 1
   wend
    a.s = prc
    outcf "sin(1)", a

   clrs
end sub

' *****************************************************************************
' OEIS A049007
sub exphpi (byval r as long =-1, byval s as long = 2)
dim c as cf, exppi as exppi_t

   print : ? "i to the power i" : ?

   c.ini: exppi.ini r, s
   do : push c, exppi.f
   loop until tm(c)
    outcf "exp(-pi/2)", c
   ''rate exppi.t

   clrs
end sub

' *****************************************************************************
' OEIS A019474
'solve Cf_w * exp(Cf_w) = 1
sub omega ()
dim as integer fl = (getsw and 2) <> 0
dim as integer t = 1
dim n as cfa, w as cf
dim cexp as exp_t

print : ? "solve w*exp(w) = 1" : ?

cvcf w, 0.6
while (t shr 1) < Cfx
   n.ini: n.x = w
   tset n, 0,1,0,1, 0,0,1,1
   cexp.ini w: w.ini '                 feed back
   do '                                Newton iteration
      push n.y, cexp.f '              (w + 1)
      push w, nextbfd(n) '                   / (exp(w) + 1)
   loop until tm(w) or (w.w = t + 3)
   push w, Inf

   if fl then
      print "omega_t ="; t
      outcf " X", w
   end if
   t shl= 1
   clrs
wend
 outcf "omega", w
end sub

' *****************************************************************************
' OEIS A002852
'Euler (1734)
'Ref.: Brent & McMillan, High-Precision Computation of Euler's Constant,
'Mathematics of Computation, 34(149), 1980, pp.305-312
sub Euler ()
dim as mpz_ptr a, b, u, v, tmp
dim as cfa m=cfa(0)
dim as cf c
dim as log_t clog
dim as ulong k, n, n2
dim as integer t
dim prc as cnv

print : ? "Euler's constant" : ?

n = int(Cfx * 0.6)
if n > 65535 then
   print "overflow gamma, n = "; n
   exit sub
end if
n2 = n * n

cvcf c, 1, n
clog.ini c: c.ini
do : push c, clog.f '                  compute -log(n)
loop until tm(c)
 prc = c.s
'' outcf "-log("+str(n)+")", c
'' rate clog.t: rate clog.t2
''exit sub

m.ini: cf2mpq m, c '                   Cf recurrence

tmp = m.u(0)
a = m.u(1): b = m.v(1)
u = m.u(3): v = m.v(3)

k = mpz_sizeinbase(b, 2) + 1 '         resize
mpz_mul_2exp  a, a, k
mpz_mul_2exp  b, b, k

''? "size_b "; k

mpz_set u, a '                         B&M, algorithm B1
mpz_set v, b

k = 0
mpz_set_si  tmp, 0
do
   mpz_add_ui  tmp, tmp, (k shl 1) + 1
   k += 1

   mpz_mul_ui  b, b, n2 '              B * n^2 / k^2
   mpz_tdiv_q  b, b, tmp

   mpz_mul_ui  a, a, n2 '             (A * n^2 / k
t = mpz_tdiv_q_ui(a, a, k)
   mpz_add  a, a, b '                             + B) / k
t = mpz_tdiv_q_ui(a, a, k)

   t = mpz_sgn(a) = 0
   t or= mpz_sgn(b) = 0
   if t then exit do '                 precision exhausted

   mpz_add  u, u, a '                  U += A
   mpz_add  v, v, b '                  V += B
loop

mpq2cf c, m, 3 '                       Euclidean division U / V
least c, prc
 outcf "-psi(1)", c
 ''rate k

   'Long partial quotients:
   'Cfx limit 97400, 97398 correct digits, 126 sec.
   'n = 97400 * 0.6 = 58440, clog.t2 = 65532

   clrs
end sub

' *****************************************************************************
' OEIS A014538
'Catalan (1865)
'Ref.: Fee, Computation of Catalan's constant using Ramanujan's formula,
'Proc. ISSAC '90, pp.157-160
sub Catalan ()
dim as mpz_ptr a, b, g, e
dim as cfa m=cfa(0)
dim as cf c
dim t as integer
dim as ulong d, k
m.ini

print : ? "Catalan's constant" : ?

a = m.u(0): g = m.u(1)
b = m.v(0): e = m.v(1)

k = int(Cfx * 3.5)

mpz_set_si  g, 1
mpz_mul_2exp  g, g, k '                1 / 2
mpz_set  a, g
mpz_set  b, g
mpz_set  e, g
mpz_mul_2exp  e, e, 1 '                unit

''? "size_g "; k

k = 0
do
   k += 1
   d = (k shl 1) + 1

   mpz_mul_ui  b, b, k '               B * k / (2k + 1)
t = mpz_tdiv_q_ui(b, b, d)

   mpz_mul_ui  a, a, k '              (A * k
   mpz_add  a, a, b '                        + B) / (2k + 1)
t = mpz_tdiv_q_ui(a, a, d)

   t = mpz_sgn(a) = 0
   if t then exit do '                 precision exhausted

   mpz_add  g, g, a '                  G += A
loop

mpq2cf c, m, 1 '                       Euclidean division G / E
 outcf "beta(2)", c
''rate k

   clrs
end sub

' *****************************************************************************
'evaluate 2(zeta(3) - 1) with Nesterenko's Cf (1996)
function zeta3cf (byref n as cfa, byref t as integer) as long
dim as ulong k, tr(1)
dim as longint h
dim as long d
zeta3cf = Inf

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   k = t shr 2
   select case t and 3
   case 0
   tr(0) = (k shl 1) + 2
   tr(1) = k * (k + 1)
   case 1
   tr(0) = (k shl 1) + 4
   tr(1) = (k + 1) * (k + 2)
   case 2
   tr(0) = (k shl 1) + 3
   tr(1) = (k + 1) * (k + 1)
   case 3

      h = clngint(k + 2) * (k + 2)
      if abs(h) > Ulng then
         print "overflow zeta3cf, ";
         print "t =";t;", id =";n.id
         t = -1: exit do
      end if

   tr(0) = (k shl 1) + 2
   tr(1) = culng(h)
   end select

   mpz_set_ui  n.u(0), tr(0)
   mpz_set_ui  n.u(2), tr(1)

   do_genr n, 0
loop
zeta3cf = d
end function

' OEIS A013631
'proven irrational in 1978 by R. Apéry
sub zeta3 ()
dim as cfa n=cfa(0)
dim as cf c
dim t as integer
n.ini

print : ? "Apéry's constant" : ?

   tset n, 0,5,0,2, 0,4,0,2 '          1 + (b d) / 2(f h)
   c.ini: t = 0
   do : push c, zeta3cf(n, t)
   loop until tm(c)
    outcf "zeta(3)", c
   ''rate t

   'Long partial quotients:
   'Cfx limit 137368, 137366 correct digits, 38 sec.
   't_max = 185355, (k + 2)^2 < maxlong 2^31 - 1

   clrs
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

   sinpi
   sin1ter
   exphpi
   omega
   Euler
   log2
   Gauss
   Catalan
   zeta3
   gamma025

print : ? "timer:"; csng(timer - tim); "s"
system
