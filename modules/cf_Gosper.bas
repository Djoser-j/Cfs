' *****************************************************************************
'Subject: R. W. Gosper's continued fraction arithmetic.
'Ref.   : HAKMEM, A.I. Lab Memo #239, M.I.T., 1972
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.06.0 with GMP library 6.1.1

#include "cfr_lib.bi"
#inclib "cfr_math"

' *****************************************************************************
'Knuth exercises and Gosper's Appendix 2
sub appxtest ()
dim as integer t, t2, fl = (getsw and 4) <> 0
dim as cfa m, n, s
dim as cf a, b, c
dim g as poly, polroot as polr_t
dim sqrt2 as sqrt_2, cexp as exp_t
dim Csqrt as sqrtN_t, quad as quad_t
dim as pair p1, p2
dim d as long
m.ini: n.ini: s.ini
p1.ini: p2.ini

print : ? "appxtest" : ?

   readcoef g, "1,0,0,-2" '            Knuth's exercise 4.5.3.13
   polroot.ini g
   c.ini: empty(a)
   do : push c, polroot.f(a)
   loop until tm(c)
    outcf "x^3 - 2 = 0", c

   tset n, 0,97,0,39, 0,-62,0,-25
   readcf n.x, "-1,5,1,1,1,2,1,2" '    4.5.3.15
   dolf(c, n)
    outcf "ex 4.5.3.15", c

   print
   cvcf c, 254/100 '                   from Gosper's Appendix 2:
    outcf "centimeters per inch", c

   print
   tset n, 0,70,0,29, 0,12,0,5 '       Cf output
   cvcf n.x, 0
   t = getsw : setsw 4
   dolf(c, n) : setsw t
    outcf "terms:", c

   print
   readcf n.x, "5,1,4,1" '             Cf input
   tset n, 0,1,0,0, 0,0,0,1
   t = getsw : setsw 4
   do : loop until nextlfd(n) > Hold
   setsw t

   print
   'Throughput
   tset n, 0,0,0,2, 0,-1,0,3 '         2 / (3 - nx)
   c.ini: sqrt2.ini
   do
      push n.x, sqrt2.f
      push c, nextlfd(n)
   loop until tm(c)
    outcf "2 / (3 - sqrt(2))", c

   'A more interesting case
   tset m, 0,1,0,-1, 0,1,0,1
   tset n, 0,-4,0,4, 0,1,0,1
   cvcf a, 1: cexp.ini a
   a.ini: b.ini
   do
      push m.x, cexp.f
      d = nextlfd(m) '                 y:= (e - 1) / (e + 1)
      push a, d: push n.x, d
      push b, nextlfd(n) '             4 * (1 - y) / (1 + y)
   loop until tm(b)
   push a, Inf
    outcf "tanh(1/2)", a
    outcf "4 / e", b

   print
   'a slightly fancier example
   tset m, 1,0,0,1, 1,0,0,-1 '        (xy + 1) / (xy - 1)
   tset n, 2,1,0,0, 1,0,1,0
   cvcf a, 1: cexp.ini a
   a.ini: b.ini: c.ini
   quad.ini 1,0,6
   do
      d = cexp.f
      push m.x, d: push m.y, d
      d = nextbfd(m) '                 x = coth(1)
      push n.x, d: push a, d
      d = quad.f '                     y = sqrt(6)
      push n.y, d: push b, d
      push c, nextbfd(n)
   loop until tm(a)
    outcf "x", a: outcf "y", b
    outcf "(2xy + x) / (xy + y)", c

   print
   cvcf c, 8.31398 '                   gas constant limits
    outcf "Rl", c
   cvcf c, 8.31466
    outcf "Ru", c

   print
   'Square roots of Cfs,
   'warmup exercise
   c.ini: quad.ini 10,0,17 '          (Knuth 4.5.3.12)
   do : push c, quad.f '               sqrt with rational argument
   loop until tm(c)
    outcf "sqrt(17/10)", c

   'The Real Thing
    outcf "coth(1)", a '               compute coth(1/2) from coth(1)
   cvcf c, 2 '                         initial guess Y0 = 2
   n.x = a: t = 1
   while (t shr 2) < Cfx
      n.y = c: m.y = c: c.ini '        feed back
      tset n, 1,0,0,-1, 0,-1,1,0 '     operations
      tset m,  0,1,1,0, 0,0,0,2
      do '                             Newton iteration
         push m.x, nextbfd(n) '        mx:= (xy - 1) / (y - x)
         push c, nextbfd(m) '          mean (mx + y) / 2
      loop until tm(c) or (c.w > t)
      push c, Inf

      if fl then
         print "t_coth ="; t
         outcf "Y", c '                approximate coth(1/2)
      end if
      t shl= 1
   wend
   if fl = 0 then
    c.s = a.s
    outcf "coth(1/2)", c
   end if

   cvcf c, "5,1" '                     Newton sqrt with Cf-argument
   Csqrt.f c, c
    outcf "sqrt(6)", c

   print
   'Non-regular Cfs
   tset p1, 1, 1
   tset n, 0,1,0,0, 0,1,0,1
   tset m, 0,4,0,0, 0,0,0,1
   a.ini: b.ini: t = 0
   do
      d = gencf(n, p1, t, 1) '         atan(1)
      push a, d: push m.x, d
      push b, nextlfd(m) '                    * 4
   loop until tm(b)
    outcf "arctan(1)", a
    outcf "pi", b

   'Conversion to decimal: see sub outcf,
   'conversion from decimal: see sub cvcf

   print
   'Approximations
   cvcf a, 0.312 '                     Joe diMolisho, over 0.312
    outcf "bl", a
   cvcf b, 0.3125
    outcf "bu", b
   a.q(4) += 1
    outcf "s", a

   cvcf a, 0.3335 '                    Knuth 4.5.3.39, avg 0.334
    outcf "bl", a
   cvcf b, 0.3345
    outcf "bu", b
   b.q(5) = Inf: b.w = 5
    outcf "s", b

   'Continued Logarithms are skipped,
   'but here is an example of usual logarithms:

   'log2(10) = (ln(5/4) + 3ln(2)) / ln(2)
   'ln(x / y) = 2atanh((x - y) / (x + y))
   print
   tset p1, -1,9
   tset p2, -1,3
   tset n, 0,1,0,0, 0,9,0,1
   tset m, 0,1,0,0, 0,3,0,1
   tset s, 0,1,3,0, 0,0,1,0
   c.ini: t = 0: t2 = 0
   do
      push s.x, gencf(n, p1, t, 1) '   atanh(1/9)
      push s.y, gencf(m, p2, t2, 1) '  atanh(1/3)
      push c, nextbfd(s)
      if t2 =-1 then push c, Inf
   loop until tm(c)
    outcf "log2(10)", c '              OEIS A028232
   rate t: rate t2
   cf2q "appr", c, 4, 1

   clrs
end sub

' *****************************************************************************
'Gosper's HAKMEM constant (1972)
sub Gosperconst ()
dim as cfa m, n, s, p, q
dim pi as pi_t, e as exp_t
dim tanhr as tanhr_t, sin3 as sin3_t
dim Csqrt as sqrtN_t
dim c as cf, d as long

print : ? "Gosper's constant" : ?

   m.ini: tset m,  0,0,0,3, 1,0,0,0 '  101B, p.39
   n.ini: tset n,  0,1,1,0, 0,0,0,1
   s.ini: tset s, 0,1,-1,0, 0,0,0,1
   p.ini: tset p,  1,0,0,0, 0,0,0,1
   q.ini: tset q,  0,1,0,0, 0,0,1,0
   pi.ini 1
   cvcf c, "1": e.ini c
   tanhr.ini 5: sin3.ini 23

   c.ini
   do
      d = pi.f
      push m.x, d: push m.y, d
      push n.x, nextbfd(m)
      push n.y, e.f
      push q.x, nextbfd(n) '           x = 3/ pi^2 + e

      push s.x, tanhr.f
      push s.y, sin3.f
      d = nextbfd(s) '                 y = tanh(sqrt(5)) - sin(3*23)
      push p.x, d: push p.y, d
      push q.y, nextbfd(p) '           y^2

      push c, nextbfd(q) '             quotient
   loop until tm(c)
   Csqrt.f c, c '                      sqrt(x / y^2)
    outcf "sqrt(3/pi^2 + e) / (tanh(sqrt(5)) - sin(69))", c

   if 1 then
      rate pi.t: rate pi.t2
      rate e.t
      rate tanhr.t: rate tanhr.cquad.t
      rate sin3.t
   end if

   'Long partial quotients:
   'Cfx 100001, 99999 correct digits, 94 sec.

   clrs
end sub

' *****************************************************************************
'Gosper's numbers in HAKMEM 239, items 97-101
sub hakmemcf ()
dim as cf a, b, c
dim n as cfa, cexp as exp_t
dim pi as pi_t, sqrt2 as sqrt_2
dim csin as sin_t, ctanh as tanh_t
dim quad as quad_t, besseli as besseli_t
dim g as string, t as integer
dim d as long, p7 as pair

print : ? "HAKMEMcf" : ?

   n.ini: tset n, 1,0,0,0, 0,0,0,1 '   Item 97, p.36
   c.ini: sqrt2.ini
   do
      d = sqrt2.f
      push n.x, d: push n.y, d
      if sqrt2.t = 1 then
         n.x.q(0) -= 1: n.y.q(0) += 1
      end if
      push c, nextbfd(n)
   loop until tm(c)
    outcf "(sqrt(2) - 1)(sqrt(2) + 1)", c

   cvcf c, "5826625309/12621792855" '  Item 98
    outcf "psi(1+x)~0", c

   print
   c.ini
   besseli.ini 1, 2,3
   do : push c, besseli.f '            Item 99
   loop until tm(c) '                  integer order only
    outcf "I-Bessel_1(2/3)/_0(2/3)", c

   c.ini
   besseli.ini 0, 2,1
   do : push c, besseli.f
   loop until tm(c)
    outcf "I-Bessel_0(2)/_1(2)", c

   print
   cvcf c, 3.14159265346744 '          101A 0), p.37
    outcf "~pi", c

   cvcf c, 10000/254 '                 1)
    outcf "inches/meter", c

   print
   tset n, 0,2,0,0, 0,0,0,1 '          2)
   cvcf c, 1: cexp.ini c
   c.ini
   do
     push n.x, cexp.f
     push c, nextlfd(n)
   loop until tm(c)
    outcf "2e", c

   tset n, 0,4,0,-2, 0,1,0,-1
   cvcf a, 2/3: cexp.ini a
   a.ini: b.ini
   do
      d = cexp.f
      push a, d: push n.x, d
      push b, nextlfd(n)
   loop until tm(b)
    outcf "exp(2/3)", a
    outcf "(4*exp(2/3) - 2) / (exp(2/3) - 1)", b

   print
   'simple programs which produce terms on demand
   cvcf a, "1"
   c.ini: cexp.ini a
   do : push c, cexp.f
   loop until tm(c)
    outcf "e", c

   c.ini: pi.ini 1
   do : push c, pi.f
   loop until tm(c)
    outcf "pi", c

   c.ini: sqrt2.ini
   do : push c, sqrt2.f
   loop until tm(c)
    outcf "sqrt(2)", c

   cvcf a, 0.5
   c.ini: csin.ini a
   do : push c, csin.f
   loop until tm(c)
    outcf "sin(0.5)", c

   p7.ini: tset p7, 7, 1
   tset n, 0,7,0,0, 0,1,0,1
   c.ini: t = 0
   do : push c, gencf(n, p7, t, 1)
   loop until tm(c)
    outcf "sqrt(7)*atan(sqrt(7))", c

   print ' Best truncations only,
   'see Appendix 2, `Approximations'
   cvcf c, 3.141592653012 '            3)
    cf2q "~pi", c, 3, 1

   print
   cvcf a, 113/36 '                    4), p.38
   cvcf b, 355/113
   d = comp(t, a, b)
   outcf "a", a
   outcf "b", b
   select case d
   case -1: print "a < b";
   case 0: print "a = b";
   case 1: print "a > b";
   end select
   print ", diff.term:"; t

   print
   cvcf c, "1/0" '                     5)
    outcf "oo", c

   cvcf c, "2,3,4,5" '                 6)
    outcf "+", c
   Cneg c
    outcf " ", c

   clrs
   Gosperconst

   print
   g = "7, 5, 1, 0,-1,-5,-1, 9" '     p.42
   cvcf c, g: Cidem c
    outcf "["+g+"] equivalent to", c

   print
   tset n, 0,29,0,12, 0,12,0,5 '       p.44, ad - bc = 1
   cvcf n.x, "1,3,5,7" '               unmodified output
   dolf(c, n)
    outcf "tail 1,3,5,7:", c

   cvcf a, 1/69
   c.ini: ctanh.ini a '                Hurwitz numbers
   do : push c, ctanh.f
   loop until tm(c)
    Cinv c: outcf "coth(1/69)", c

   c.ini: quad.ini 1,0,105
   do : push c, quad.f
   loop until tm(c)
    outcf "sqrt(105)", c

   'item 101C : see Appendix 2, `Joe diMolisho'

   clrs
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

   appxtest
   hakmemcf

print : ? "timer:"; csng(timer - tim); "s"
system
