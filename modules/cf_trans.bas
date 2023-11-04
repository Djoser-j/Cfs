' *****************************************************************************
'Subject: Continued fraction arithmetic demo.
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#include "cfr_lib.bi"
#inclib "cfr_math"

dim shared crnd as rnd_t '             pseudo-random Cf generator
crnd.seed 0

' *****************************************************************************
'run all transcendental functions with Cf arguments
sub cftrans ()
dim as cfa n
dim as cf a, b, c, e
dim cquad as quad_t
dim clog as log_t, cexp as exp_t
dim ctan as tan_t, ctanh as tanh_t
dim csin as sin_t, csinh as sinh_t
dim ccos as cos_t, ccosh as cosh_t
dim catan as atan_t, catanh as atanh_t
dim casin as asin_t, casinh as asinh_t
dim cacos as acos_t, cacosh as acosh_t
dim cpow as pow_t, sqrt as sqrt_t
dim gd as gd_t, gdinv as gdinv_t
dim t as long, i as integer

print : ? "Cf_trans";

#macro cfun(id, msg)
   c.ini: id.ini a
   do
      push c, id.f
   loop until tm(c)
   least c, a.s
    outcf msg, c
   if id.t > 0 then print "k ";id.k;", ";
   rate id.t
#endmacro

n.ini
for i = 1 to 10
   print : ? "____________________________________________"
   print "i ="; i : ?

if i = 1 then
   cvcf e, 80 / 81
else
   e = b
end if

a.ini
select case i
case 1
   cquad.ini 1,-1, 1
   do : push a, cquad.f
   loop until tm(a)
case 2
   crnd.ini '                          pseudo-random
   do : push a, crnd.f
   loop until tm(a) 'a.w = 46
   push a, Inf
case 3
   cquad.ini 29,0,8 '                  Knuth
   do : push a, cquad.f
   loop until tm(a)
case 4
   cquad.ini 43,214,307
   do : push a, cquad.f
   loop until tm(a)
case 5
   cquad.ini -25,52,26
   do : push a, cquad.f
   loop until tm(a)
case 6
   cvcf c, "1": ctan.ini c
   do : push a, ctan.f
   loop until tm(a)
   Cinv a
case 7
   cvcf c, "1": ctanh.ini c
   do : push a, ctanh.f
   loop until tm(a)
case 8
   cvcf a, 1 / 32
case 9
   cvcf a, -2323 / 3070
case else
   cvcf a, 11
end select

   ''Cinv a
   ''Cneg a
    outcf "cf(z)", a
   b = a

   clrs
   print : ? "____________________________________________"
   cfun(sqrt, "sqrt(z)")
   print
   a = c
   n.x = a: n.y = a
   tset n, 1,0,0,0, 0,0,0,1
   dobf(c, n)
   least c, a
    outcf "sqrt(z)^2", c

   print
   a = b
   n.x = a: n.y = a
   tset n, 1,0,0,0, 0,0,0,1
   dobf(c, n)
   least c, a
    outcf "z^2", c
   print
   a = c
   cfun(sqrt, "sqrt(z^2)")

   print : ? "____________________________________________"
   a = b
   cfun(cexp, "exp(z)")
   print
   a = c
   cfun(clog, "log(z)")
   print
   a = b
   cfun(clog, "log(z)")
   print
   a = c
   cfun(cexp, "exp(z)")

   print : ? "____________________________________________"
   cpow.f c, b, e
    outcf "z ^ e", c
   rate cpow.clog.t
   rate cpow.cexp.t
   print
   Cinv e
   cpow.f c, c, e
    outcf "c ^ -e", c
   rate cpow.clog.t
   rate cpow.cexp.t

   clrs

   print : ? "____________________________________________"
   a = b
   cfun(csin, "sin(z)")
   print
   casin.f c, c
    outcf "asin(z)", c
   rate casin.catan.t
   print
   casin.f a, b
    outcf "asin(z)", a
   rate casin.catan.t
   print
   cfun(csin, "sin(z)")

   print : ? "____________________________________________"
   a = b
   cfun(ccos, "cos(z)")
   print
   cacos.f c, c
    outcf "acos(z)", c
   rate cacos.catan.t
   print
   cacos.f a, b
    outcf "acos(z)", a
   rate cacos.catan.t
   print
   cfun(ccos, "cos(z)")

   print : ? "____________________________________________"
   a = b
   cfun(ctan, "tan(z)")
   print
   a = c
   cfun(catan, "atan(z)")
   print
   a = b
   cfun(catan, "atan(z)")
   print
   a = c
   cfun(ctan, "tan(z)")

   clrs

   print : ? "____________________________________________"
   a = b
   cfun(csinh, "sinh(z)")
   print
   casinh.f c, c
    outcf "asinh(z)", c
   rate casinh.csqrt.t
   rate casinh.clog.t
   print
   casinh.f a, b
    outcf "asinh(z)", a
   rate casinh.csqrt.t
   rate casinh.clog.t
   print
   cfun(csinh, "sinh(z)")

   print : ? "____________________________________________"
   a = b
   cfun(ccosh, "cosh(z)")
   print
   cacosh.f c, c
    outcf "acosh(z)", c
   rate cacosh.csqrt.t
   rate cacosh.clog.t
   print
   cacosh.f a, b
    outcf "acosh(z)", a
   rate cacosh.csqrt.t
   rate cacosh.clog.t
   print
   cfun(ccosh, "cosh(z)")

   print : ? "____________________________________________"
   a = b
   cfun(ctanh, "tanh(z)")
   print
   catanh.f c, c
    outcf "atanh(z)", c
   rate catanh.clog.t
   print
   catanh.f a, b
    outcf "atanh(z)", a
   rate catanh.clog.t
   print
   cfun(ctanh, "tanh(z)")

   print : ? "____________________________________________"
   a = b
   gd.f c, a
    outcf "gd(z)", c
   rate gd.cexp.t
   rate gd.catan.t
   print
   gdinv.f c, c
    outcf "gd^-1(z)", c
   rate gdinv.csin.t
   rate gdinv.clog.t
   print
   gdinv.f a, b
    outcf "gd^-1(z)", a
   rate gdinv.csin.t
   rate gdinv.clog.t
   print
   gd.f c, a
    outcf "gd(z)", c
   rate gd.cexp.t
   rate gd.catan.t

   clrs
next i
end sub

'partial quotients frequency count
sub pqfreq ()
const ct as integer = 7000000
dim c as cf
dim as integer d, i, r(65)

print : ? "pq frequency"

randomize 11
d = rnd * 1000000
print : ? "s "; d
crnd.seed d

for i = 1 to ct
   d = crnd.f

   select case d
   case is < 64
      r(d) += 1
   case is < 4096
      r(64) += 1
   case else
      r(65) += 1
   end select
next i
for i = 0 to 65
   print using " &: #.#####"; i; r(i) * 10 / ct
next i

for i = 1 to 10
   c.ini: crnd.ini
   do
      push c, crnd.f
   loop until c.w = 35
   push c, Inf
    outcf "c", c
next i
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

cftrans
pqfreq

print : ? "timer:"; csng(timer - tim); "s"
system
