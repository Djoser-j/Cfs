' *****************************************************************************
'Subject: Continued fraction arithmetic demo.
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#include "cfr_lib.bi"
#inclib "cfr_math"

'applications with (Gaussian) rational arguments
' *****************************************************************************
'principal value of Log(x + yi)
sub Cxlog (byval x as long, byval y as long)
dim clog as logr_t, cpi as pi_t
dim catan as atan_t, n as cfa
dim as cf a, b
dim g as string, sw as integer
dim d as long, s as longint

s = clngint(x) * x + clngint(y) * y
if s = 0 or (s shr 4) > Ulng then
   print "overflow Cxlog, s ="; s
   exit sub
end if

g = "log("+str(x)+" + i*"+str(y)+")"
cvcf a, s, 1: clog.ini a

sw = x < 0
sw or= x = 0 and y < 0
if sw then
   d = iif(y < 0,-1, 1)
   x = -x: y = -y
   n.ini: tset n, 0,1,d,0, 0,0,0,1 '   atan +/- pi
end if
cvcf a, y, x
catan.ini a

cpi.ini 1
a.ini: b.ini
do
   push a, clog.f '                    log(x^2 + y^2)/ 2

   d = catan.f '                       i*atan(y / x)
   if sw then
      push n.x, d
      push n.y, cpi.f
      d = nextbfd(n)
   end if
   push b, d
loop until tm(a) and tm(b)
 print g
 outcf "re_", a
 outcf "im_", b
clrs
end sub

'greatest common divisor
function gcd (byval a as integer, byval b as integer) as integer
dim as integer c
while b
   c = a mod b
   a = b: b = c
wend
gcd = abs(a)
end function

#macro gcdred(id)
   if w = 0 then exit sub
   if w < 0 then w = -w: x = -x: y = -y
   g = id + "(("+str(x)+" + i*"+str(y)+") / "+str(w)+")"
   u = y: y = w: v = w
   d = gcd(x, y)
   if d > 1 then x \= d: y \= d
   d = gcd(u, v)
   if d > 1 then u \= d: v \= d
#endmacro

'exp(pi * (x + yi) / w)
sub Cxexp (byval x as long, byval y as long, byval w as long)
dim cexppi as exppi_t
dim ccospi as cospi_t, csinpi as sinpi_t
dim as cfa m, n
dim as cf a, b
dim as string g
dim as long d, u, v

gcdred("exppi")

m.ini: tset m, 1,0,0,0, 0,0,0,1 '      real,
n.ini: tset n, 1,0,0,0, 0,0,0,1 '           imaginary parts
cexppi.ini x, y
ccospi.ini u, v: csinpi.ini u, v
a.ini: b.ini
do
   d = cexppi.f
   push m.x, d
   push m.y, ccospi.f
   push n.x, d
   push n.y, csinpi.f
   push a, nextbfd(m) '                exp(pi*x/y) * cos(pi*u/v)
   push b, nextbfd(n) '                exp(pi*x/y) * i*sin(pi*u/v)
loop until tm(a) and tm(b)
 print g
 outcf "re_", a
 outcf "im_", b
clrs
end sub

'sin((x + yi) / w)
sub Cxsin (byval x as long, byval y as long, byval w as long)
dim csinh as sinh_t, csin as sin_t
dim ccosh as cosh_t, ccos as cos_t
dim as cfa n
dim as cf a, b, c
dim as string g
dim as long d, u, v

gcdred("sin")

cvcf a, x, y
csin.ini a
ccos.ini a
cvcf b, u, v
ccosh.ini b
csinh.ini b

n.ini: c.ini
tset n, 1,0,0,0, 0,0,0,1
do
   push n.x, csin.f
   push n.y, ccosh.f
   push c, nextbfd(n)
loop until tm(c)
 print g
 outcf "re_", c

tset n, 1,0,0,0, 0,0,0,1
c.ini
do
   push n.x, ccos.f
   push n.y, csinh.f
   push c, nextbfd(n)
loop until tm(c)
 outcf "im_", c
clrs
end sub

' *****************************************************************************
'complete elliptic integrals with rational parameter m = x / y
sub Ellint (byval x as long, byval y as long)
dim as integer fl = (getsw and 2) <> 0
dim as integer i, j, sw, t = 1
dim g as string, cpi as pi_t
dim csqrt as sqrtN_t
dim as cfa m, n, p, q, r
dim as cf a, b, s1, s2
dim as long d
dim as cnv prc, prs
prc = cnvmax: prs = prc

sw = (x < 0 or y < 0)
if sw or (x > y) then '                parameter m = sin(fi)^2,
   print "illegal parameter m: ellint" '   modular angle fi
   exit sub '                          eccentricity k = sin(fi)
end if
if x = y then
   cvcf a, 1, 0
   cvcf b, 1, 1
    outcf "K(1)", a
    outcf "E(1)", b
   exit sub
end if

cvcf a, 1 '                            a:= 1
cvcf b, y - x, y '                     b:= cos(fi)^2 = 1 - m
cvcf s1, x, y shl 1 '                  c:= sin(fi)^2 / 2 = m / 2
cvcf s2, 0

m.ini: n.ini: r.ini
if x = 0 then goto skip
p.ini: q.ini
do
   m.x = a: n.x = a
   csqrt.f n.y, b
   m.y = n.y: least prc, n.y
   p.x = s1: q.x = s2
   tset m, 0,1,1,0, 0,0,0,2 '          arithmetic,
   tset n, 1,0,0,0, 0,0,0,1 '          geometric mean
   tset r, 1,0,0,0, 0,0,0,1
   tset p, 0,1,t,0, 0,0,0,1
   tset q, 0,1,t,0, 0,0,0,1
   a.ini: s1.ini
   b.ini: s2.ini
   do
      d = nextbfd(m) '                (a + sqrt(b)) / 2
      push a, d

      push r.x, d: push r.y, d
      push p.y, nextbfd(r) '           a ^ 2
      push s1, nextbfd(p) '            Sum(a_t ^ 2)

      d = nextbfd(n) '                 a * sqrt(b)
      push b, d

      push q.y, d
      push s2, nextbfd(q) '            Sum(b_t ^ 2)

      sw = tm(a) and tm(s1)
      sw and= tm(b) and tm(s2)
   loop until sw

   sw = comp(i, n.x, a)
   sw = abs(i - j) < 3: j = i '        convergence at i
   sw or= i > Cfx - 2 '                fully equal Cfs
   sw and= (t shl 1) > Cfx '           skip initial iterations

   least prc, a
   least prs, s1
   least prs, s2

   if fl then
      print "t, i:"; t;" ";i
      rawcf " (a + b)/2", a, 0
      rawcf " a*b", b, 0
      outcf " s1", s1
      outcf " s2", s2
   end if
   t shl= 1
loop until sw

skip:
m.x = a
r.x = s1: r.y = s2
tset m,  0,0,1,0, 0,1,0,0 '            K(m), A&S 17.6.3
tset r, 0,1,-1,0, 0,0,0,1
tset n, -1,1,0,0, 0,0,0,1 '            E(m), A&S 17.6.4
a.ini: b.ini: cpi.ini 2
do
   push m.y, cpi.f
   d = nextbfd(m) '                    pi/2 / a_t
   push a, d

   push n.x, d
   push n.y, nextbfd(r) '              r = Sum(a_t ^ 2) - Sum(b_t ^ 2)
   push b, nextbfd(n) '                K - K * r
loop until tm(a) and tm(b)
least a, prc
g = "("+str(x)+"/"+str(y)+")"
 outcf "K" + g, a
least b, prc
least b, prs
 outcf "E" + g, b
clrs
end sub

'x / y = z, solve Cf_w * exp(Cf_w) = z
sub Lambertw (byval x as long, byval y as long)
dim as integer fl = (getsw and 2) <> 0
dim as integer i, j, sw, t = 1
dim w as cf, cexp as exp_t
dim as cfa n, s, q
dim as double u, z
dim g as string

sw = (y = 0)
if sw = 0 then
   z = cdbl(x) / y
   sw = z < -0.36787944117 '          -1/e
end if
if sw then
   print " illegal argument lambertw"
   exit sub
elseif x = 0 then
   cvcf w, 0: goto skip
end if

if z < 3 then
   cvcf w, sgn(z) * 0.5
else
   u = log(z)
   cvcf w, u - log(u)
end if

do
   n.x = w: n.y = w: q.y = w '         feed back
   cexp.ini w: w.ini
   n.ini: tset n, 1,0,0,0, 0,0,0,1
   s.ini: tset s, y,0,0,x, 0,0,y,0
   q.ini: tset q, 0,1,0,0, 0,0,1,1
   do '                                Newton iteration
      push s.x, nextbfd(n) '          (w^2
      push s.y, cexp.f
      push q.x, nextbfd(s) '              + z * exp(-w))
      push w, nextbfd(q) '                              / w + 1
   loop until tm(w) or (w.w = t + 3)
   push w, Inf

   sw = comp(i, n.x, w)
   sw = abs(i - j) < 3: j = i '        convergence at i
   sw or= i > Cfx - 2 '                fully equal Cfs
   sw and= i > Cfx shr 1 '             skip initial iterations

   if fl then
      print "lambertw_t ="; t
      outcf " X", w
   end if
   t shl= 1
   clrs
loop until sw

skip:
g = "("+str(x)+"/"+str(y)+")"
 outcf "lambertw" + g, w
end sub

' *****************************************************************************
'run all transcendental functions with rational arguments
sub alltrans (byval r as integer)
dim as cf a, b, c
dim f as poly
dim polroot as polr_t
dim sqrtN as sqrtN_t
dim csqrt as sqrt_t, cquad as quad_t
dim clog as log_t, cexp as exp_t
dim ctanpi as tanpi_t
dim csinpi as sinpi_t
dim ccospi as cospi_t
dim ctan as tan_t, ctanh as tanh_t
dim csin as sin_t, csinh as sinh_t
dim ccos as cos_t, ccosh as cosh_t
dim catan as atan_t, catanh as atanh_t
dim casin as asin_t, casinh as asinh_t
dim cacos as acos_t, cacosh as acosh_t
dim gd as gd_t, gdinv as gdinv_t
dim cagm as agm_t, cpow as pow_t
dim cbesselj as besselj_t
dim cbesseli as besseli_t
dim prc as cnv
dim as integer p, q, t
dim as integer x =-1, u = 0
dim as integer y = r, v = 1
dim g as string, d as long

print : ? "alltrans"

#macro cfun(id, msg)
   c.ini: id.ini a
   do
      push c, id.f
   loop until tm(c)
   least c, a.s
    outcf msg, c
#endmacro

do
   print
   Cxlog x, y
   Cxexp x-y, 2*x-y, y
   Cxsin x-y, 2*x-y, y
   Ellint x, y
   Lambertw x, y

   cvcf a, x / y
   g = str(x)+"/"+str(y)

   cvcf b, 1
   cagm.f c, b, a
    outcf "agm(1,"+g+")", c

   g = "("+ g +")"

   '5 square root methods
   cfun(csqrt, "sqrt"+g)

   sqrtN.f c, a
    outcf "sqrtN"+g, c

   cvcf b, 1/2
   cpow.f c, a, b
    outcf g+"^(1/2)", c

   c.ini: cquad.ini y, 0, x
   do : push c, cquad.f
   loop until tm(c)
    outcf "qx^2 - p = 0", c

   readcoef f, str(y)+",0,"+str(-x)
   polroot.ini f
   c.ini: empty(b)
   do : push c, polroot.f(b)
   loop until tm(c)
    outcf "poly f(x) = 0", c


   cfun(clog, "log"+g)
   cfun(cexp, "exp"+g)

   c.ini: csinpi.ini x, y
   do : push c, csinpi.f
   loop until tm(c)
    outcf "sinpi" + g, c

   c.ini: ccospi.ini x, y
   do : push c, ccospi.f
   loop until tm(c)
    outcf "cospi" + g, c

   c.ini: ctanpi.ini x, y
   do : push c, ctanpi.f
   loop until tm(c)
    outcf "tanpi" + g, c

   cfun(csin, "sin"+g)
   cfun(ccos, "cos"+g)
   cfun(ctan, "tan"+g)

   clrs

   cfun(csinh, "sinh"+g)
   cfun(ccosh, "cosh"+g)
   cfun(ctanh, "tanh"+g)

   cfun(catan, "atan"+g)
   catanh.f c, a
    outcf "atanh"+g, c

   casin.f c, a
    outcf "asin"+g, c
   casinh.f c, a
    outcf "asinh"+g, c
   cacos.f c, a
    outcf "acos"+g, c
   cacosh.f c, a
    outcf "acosh"+g, c

   gd.f c, a
    outcf "gd"+g, c
   gdinv.f c, a
    outcf "gd^-1"+g, c

   clrs

   c.ini: d = 1
   cbesselj.ini d, x, y
   do : push c, cbesselj.f
   loop until tm(c)
    outcf "J-Bessel" + g, c

   c.ini
   cbesseli.ini d, x, y
   do : push c, cbesseli.f
   loop until tm(c)
    outcf "I-Bessel" + g, c

   p = sgn(x + .001) * y: q = abs(x)
   empty(b): if q then cvcf b, p / q
   g = "("+str(x)+"/"+str(y)+")^"
   g += "("+str(p)+"/"+str(q)+")"
   cpow.f c, a, b
    outcf g, c

   g = "("+str(p)+"/"+str(q)+")^"
   g += "("+str(x)+"/"+str(y)+")"
   cpow.f c, b, a
    outcf g, c

   clrs

   d = int((r + y) / v) '              Farey sequence
   t = d * u - x
   x = u: u = t
   t = d * v - y
   y = v: v = t

loop until x / y > 1.5

print
dim n as cfa
n.ini
cvcf a, 1/2
cvcf b, 1
cagm.f n.x, a, b
prc = n.x.s
tset n, 0,2,0,0, 0,0,0,1
dolf(c, n)
least c, prc
outcf "2*agm(1/2,1)", c

cvcf a, 1
cvcf b, 2
cagm.f c, a, b
outcf "    agm(1,2)", c

print
Ellint 80, 81

end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

alltrans 5 '                           order 5 Farey sequence

print : ? "timer:"; csng(timer - tim); "s"
system
