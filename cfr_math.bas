' *****************************************************************************
'Subject: Generate algebraic and transcendental functions from Cfs.
'Refs.  : Abramowitz & Stegun, Handbook of Mathematical Functions, 1964
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#include once "cfr_lib.bi"


' *****************************************************************************
'declare and initialize prc
#macro initprc(prc, c)
   dim prc as cnv = cnvmax
   if qirr(c) then prc = c.s
#endmacro

'Newton square root and arithmetic-geometric mean
' *****************************************************************************
'Cf_r:= Cf_a ^ (1 / 2)
sub sqrtN_t.f (byref r as cf, byref a as cf) export
dim as integer fl = (getsw and 2) <> 0
dim as integer sw, t = 1
initprc(s, a)
if @a <> @r then r = a

if is0(r) then exit sub
with r
   if .q(0) < 0 then
      print " illegal argument sqrtN_t"
      empty(r): exit sub
   end if
   sw = (.q(0) = 0)
   if sw then Cinv r
   '
   n.y = r '                           radicand R > 1
  .q(0) = int(sqr(.q(0) - (.q(1) = 1)))
  .q(1) = Inf '                        initial guess X
end with
m.ini: n.ini

if fl then
   rawcf "R", n.y, 0
   rawcf "X0", r, 0
end if

while (t shr 1) < Cfx
   n.x = r: m.x = r: r.ini '           feed back
   tset n, 0,0,1,0, 0,1,0,0
   tset m, 0,1,1,0, 0,0,0,2
   do '                                Newton iteration
      push m.y, nextbfd(n) '           y:= R/ X
      push r, nextbfd(m) '             mean (X + y)/ 2
   loop until tm(r) or (r.w = t + 3)
   push r, Inf

   if fl then
      print "sqrtN_t ="; t
      rawcf " X", r, 0
   end if
   t shl= 1
wend
if sw then Cinv r
least r, s
off8 -= 2
end sub

'Cf_a:= agm(Cf_g, Cf_h)
sub agm_t.f (byref a as cf, byref g as cf, byref h as cf) export
dim as integer fl = (getsw and 2) <> 0
dim as integer i, j, sw, t = 1
initprc(s, g)
least s, h

if is0(g) or is0(h) then
   cvcf a, 0: exit sub
end if
if g.q(0) < 0 or h.q(0) < 0 then
   print " illegal argument agm_t"
   empty(a): exit sub
end if

m.ini: n.ini
if @g <> @a then a = g
b = h
do
   m.x = a: m.y = b
   n.x = a: n.y = b
   tset m, 0,1,1,0, 0,0,0,2 '          arithmetic,
   tset n, 1,0,0,0, 0,0,0,1 '          geometric mean
   a.ini: b.ini
   do
      push a, nextbfd(m)
      push b, nextbfd(n)
   loop until tm(a) and tm(b)

   sw = comp(i, n.x, a)
   sw = abs(i - j) < 3: j = i '        convergence at i
   sw or= i > Cfx - 2 '                fully equal Cfs
   sw and= (t shl 1) > Cfx '           skip initial iterations

   least s, a
   csqrt.f b, b
   least s, b

   if fl then
      print "agm_t ="; t
      rawcf " (a + b)/2", a, 0
      rawcf " sqrt(a*b)", b, 0
   end if
   t shl= 1
loop until sw
least a, s
off8 -= 2
end sub

'building transcendental functions
' *****************************************************************************
'print trafo
sub prnt_ (byref n as cfa, byval sw as integer, byval t as integer)
dim as zstring ptr s
dim as mpz_ptr tp(2) = {n.u(0), n.u(2), n.v(0)}
dim as string g = " <-tr("+str(t)+") "
dim as integer i, j = 1 - sw
for i = 0 to j
   s = allocate(mpz_sizeinbase(tp(i), 10) + 2)
   mpz_get_str  s, 10, tp(i)
   g += " "+*s
   g += iif(i = 1,","," ")
next i
if sw = 0 then g += " 1 "
print g+" 0"
end sub

'n.u(1) / n.v(1) = q, return integer(log_2(q))
'round down (sw = 0) or round half up (sw = 1)
function qsize (byref n as cfa, byval sw as integer) as integer
const alog2 as double = 1.44269504089
dim t(1) as const double = {0, 0.5}
dim as clong h, k
dim as double w = mpz_get_d_2exp(@k, n.u(1))
   w /= mpz_get_d_2exp(@h, n.v(1)): k -= h
   if k < 0 then return 0
   if k > mexp then return -1
   w = abs(w * pw2(k))
   if w < 1 then return 0
qsize = int(log(w) * alog2 + t(sw))
end function

#macro kerror(id)
if k = -1 then
   print " overflow "; id
   t = -1: exit sub
end if
#endmacro

'u / v:= u / 2v
sub halve overload (byref n as cfa)
with n
   if mpz_odd_p(.u(1)) then
      mpz_mul_2exp .v(1), .v(1), 1
   else
      mpz_tdiv_q_2exp .u(1), .u(1), 1
   end if
end with
end sub

'The below functions output a stream of Cf digits (lazy evaluation).
'Use func.ini to initialize member variables and start a new stream.
' *****************************************************************************
'solve ax^2 + bx - c = 0, a <> 0 and x >= 0
'ref: H. Cohen, Computational algebraic number theory, 5.7.2
sub quad_t.ini (byval a as long, byval b as long, byval c as long) export
dim as integer sw = (a = 0)
dim as double disc = cdbl(b) * b
   disc += cdbl(a) * c * 4
   sw or= (disc < 0 or disc > Flint)
   if sw then
      print " illegal argument quad_t"
      t = -1: exit sub
   end if

   i = 0: j = 1: t = 0
   rD = int(sqr(disc)) '               initialize
   p(j) = b: q(j) = a shl 1
   d = (rD - p(j)) \ q(j)
   p(i) = d * q(j) + p(j)
   disc -= cdbl(p(i)) * p(i)
   q(i) = clng(disc / q(j))
end sub

function quad_t.f as long export
if t = 0 then t = 1: return d '        first digit
if q(i) = 0 then t = -1 '              square_x: finite Cf
if t < 0 then return Inf
   swap i, j: t += 1
   d = (rD + p(j)) \ q(j) '            next digit
   p(i) = d * q(j) - p(j) '            recurrence
   q(i) -= d * (p(i) - p(j))
f = d
end function

'square root(2)
sub sqrt_2.ini export
t = 0
end sub

function sqrt_2.f as long export
f = 2 + (t = 0)
t += 1
end function

'pseudo-random partial quotients
sub rnd_t.seed (byval s as ushort) export
dim as integer d, i
' D = 9007199254710721
' h(D) = 1, period 273890412
if (s and 1) = 0 then
cgenr.ini 1,94906265,29615124
else
cgenr.ini 46791780,15789761,46791780
end if

for i = 1 to s
   d = cgenr.f
next i
end sub

sub rnd_t.ini export
with cgenr
  .t = (.p(.i) xor .q(.i)) and 2
  .t shr= 1: .d = 0
end with
end sub

function rnd_t.f as long export
f = cgenr.f
end function

'Elementary transcendental functions f(a) for Cf_a.
'The register pointers are re-assigned with
'each initialization to save mpz-resources.
'For Cfx larger than, say, 1000 these functions quickly slow down,
'unless the Cf-argument is a small rational.
' *****************************************************************************
'evaluate sqrt(1 + 2x/y)
function sqrtcf (byref n as cfa, byref q as pair,_
 byref t as integer) as long
dim as integer fl = (getsw and 4) <> 0
dim d as long

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   d = (t shl 1) + 1
   mpz_set_ui  n.u(0), 2 '             1  x  3x  x  5x 3x  7x 5x
   mpz_mul_ui  n.u(2), q.a, d '        -- -- -- --- -- --- -- ---...
   '                                   1- y+ 2+ 3y+ 2+ 5y+ 2+ 7y+
   if fl then prnt_ n, 0, t
   do_genr n, 0

   mpz_mul_ui  n.u(0), q.b, d '        (2t + 1)y

   mpz_mul_2exp  tmp, q.a, 1
   mpz_sub  n.u(2), n.u(2), tmp '      (2t - 1)x

   if fl then prnt_ n, 0, t
   do_genr n, 0
loop
sqrtcf = d
end function

sub sqrt_t.ini (byref a as cf) export
dim as long u0 = 0, u1 '               a >= 0
dim j as integer

   sw = (a.q(0) = 0): t = 0
   if a.q(0) < 0 then '                undefined
      print " illegal argument sqrt_t"
      t = -1
   end if
   if is0(a) then t = -1
   if t < 0 then exit sub

   m.ini: cf2mpq m, a
   if sw then swap m.u(1), m.v(1) '    Cf >= 1

   k = qsize(m, 1)
   kerror("sqrt_t")
   for j = 1 to k
      halve m : next

   with m
      mpz_sub .u(1),.u(1),.v(1)
      ' u:= x - y, v:= 2y              evaluate sqrt(1 + 2x/y)
      q.ini: halve m
      mpz_set  q.a,.u(1)
      mpz_set  q.b,.v(1)

      ' m:= 0,v,0,1, 0,v - u,0,1
      swap .u(1),.v(1)
      mpz_sub .v(1),.u(1),.v(1)
      mpz_set_ui .u(3), 1
   end with

   u1 = 1 shl (k \ 2)
   fl = (k and 1) = 1
   if fl then swap u0, u1 '            odd power k
   n.ini: sqrt2.ini
   tset n, u0,u1,0,0, 0,0,0,1 '        sqrt(a) * sqrt(2)^ k
end sub

function sqrt_t.f as long export
if sw then sw = 0: return 0 '          1 / sqrt(a)
if t < 0 then return Inf
d = sqrtcf(m, q, t)
if fl then
   push n.x, d
   push n.y, sqrt2.f
   d = nextbfd(n)
elseif k then
   push n.x, d
   d = nextlfd(n)
end if
f = d
end function

' *****************************************************************************
'generic Cf for transcendental functions,
'set switch to evaluate tanh (sw = 0) or atan (sw = 1)
function gencf (byref n as cfa, byref q as pair,_
 byref t as integer, byval sw as integer) as long export
dim as integer fl = (getsw and 4) <> 0
dim t2 as longint
dim d as long

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1
   '                                   x  x2  x2  x2  x2
   d = (t shl 1) + 1 '                 -- --- --- --- ---...   tanh A&S 4.5.70
   mpz_mul_ui  n.u(0), q.b, d '        y+ 3y+ 5y+ 7y+ 9y+

   if sw = 0 then '                    x  x2  4x2 9x2 16x2
      mpz_set  n.u(2), q.a '           -- --- --- --- ----...  atan A&S 4.4.43
   else '                              y+ 3y+ 5y+ 7y+  9y+
      t2 = clngint(t) * t
      if abs(t2) > Ulng then
         print "overflow gencf, ";
         print "t =";t;", id =";n.id
         t = -1: exit do
      end if
      mpz_mul_ui  n.u(2), q.a, t2
   end if

   if fl then prnt_ n, 0, t
   do_genr n, 0
loop
gencf = d
end function

'evaluate atan(1 / y) (sw = 0) or atanh(1 / y) (sw = 1)
function acotcf (byref n as cfa, byval y as long, byref t as integer,_
 byval sw as integer) as long
dim as integer i, fl = (getsw and 4) <> 0
dim as longint k(1)
dim as long d

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   d = (t shl 1) + 1
   k(0) = clngint(y) * d '             1   1   4   9  16
   k(1) = clngint(t) * t '             -- --- --- --- ---    atan(1 / y)
   mpz_set_ui  n.u(0), k(0) '          y+ 3y+ 5y+ 7y+ 9y+...
   mpz_set_ui  n.u(2), k(1)
   if sw then mpz_neg n.u(2), n.u(2) ' atanh(1 / y)

   for i = 0 to 1
      if abs(k(i)) > Ulng then
         print "overflow acotcf, ";
         print "t =";t;", id =";n.id
         t = -1: exit do
      end if
   next i

   if fl then prnt_ n, 0, t
   do_genr n, 0
loop
acotcf = d
end function

'logarithm and exponential
' *****************************************************************************
sub log_t.ini (byref a as cf) export
dim as integer j, sw '                 a > 0
   if leq0(a) then '                   undefined
      print " illegal argument log_t"
      t = -1: exit sub
   end if

   sw = (a.q(0) = 0)
   m.ini: cf2mpq m, a

   t = iif(fl, 1, 2)
   if sw then
      swap m.u(1), m.v(1): t = -t '    Cf >= 1
   end if
   k = qsize(m, 1)
   kerror("log_t")

   for j = 1 to k
     halve m : next

   with m
      q.ini
      ' u / v:= (x - y) / (x + y)      evaluate atanh(u / v)
      ' m:= 0,u,0,0, 0,v,0,1
      mpz_add  q.b, .u(1),.v(1)
      mpz_sub .u(1),.u(1),.v(1)

      mpz_mul  q.a, .u(1),.u(1)
      mpz_neg  q.a, q.a '             -u^2
      mpz_set .v(1), q.b '             v

      if mpz_sgn(q.a) = 0 then
         tset m, 0,k,0,0, 0,3,0,1 '    power of two
         mpz_set_si  q.a,-1
         mpz_set_ui  q.b, 3: k = 0
      end if

      if k then
         n.ini: s.ini
         tset n, 0,k,0,0, 0,3,0,1 '    + k * atanh(1/3)
         tset s, 0,t,t,0, 0,0,0,1 '                    * 2
      else
         if fl = 0 then
            mpz_mul_2exp .u(1),.u(1), 1
         end if
         if sw then mpz_neg .u(1),.u(1)
      end if
   end with
   t = 0: t2 = 0
end sub

function log_t.f as long export
if t < 0 then return Inf
d = gencf(m, q, t, 1)
if k then
   push s.y, d
   push s.x, acotcf(n, 3, t2, 1)
   if t2 = -1 then t = -1
   d = nextbfd(s)
end if
f = d
end function

constructor logr_t export
this.fl = -1
end constructor

sub exp_t.ini (byref a as cf) export
dim j as integer
   if a.q(0) = Inf then
      t = -1: exit sub
   end if

   m.ini: cf2mpq m, a '                      2x   x2  x2  x2
   '                                   1 + ------ --- --- ---
   k = qsize(m, 0) '                       (y-x)+ 3y+ 5y+ 7y+...
   kerror("exp_t")
   if k > 0 then k -= 1
   for j = 0 to k '                    evaluate exp(2x/y): offset 0
      halve m : next
   ' u:= x: v:= 2y

   with m
      q.ini
      mpz_mul  q.a,.u(1),.u(1) '       x^2
      mpz_set  q.b,.v(1) '             y

      ' m:= 0,u + v,0,1, 0,v - u,0,1
      mpz_add  tmp,.u(1),.v(1)
      mpz_sub .v(1),.v(1),.u(1)
      mpz_set .u(1), tmp
      mpz_set_ui .u(3), 1
   end with

   k -= 1: t = 0
   if k > -1 then
      redim n(k)
      for j = 0 to k '                 x = y = exp(z)
         n(j).ini
         tset n(j), 1,0,0,0, 0,0,0,1 ' exp(2z):= xy
      next j
   end if
end sub

function exp_t.f as long export
if t < 0 then return Inf
dim j as integer
d = gencf(m, q, t, 0)
for j = 0 to k
   push n(j).x, d
   push n(j).y, d
   d = nextbfd(n(j))
next j
f = d
end function

'hyperbolic and circular functions
' *****************************************************************************
sub tanh_t.ini (byref a as cf) export
dim j as integer
   if a.q(0) = Inf then
      t = -1: exit sub
   end if

   m.ini: cf2mpq m, a
   k = qsize(m, 0)
   kerror("tanh_t")

   if fl then
      if k > 0 then k -= 1
      s.ini: halve m '                 x = y = tanh a / 2b
      ' u:= x: v:= 2y
      ' m:= 0,u,0,0, 0,v,0,1
      if fl = 1 then '                 A&S 4.5.31
         tset s, 0,1,1,0, -i,0,0,1 '   sinh := (x + y) / (1 - xy)
      else '                           A&S 4.5.32
         tset s, i,0,0,1, -i,0,0,1 '   cosh := (1 + xy) / (1 - xy)
      end if
   end if

   for j = 1 to k
      halve m : next

   q.ini
   mpz_mul q.a, m.u(1), m.u(1) '       x^2
   mpz_set q.b, m.v(1) '               y
   if i = -1 then mpz_neg q.a, q.a '  -x^2 for tan

   k -= 1: t = 0
   if k > -1 then
      redim n(k)
      for j = 0 to k '                 x = y = tanh(z)
         n(j).ini '                    A&S 4.5.33
         tset n(j), 0,1,1,0, i,0,0,1 ' tanh(2z):= (x + y) / (1 + xy)
      next j
   end if
end sub

function tanh_t.f as long export
if t < 0 then return Inf
dim j as integer
d = gencf(m, q, t, 0)
for j = 0 to k
   push n(j).x, d
   push n(j).y, d
   d = nextbfd(n(j))
next j
if fl then
   push s.x, d
   push s.y, d
   d = nextbfd(s)
end if
f = d
end function

constructor sinh_t export
this.fl = 1
end constructor

constructor cosh_t export
this.fl = 2
end constructor

constructor tan_t export
this.i = -1
end constructor

constructor sin_t export
this.i = -1
end constructor

constructor cos_t export
this.i = -1
end constructor

'inverse circular functions
' *****************************************************************************
'pi with Machin's 1706 formula
'set k = 1, 2 or 4 for pi, pi/2, pi/4
sub pi_t.ini (byval k as long) export
   if k <> 2 then
      k = iif(k = 4, 1, 4)
   end if
   m.ini: n.ini: s.ini
   tset m,  0,4,0,0, 0,5,0,1 '         4*atan(1/5)
   tset n, 0,-1,0,0, 0,239,0,1 '                  - atan(1/239)
   tset s,  0,k,k,0, 0,0,0,1 '                                 * k
   t = 0: t2 = 0
end sub

function pi_t.f as long export
if t < 0 then return Inf
do
   push s.x, acotcf(m, 5, t, 0)
   push s.y, acotcf(n, 239, t2, 0)
   d = nextbfd(s)
loop until d > Hold
f = d
end function

sub atan_t.ini (byref a as cf) export
   t = 0
   if a.q(0) = Inf then '              return pi/2
      m.ini: q.ini
      tset m, 0,0,0,0, 0,1,0,1
      mpz_set_ui  q.a, 0
      mpz_set_ui  q.b, 1

      sw = -1: k = 2
      n.ini: cpi.ini k
      tset n, 0,1,1,0, 0,0,0,1

   else
dim as integer i, sg = (a.q(0) < 0)
dim as double z = iif(sg,-70/29, 29/70)
dim c as cf
      for i = 0 to 1 '                 select trafo
         cvcf c, z
         if comp(k, a, c) < 0 then exit for
         z = 1 / z
      next i
      sw = iif(sg,(i + 2) and 3, i)

      m.ini: cf2mpq m, a
      with m
         select case sw '              Möbius transformations
         case 1
         ' m:= 0,x - y,0,0, 0,y + x,0,1
            mpz_sub  tmp,.u(1),.v(1)
            mpz_add .v(1),.v(1),.u(1) ' S0
            mpz_set .u(1), tmp
         case 2
         ' m:= 0,-y,0,0, 0,x,0,1
            if mpz_sgn(.u(1)) > 0 then
               mpz_neg .v(1),.v(1) '    S-
            else
               mpz_abs .u(1),.u(1)
            end if
            swap .u(1),.v(1)
         case 3
         ' m:= 0,x + y,0,0, 0,y - x,0,1
            mpz_add  tmp,.u(1),.v(1)
            mpz_sub .v(1),.v(1),.u(1) ' Soo
            mpz_set .u(1), tmp
         end select

         q.ini
         mpz_mul  q.a,.u(1),.u(1) '    u^2
         mpz_set  q.b,.v(1) '          v
      end with

      n.ini
      d = iif(sg,-1, 1)
      tset n, 0,1,d,0, 0,0,0,1 '       +/- pi/k
      k = iif(sw and 1, 4, 2)
      cpi.ini k
   end if
end sub

function atan_t.f as long export
if t < 0 then return Inf
d = gencf(m, q, t, 1)
if sw then
   push n.x, d
   push n.y, cpi.f
   d = nextbfd(n)
end if
f = d
end function

'Eager evaluation for arcsin, arccos, inverse hyperbolic,
'Gudermannian and power functions.
' *****************************************************************************
'save and reset mpz stack positions
#define getoff b2 = off2, b8 = off8
#define setoff off2 = b2: off8 = b8

sub asin_t.f (byref r as cf, byref a as cf) export
dim as integer sw, getoff '            abs(a) <= 1
dim as long d
initprc(s, a)

sw = iif(fl, 2, 0) '                  cos or sin
if a.q(0) < 0 then sw or= 1 '         Cf_a < 0

n.ini: n.x = a: n.y = a
if fl then
   tset n, -1,0,0,1, 1,0,0,0 '        (1 - x^2) / x^2
else
   tset n, 1,0,0,0, -1,0,0,1 '         x^2 / (1 - x^2)
end if
dobf(r, n)

if r.q(0) > Inf and r.q(0) < 0 then
   print " illegal argument ";
   print iif(fl, "acos_t", "asin_t")
   catan.t = 0: cpi.t = 0
   empty(r): off8 = b8: exit sub
end if

csqrt.f r, r
least s, r

tset n, 0,-1,1,0, 0,0,0,1
catan.ini r
r.ini: cpi.ini 1
do
   d = catan.f
   if sw = 1 then '                   -atan
      push n.x, d
      d = nextlfd(n)
   elseif sw = 3 then '                pi - atan
      push n.x, d
      push n.y, cpi.f
      d = nextbfd(n)
   end if
   push r, d
loop until tm(r)
least r, s
setoff
end sub

constructor acos_t export
this.fl = -1
end constructor

'inverse hyperbolic functions
' *****************************************************************************
sub atanh_t.f (byref r as cf, byref a as cf) export
dim as integer getoff '                abs(a) < 1
initprc(s, a)

n.ini: n.x = a
tset n, 0,1,0,1, 0,-1,0,1 '           (x + 1) / (1 - x)
dolf(r, n)
least s, r

if leq0(r) then
   print " illegal argument atanh_t"
   clog.t = 0: off8 = b8
   empty(r): exit sub
end if

clog.ini r: r.ini
do
   push r, clog.f '                    log((x+1)/(1-x)) / 2
loop until tm(r)
least r, s
setoff
end sub

sub asinh_t.f (byref r as cf, byref a as cf) export
dim as integer getoff '                a >= 1 for acosh
dim as long d
initprc(s, a)

csqrt.t = 0: clog.t = 0
if fl and (a.q(0) < 1) then
   print " illegal argument acosh_t"
   empty(r): exit sub
end if

n.ini: n.x = a: n.y = a
d = iif(fl,-1, 1)
tset n, 1,0,0,d, 0,0,0,1 '             x^2 +/- 1
dobf(r, n)
least s, r

if r.q(0) < 0 then
   print " illegal argument acosh_t"
   empty(r): off8 = b8: exit sub
end if

tset n, 0,1,1,0, 0,0,0,1 '             x + sqrt(x^2 +/- 1)
csqrt.ini r: r.ini
do
   push n.y, csqrt.f
   push r, nextbfd(n)
loop until tm(r)
least s, r

clog.ini r: r.ini
do
   push r, clog.f
loop until tm(r)
least r, s
setoff
end sub

constructor acosh_t export
this.fl = -1
end constructor

'Gudermannian function
' *****************************************************************************
'2atan(exp(a)) - pi/2
sub gd_t.f (byref r as cf, byref a as cf) export
dim as integer getoff
initprc(s, a)

if is0(a) then
   cexp.t = 0: catan.t = 0
   cpi.t = 0: r = a: exit sub
end if

cexp.ini a: r.ini
do
   push r, cexp.f
loop until tm(r)
least s, r

cpi.ini 2: n.ini
tset n, 0,2,-1,0, 0,0,0,1
catan.ini r: r.ini
do
   push n.x, catan.f
   push n.y, cpi.f
   push r, nextbfd(n)
loop until tm(r)
least r, s
setoff
end sub

'inverse Gudermannian
'log(sec(a) + tan(a))
sub gdinv_t.f (byref r as cf, byref a as cf) export
dim as integer getoff '                abs(a) < pi/2
initprc(s, a)

r.ini: n.ini
tset n, 0,1,0,1, 0,0,1,0
csin.ini a: ccos.ini a
do
   push n.x, csin.f
   push n.y, ccos.f
   push r, nextbfd(n)
loop until tm(r)
least s, r

if leq0(r) then
   print " illegal argument gdinv_t"
   clog.t = 0: empty(r)
   setoff : exit sub
end if

clog.ini r: r.ini
do
   push r, clog.f
loop until tm(r)
least r, s
setoff
end sub

'power function
' *****************************************************************************
sub pow_t.f (byref r as cf, byref a as cf, byref e as cf) export
dim as integer getoff '                a > 0
initprc(s, a)
least s, e

if leq0(a) or (e.q(0) = Inf) then
   print " illegal argument pow_t"
   clog.t = 0: cexp.t = 0
   empty(r): exit sub
end if

n.ini: n.x = e
tset n, 1,0,0,0, 0,0,0,1 '             e * log(a)
r.ini: clog.ini a
do
   push n.y, clog.f
   push r, nextbfd(n)
loop until tm(r)
least s, r

cexp.ini r: r.ini
do
   push r, cexp.f
loop until tm(r)
least r, s
setoff
end sub

'(resume lazy evaluation)
'two adaptations for fast computing of Gosper's constant
' *****************************************************************************
'x / y:= x / 2y
sub halve (byref x as long, byref y as long)
if (x and 1) = 0 then
   x shr= 1
else
   y shl= 1
end if
end sub

'sin(3k)
sub sin3_t.ini (byval k as long) export
if abs(k) > 32767 then t =-1: exit sub
dim as long v = 1
   h.ini: s.ini: halve k, v
   tset h, 0,k,0,0, 0,v,0,1 '          x/y = k and tan(k/2) = z
   tset s, 0,1,1,0, 1,0,0,1 '          sin(k) = 2z / (z^2 + 1)

   m.ini: n.ini '                      tripling formula
   tset m,-4,0,0,3, 0,0,0,1 '          3 - 4*sin(k)^2
   tset n, 1,0,0,0, 0,0,0,1

   q.ini: k = abs(k): t = 0
   mpz_set_si q.a,-k
   mpz_mul_ui q.a, q.a, k '           -k^2
   mpz_set_ui q.b, v
end sub

function sin3_t.f as long export
if t < 0 then return Inf
do
   d = gencf(h, q, t, 0) '             tan(k/2)
   push s.x, d: push s.y, d
   d = nextbfd(s) '                    sin(k)
   push m.x, d: push m.y, d
   push n.x, d: push n.y, nextbfd(m)
   d = nextbfd(n) '                    sin(3k) = sin(k) * m
loop until d > Hold
f = d
end function

'tanh(square root(k))
sub tanhr_t.ini (byval k as short) export
if k < 1 then t =-1: exit sub
   m.ini: n.ini '                      1  k  k  k
   tset m, 0,1,0,0, 0,1,0,1 '          -- -- -- --
   tset n, 1,0,0,0, 0,0,0,1 '          1+ 3+ 5+ 7+...

   q.ini: t = 0
   mpz_set_ui q.a, k
   mpz_set_ui q.b, 1
   cquad.ini 1, 0, k
end sub

function tanhr_t.f as long export
if t < 0 then return Inf
   push n.x, gencf(m, q, t, 0) '       tanh(sqrt(k)) / sqrt(k)
   push n.y, cquad.f '                                        * sqrt(k)
f = nextbfd(n)
end function

'Elementary transcendental functions f(z * pi) for
'rational z = x / y, where x and y are small integers.
' *****************************************************************************
'adjust signs and check arguments
#macro isvalid
   if y < 0 then y = abs(y): x = -x
   t = iif(abs(x) > 32767,-1, 0)
   if y = 0 or y > 32767 then t =-1
   if t < 0 then exit sub
#endmacro

'evaluate exp(x * pi / y * 2)
function exppcf (byref n as cfa, byval x2 as long, byval y as long,_
 byref t as integer) as long
dim as integer fl = (getsw and 4) <> 0
dim as longint k
dim as long d

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   d = (t shl 1) + 1 '                   2z   z2+1 z2+4 z2+9
   k = clngint(y) * d '  z = 2x/y, 1 + ------ ---- ---- ----   exp(z*pi/2) A&S 4.2.42
   if abs(k) > Ulng then '             (1-z)+  3+   5+   7+...
      print "overflow exppcf, ";
      print "t =";t;", id =";n.id
      t = -1: exit do
   end if

   with n
      d = y * t
      mpz_set_ui .u(0), d
      mpz_set_ui .u(2), x2 '           tr(1) = (2x + ty*i)(2x - ty*i)
      mpz_addmul_ui .u(2),.u(0), d
      mpz_set_ui .u(0), k '            tr(0) = (2t + 1) * y
   end with

   if fl then prnt_ n, 0, t '          tr(2) = 1, tr(3) = 0
   do_genr n, 0
loop
exppcf = d
end function

sub exppi_t.ini (byval x as long, byval y as long) export
   isvalid
   halve y, x
   if abs(x) > 32767 then t =-1: exit sub
   n.ini: a2 = x * x: b = y
   tset n, 0,x + y,0,1, 0,y - x,0,1
end sub

function exppi_t.f as long export
if t < 0 then return Inf
f = exppcf(n, a2, b, t)
end function

'circular functions
' *****************************************************************************
'evaluate tan(x * pi / y * 4)
function tanpcf (byref n as cfa, byval x4 as long, byval y as long,_
 byref t as integer) as long
dim as integer i, fl = (getsw and 4) <> 0
dim as longint k(2), l(2), h
dim as long d
l(0) = Ulng
l(1) = Ulng shr 1: l(2) = l(1) '       max signed long

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   d = (t shl 1) + 1
   k(0) = clngint(y) * d

   h = clngint(y) * t '                z  1-z2 4-z2 9-z2
   k(1) = h + x4 '           z = 4x/y, -- ---- ---- ----   tan(z*pi/4) A&S 4.3.95
   k(2) = h - x4 '                     1+  3+   5+   7+...

   for i = 0 to 2
      if abs(k(i)) > l(i) then
         print "overflow tanpcf, ";
         print "t =";t;", id =";n.id
         t = -1: exit do
      end if
   next i

   with n
      mpz_set_ui .u(0), k(0)
      mpz_mul_ui .u(0),.u(0), y '      tr(0) = (2t + 1) * y^2
      mpz_set_si .u(2), k(1)
      mpz_mul_si .u(2),.u(2), k(2) '   tr(1) = (ty + 4x)(ty - 4x)
      mpz_set_ui .v(0), y
      mpz_mul_ui .v(0),.v(0), y '      tr(2) = y^2, tr(3) = 0
   end with

   if fl then prnt_ n,-1, t
   do_genr n,-1
loop
tanpcf = d
end function

sub tanpi_t.ini (byval x as long, byval y as long) export
   isvalid
   a = x: b = y: a mod= b
   if (a shl 1) > b then a -= b '      least residue mod 1
   halve b, a: halve b, a '            4x / y
   n.ini: tset n, 0,a,0,0, 0,b,0,b
end sub

function tanpi_t.f as long export
if t < 0 then return Inf
f = tanpcf(n, a, b, t)
end function

sub sinpi_t.ini (byval x as long, byval y as long) export
   isvalid
   a = x: b = y shl 1
   a mod= b: if a > y then a -= b '    least residue mod 2
   b = y: halve b, a '                 2x / y
   m.ini: n.ini
   tset m, 0,a,0,0, 0,b,0,b '          x = y = tan a / 2b
   if fl = 1 then '                   (x + y) / (xy + 1)
      tset n, 0,1,1,0, 1,0,0,1 '       A&S 4.3.24
   else '                             (1 - xy) / (1 + xy)
      tset n,-1,0,0,1, 1,0,0,1 '       A&S 4.3.25
   end if
end sub

function sinpi_t.f as long export
if t < 0 then return Inf
do
   d = tanpcf(m, a, b, t)
   push n.x, d: push n.y, d
   d = nextbfd(n)
loop until d > Hold
f = d
end function

constructor cospi_t export
this.fl = 2
end constructor

' *****************************************************************************
'ratio of Bessel functions of integer order I_v(x/y) / I_v-1(x/y)
function besscf (byref n as cfa, byval x2 as long, byval y as long,_
 byval v as integer, byref t as integer) as long
dim as integer fl = (getsw and 4) <> 0
dim as longint k
dim as long d

do
   if t then
      d = linfd(n)
      if d > Hold then exit do
   end if
   t += 1

   d = (v + t) shl 1
   k = clngint(y) * d
   if abs(k) < Ulng then '              x      x2       x2
      mpz_set_ui  n.u(0), k '          ---- -------- --------    I_v / I_v-1
      mpz_set_si  n.u(2), x2 '         2vy+ 2(v+1)y+ 2(v+2)y+...

   else
      print "overflow besscf, ";
      print "t =";t;", id =";n.id
      t = -1: exit do
   end if

   if fl then prnt_ n, 0, t
   do_genr n, 0
loop
besscf = d
end function

sub besseli_t.ini (byval v as integer,_
 byval x as long, byval y as long) export
   isvalid
   a2 = i * x * x: b = y
   n.ini: u = abs(v): t = 0
   tset n, 0,x,0,0, 0,2*u*y,0,1
end sub

function besseli_t.f as long export
if t < 0 then return Inf
f = besscf(n, a2, b, u, t)
end function

constructor besselj_t export
this.i = -1
end constructor

' *****************************************************************************
'root finding with Lagrange's 1767 Cf method
'put in Cf_c the first few partial quotients of the sought-for root
function polr_t.f (byref c as cf) as long export
dim as integer fl = (getsw and 2) <> 0
dim as integer i, j, k, sg
dim d as long, q as double
if t < 0 then return Inf
t += 1

sg = mpz_sgn(h(n))
if sg < 0 then
   for j = 0 to n
      mpz_neg h(j), h(j) '             make positive h.n
   next j
elseif sg = 0 then '                   rational root found,
   t = -1: return Inf '                coefficient h.n = 0
end if

if fl then '                           print all poly's
dim as clong e
dim s as zstring ptr
   if t = 1 then print "cpolroot"
   if (getsw and 16) <> 0 then
      for j = n to 0 step -1
         q = mpz_get_d_2exp(@e, h(j))
         print str(q);"b+";str(e);" ";
      next j
   else
      for j = n to 0 step -1
         s = allocate(mpz_sizeinbase(h(j), 10) + 2)
         mpz_get_str  s, 10, h(j)
         print *s;" ";
      next j
   end if
   print
end if

if c.r < c.w then '                    read Cf_c terms
   d = popd(c)
   if d <> 0 then '                    initial shifts
      mpz_set_si  tmp, d
      for i = 0 to n - 1
         for j = n - 1 to i step -1
            mpz_addmul h(j), h(j + 1), tmp
         next j
      next i
   end if
   goto recipr
end if

sg = 0: i = mpz_sgn(h(0))
for j = 1 to n
   k = i: i = mpz_sgn(h(j)) '          coefficient sign changes
   sg -= (k xor i) = -2
next j

q = 0
if sg = 1 then '                       reduced poly
   mpz_tdiv_q  tmp, h(n - 1), h(n)
   mpz_neg  tmp, tmp '                 sum of roots

   q = mpz_get_d(tmp)
   if q > 4 then

      if q > Maxd then
         i = mpz_sizeinbase(h(n - 1), 2)
         j = mpz_sizeinbase(h(n), 2)
         if abs(i - j) > mexp then '   declare infinite
            t = -1: return Inf
         end if

         q = Maxd '                    piecewise transmission
         mpz_set_si  tmp, Maxd
      end if

      for i = 0 to n - 1 '             prelim shift
         for j = n - 1 to i step -1
            mpz_addmul h(j), h(j + 1), tmp
         next j
      next i

   else
      q = 0
   end if
end if

for d = clng(q) to Maxd - 1 '          search
   mpz_set_si  tmp, 0
   for j = 0 to n
      mpz_add  tmp, tmp, h(j) '        evaluate h(1)
   next j

   sg = mpz_sgn(tmp) > 0
   if sg then exit for '               h(1) > 0

   for i = 0 to n - 1
      for j = n - 1 to i step -1 '     substitute x:= y + 1
         mpz_add h(j), h(j), h(j + 1)
      next j
   next i
next d

recipr:
i = (n - 1) shr 1
for j = 0 to i
   swap h(j), h(n - j) '               reciprocate
next j

t0 -= (d = 0)
sg = t0 > mexp and (t0 and 1) = 0
if sg then t = -1 '                    no real root or Cf truncated
if d > 0 and d < Maxd then t0 = 0 '    reset

f = d
end function

constructor polr_t export
redim h(dmx)
end constructor

