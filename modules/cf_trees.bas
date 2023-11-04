' *****************************************************************************
'Subject: Expression trees for transcendental functions.
'Ref.   : P. J. Potts, Exact Real Arithmetic using Möbius Transformations,
'         PhD thesis Imperial College London, 1998, Chapter 10, pp.135-156.
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1


'Provide a large number of state registers in cfr_lib.bi (e.g. Cfx * 5)

#include "cfr_lib.bi"
#inclib "cfr_math"

const Trx = (Regs - 8) shr 2 '         expression tree height


'Potts' expression trees
' *****************************************************************************
type sqrt_l
declare function f (byval d as long) as long
   as cfa n(Trx)
   as integer t = 0
end type

type log_l
declare function f (byval d as long) as long
   as cfa n(Trx)
   as integer t = 0
end type

type exp_l
declare function f (byval d as long) as long
   as cfa m, n(Trx)
   as integer t = 0
end type

type tanh_l
declare function f (byval d as long) as long
   as cfa m, n(Trx)
   as integer fl = 1, t = 0
end type

type tan_l extends tanh_l
declare constructor ()
end type

type sinh_l
declare function f (byval d as long) as long
   as cfa m, n
   as tanh_l ctanh
   as integer i = 1, t = 0
end type

type sin_l extends sinh_l
declare constructor ()
end type

type atan_l
declare function f (byval d as long) as long
   as cfa h, m, n(Trx)
   as pi_t cpi
   as integer sw, t = 0
end type

type arctan_l
declare function f (byval d as long) as long
   as cfa h, m, n(Trx)
   as pi_t cpi
   as integer sw, t = 0
end type

type asin_l
declare function f (byval d as long) as long
   as cfa n
   as sqrt_l csqrt
   as atan_l catan
   as integer t = 0
end type

type atanh_l
declare function f (byval d as long) as long
   as cfa m, n(Trx)
   as integer t = 0
end type

type asinh_l
declare function f (byval d as long) as long
   as cfa m, n
   as sqrt_l csqrt
   as log_l clog
   as integer t = 0
end type

type modpi4_l
declare function f (byref q as long, byval d as long) as long
   as cfa n
   as pi_t cpi
   as integer sw, t = 0
end type


' *****************************************************************************
#macro mainloop(id)
   push n(Trx).x, d
   for i = Trx - 1 to 0 step -1
      push n(i).x, d
      push n(i).y, nextbfd(n(i+1))
   next i
   id = nextbfd(n(0))
#endmacro

function sqrt_l.f (byval d as long) as long export
dim i as integer
f = Inf

if t < 1 then
   if t < 0 then exit function

   for i = 0 to Trx
      n(i).ini
      tset n(i), 1,2,1,0, 0,1,2,1
   next i
   empty(n(Trx).y)
end if

mainloop(sqrt_l.f)
t += 1
end function

function log_l.f (byval d as long) as long export
dim as integer i, j
f = Inf

if t < 1 then
   if t < 0 then exit function

   n(0).ini
   tset n(0), 1,1,-1,-1, 0,1,1,0
   for i = 1 to Trx
      n(i).ini: j = (i shl 1) + 1
      tset n(i), i,j,i+1,0, 0,i+1,j,i
   next i
   empty(n(Trx).y)
end if

mainloop(log_l.f)
t += 1
end function

function exp_l.f (byval d as long) as long export
dim as integer i, j
f = Inf

if t < 1 then
   if t < 0 then exit function

   m.ini: tset m, 0,1,0,1, 0,-1,0,1 '  Soo{a} = (a+1)/(1-a)

   for i = 0 to Trx
      n(i).ini: j = (i shl 1) + 1
      tset n(i), j+1,j,j-1,j, j,j-1,j,j+1
   next i
   empty(n(Trx).y)
end if

push m.x, d
d = nextlfd(m)
mainloop(exp_l.f)
t += 1
end function

function tanh_l.f (byval d as long) as long export
dim as integer i, j
f = Inf

if t < 1 then
   if t < 0 then exit function

   m.ini: tset m, 0,1,0,1, 0,-1,0,1 '  Soo{a}

   n(0).ini
   if fl = -1 then
      tset n(0), 1,1,-1,-1, 2,0,0,2 '  tangent
      for i = 1 to Trx
         n(i).ini: j = (i shl 1) + 1
         tset n(i), j,j-2,j,j+2, j+2,j,j-2,j
      next i
   else
      tset n(0), 1,1,-1,-1, 0,2,2,0 '  hyperbolic tangent
      for i = 1 to Trx
         n(i).ini: j = (i shl 1) + 1
         tset n(i), j-2,j,j+2,j, j,j+2,j,j-2
      next i
   end if
   empty(n(Trx).y)
end if

push m.x, d
d = nextlfd(m)
mainloop(tanh_l.f)
t += 1
end function

constructor tan_l () export
this.fl = -1
end constructor

function sinh_l.f (byval d as long) as long export
f = Inf

if t < 1 then
   if t < 0 then exit function

   m.ini: n.ini
   tset m, 0,1,0,0, 0,0,0,2 '          x = y = tanh(a / 2)
   tset n, 0,1,1,0, -i,0,0,1 '         sinh:= (x + y) / (1 - xy)
   ctanh.t = 0: ctanh.fl = i
end if

push m.x, d
d = ctanh.f(nextlfd(m))
push n.x, d
push n.y, d
sinh_l.f = nextbfd(n)
t += 1
end function

constructor sin_l () export
this.i = -1
end constructor

function atan_l.f (byval d as long) as long export
dim as integer i, j
f = Inf

if t < 1 then
   if d > Hold then
      sw = (d < -1 or d > 0): t = 0
   else
      atan_l.f = d: t = -1
   end if
   if t < 0 then exit function

   i = 1
   if sw then
      h.ini: i = -1: cpi.ini 2
      tset h, 0,-1,1,0, 0,0,0,1
   end if
   m.ini: tset m, 0,1,0,1, 0,-i,0,i '  Soo{a} or Soo{1/a}

   n(0).ini
   tset n(0), 1,1,-1,-1, 2,0,0,2
   for i = 1 to Trx
      n(i).ini: j = (i shl 1) + 1
      tset n(i), j,i,0,i+1, i+1,0,i,j
   next i
   empty(n(Trx).y)
end if

push m.x, d
d = nextlfd(m)
mainloop(d)
if sw then '                           pi/2 - arccot
   push h.x, d
   push h.y, cpi.f
   d = nextbfd(h)
end if
atan_l.f = d
t += 1
end function

function arctan_l.f (byval d as long) as long export
dim as integer i, j
dim k as long
f = Inf

if t < 1 then
   if d > Hold then
      sw = (d < -1 or d > 0): t = 0
   else
      arctan_l.f = d: t = -1
   end if
   if t < 0 then exit function

   if sw then
      m.ini: h.ini: cpi.ini 2
      tset m, 0,0,0,1, 0,1,0,0 '       1 / a
      tset h, 0,-1,1,0, 0,0,0,1 '      pi/2 - arccot
   end if

   for i = 0 to Trx
      n(i).ini
      k = i + 1: j = (i shl 1) + 1
      tset n(i), 0,1,0,0, k*k,0,0,j
   next i
   empty(n(Trx).y)
end if

if sw then
   push m.x, d
   d = nextlfd(m)
   mainloop(d)
   push h.x, d
   push h.y, cpi.f
   d = nextbfd(h)
else
   mainloop(d)
end if
arctan_l.f = d
t += 1
end function

function asin_l.f (byval d as long) as long export
f = Inf

if t < 1 then
   if d > Hold then
      t = (d < -1 or d > 0)
   else
      asin_l.f = d: t = -1
   end if
   if t < 0 then exit function

   n.ini: tset n, 1,0,0,0, -1,0,0,1 '  atan(sqrt(a^2 / (1-a^2)))
   csqrt.t = 0: catan.t = 0
end if

push n.x, d
push n.y, d
d = csqrt.f(nextbfd(n))
asin_l.f = catan.f(d)
t += 1
end function

function atanh_l.f (byval d as long) as long export
dim as integer i, j
f = Inf

if t < 1 then
   if t < 0 then exit function

   m.ini: tset m, 0,1,0,1, 0,-1,0,1 '  Soo{a}

   n(0).ini
   tset n(0), 1,1,-1,-1, 0,2,2,0
   for i = 1 to Trx
      n(i).ini: j = (i shl 1) + 1
      tset n(i), i,j,i+1,0, 0,i+1,j,i
   next i
   empty(n(Trx).y)
end if

push m.x, d
d = nextlfd(m)
mainloop(atanh_l.f)
t += 1
end function

function asinh_l.f (byval d as long) as long export
f = Inf

if t < 1 then
   if t < 0 then exit function

   m.ini: n.ini
   tset m, 1,0,0,1, 0,0,0,1 '          a^2 + 1
   tset n, 0,1,1,0, 0,0,0,1 '          log(a + sqrt(1+a^2))
   csqrt.t = 0: clog.t = 0
end if

push m.x, d
push m.y, d
push n.x, d
push n.y, csqrt.f(nextbfd(m))
d = nextbfd(n)
asinh_l.f = clog.f(d)
t += 1
end function

function modpi4_l.f (byref q as long, byval d as long) as long export
dim as integer wx, wy
dim h as long
f = Inf

if t < 1 then
   if t < 0 then exit function

   n.ini: cpi.ini 4
   tset n, 0,1,0,0, 0,0,1,0
   sw = -1
end if

h = cpi.f
push n.x, d
push n.y, h
if sw then
   q = nextbfd(n) '                    q:= floor[x / pi/4]
   t += 1
   if q <= Hold then return q

   wx = n.x.w: wy = n.y.w
   tset n, 0,1,-q,0, 0,0,0,1 '         x - q * pi/4
   n.x.w = wx: n.y.w = wy: sw = 0
end if

modpi4_l.f = nextbfd(n)
t += 1
end function

'module
' *****************************************************************************
sub treetest ()
dim as cfa n
dim as cf a, b, c
dim as poly g
dim cquad as quad_t
dim sqrt as sqrt_l, cpi as pi_t
dim clog as log_l, cexp as exp_l
dim ctanh as tanh_l, ctan as tan_l
dim csinh as sinh_l, csin as sin_l
dim catanh as atanh_l, catan as atan_l
dim carctan as arctan_l
dim casinh as asinh_l, casin as asin_l
dim cmodpi4 as modpi4_l
dim as long d, u, v, w, t
dim i as integer
dim prc as cnv

#macro cfun(rd, vd, pd)
   clrs
   print
   b.ini: c.ini
   cquad.ini u, v, w
   rd.t = 0: vd.t = 0
   do
      d = rd.f(cquad.f)
      push b, d
      push c, vd.f(d)
   loop until tm(b) and tm(c)
   least c, pd.s
#endmacro

for i = 0 to 5
print "i ="; i : ?
   select case i
   case 0: u = 1: v =-1: w = 1
   case 1: u =-2: v = 8: w = 7
   case 2: u =-25: v = 52: w = 26
   case 3: u = 9424900: v = 0: w = 5396329
   case 4: u = 29: v = 0: w = 8
   case else: u = 64: v = 0: w = 1
   end select

   prc = cnvmax
   a.ini: cquad.ini u, v, w
   do : push a, cquad.f
   loop until tm(a)
    outcf "cf(z)", a
   least prc, a

   print
   b.ini: c.ini
   cquad.ini u, v, w
   n.ini: sqrt.t = 0
   tset n, 1,0,0,0, 0,0,0,1
   do
      d = sqrt.f(cquad.f)
      push b, d
      push n.x, d
      push n.y, d
      push c, nextbfd(n)
   loop until tm(b) and tm(c)
   least b, prc
    outcf "sqrt(z)", b
   rate sqrt.t
   least c, b
    outcf "sqrt(z)^2", c

   cfun(cexp, clog, b)
    outcf "exp(z)", b
   rate cexp.t
    outcf "log(exp(z))", c

   cfun (csin, casin, b)
    outcf "sin(z)", b
   rate csin.t
    outcf "asin(sin(z))", c

   clrs
   print
   b.ini: c.ini
   cquad.ini u, v, w
   cmodpi4.t = 0: ctan.t = 0
   do : d = cmodpi4.f(t, cquad.f)
   loop while d = Hold
   n.ini
   select case t and 3
   case 0: tset n, 0,1,0, 0, 0, 0,0,1 ' S+
   case 1: tset n, 0,1,0, 1, 0,-1,0,1 ' Soo
   case 2: tset n, 0,0,0,-1, 0, 1,0,0 ' S-
   case 3: tset n, 0,1,0,-1, 0, 1,0,1 ' S0
   end select
   do
      push b, d
      push n.x, ctan.f(d)
      push c, nextlfd(n)
      d = cmodpi4.f(t, cquad.f)
   loop until tm(b) and tm(c)
    print "sw "; t and 3
   least b, prc
    outcf "modpi4(z)", b
   least c, b
    outcf "tan(modpi4(z))", c
   rate ctan.t

   cfun (ctan, catan, b)
    outcf "tan(z)", b
   rate ctan.t
    outcf "atan(tan(z))", c

   cfun (ctan, carctan, b)
    outcf "tan(z)", b
   rate ctan.t
    outcf "arctan(tan(z))", c

   clrs
   'compare Potts' 1996 tensor sequence...
   print
   c.ini: carctan.t = 0
   cquad.ini u, v, w
   do
      push c, carctan.f(cquad.f)
   loop until tm(c)
   least c, prc
    outcf "arctan(z)", c
   rate carctan.t

   clrs
   'and 1998 expression tree methods
   print
   c.ini: catan.t = 0
   cquad.ini u, v, w
   do
      push c, catan.f(cquad.f)
   loop until tm(c)
   least c, prc
    outcf "atan(z)", c
   rate catan.t

   cfun (csinh, casinh, b)
    outcf "sinh(z)", b
   rate csinh.t
    outcf "asinh(sinh(z))", c

   cfun (ctanh, catanh, b)
    outcf "tanh(z)", b
   rate ctanh.t
    outcf "atanh(tanh(z))", c

   clrs
   print : ? "_______________________________"
next i
end sub

sub trantest ()
dim as cfa n
dim as cf a, c
dim clog as log_t, cexp as exp_t
dim ctan as tan_t, ctanh as tanh_t
dim csin as sin_t, csinh as sinh_t
dim catan as atan_t, catanh as atanh_t
dim casin as asin_t, casinh as asinh_t
dim cquad as quad_t, csqrt as sqrt_t
dim i as integer, prc as cnv
dim as long u, v, w

#macro cfun(id, ag, msg)
   prc = ag.s
   c.ini: id.ini ag
   do
      push c, id.f
   loop until tm(c)
   least c, prc
    outcf msg, c
#endmacro

for i = 0 to 5
   print "i ="; i : ?

   select case i
   case 0: u = 1: v =-1: w = 1
   case 1: u =-2: v = 8: w = 7
   case 2: u =-25: v = 52: w = 26
   case 3: u = 9424900: v = 0: w = 5396329
   case 4: u = 29: v = 0: w = 8
   case else: u = 64: v = 0: w = 1
   end select

   a.ini: cquad.ini u,v,w
   do
      push a, cquad.f
   loop until tm(a)
    outcf "cf(z)", a
   least a, cnvmax

   print
   cfun(csqrt, a, "sqrt(z)")
   rate csqrt.t: prc = c.s
   n.ini: n.x = c: n.y = c
   tset n, 1,0,0,0, 0,0,0,1
   dobf(c, n)
   least c, prc
    outcf "sqrt(z)^2", c

   print
   cfun(cexp, a, "exp(z)")
   rate cexp.t
   cfun(clog, c, "log(exp(z))")

   print
   cfun(csin, a, "sin(z)")
   rate csin.t
   casin.f c, c
    outcf "asin(sin(z))", c

   ? : ? "sw  "
   ? "modpi4(z)"
   ? : ? "tan(modpi4(z))"
   ? : ? "Cfx"
   ? : ? : ? : ? : ? : ?

   print
   cfun(ctan, a, "tan(z)")
   rate ctan.t
   cfun(catan, c, "arctan(tan(z))")

   print
   cfun(catan, a, "arctan(z)")
   rate catan.t

   ? : ? : ? : ?

   print
   cfun(csinh, a, "sinh(z)")
   rate csinh.t
   casinh.f c, c
    outcf "asinh(sinh(z))", c

   print
   cfun(ctanh, a, "tanh(z)")
   rate ctanh.t
   catanh.f c, c
    outcf "atanh(tanh(z))", c

   clrs
   print : ? "_______________________________"
next i
end sub

'logistic map x <- rx * (1 - x)
'chaotic orbit r = 4 (bit-shift map)
sub logistic (byval ix as const short)
dim as integer i, t
dim as cfa m, n(ix)
dim as cf a, c
dim as long d

print : ? "logistic map" : ?

cvcf a, 43/64
outcf "a", a

setvw 1 + Cfx \ 2
m.ini: c = a
for i = 0 to ix
   m.x = c: m.y = c
   tset m, -4,4,0,0, 0,0,0,1
   dobf(c, m)
next i
outcf "i_"+str(i), c

clrs
setvw Cfx
if ix < 1 then exit sub

if Regs < ix+2 then
   print : ? "recursion demo needs more state registers:"
   print "set Regs >"; ix+1; " in cfr_lib.bi and recompile.": ?

else
   for i = 0 to ix
      n(i).ini
      tset n(i), -4,4,0,0, 0,0,0,1
   next i

   print "recursive: true orbit"
   n(0).x = a
   n(0).y = a
   c.ini: t = 1
   do
      i = 0
      do
         i += 1
         d = nextbfd(n(i-1))
         push n(i).x, d
         push n(i).y, d
         if d > Hold then '            current depth
            if i = t then t = i + 1
         end if
      loop while (t > i) and (i < ix)
      push c, nextbfd(n(i))
   loop until tm(c)
   outcf "t_"+str(t), c

end if

clrs
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

if Regs < Cfx shl 2 then
   print "To produce valid results this program needs more state registers:"
   print "set Regs ~ Cfx * 5 (or define trees) in cfr_lib.bi and recompile."
   print : ? "No trees, applying defaults..." : ?

   trantest
else
   treetest
end if

'demo applies to Cfx = 101
logistic(142)
'50 correct partial quotients left
logistic(283)
'with 284 iterations the true orbit is lost

print : ? "timer:"; csng(timer - tim); "s"
system
