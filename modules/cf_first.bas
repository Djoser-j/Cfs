' *****************************************************************************
'Subject: Continued fraction arithmetic tryout
'Code   : FreeBasic 1.06.0 with GMP library 6.1.1

#include "cfr_lib.bi"
#inclib "cfr_math"

#define no_BR

'integer powers and cube root
' *****************************************************************************
#ifdef __BR
'Cf_a:= Cf_g ^ k, integer k
'breadth-first
sub Cpwr (byref a as cf, byref g as cf, byval k as short)
dim as integer fl = (getsw and 2) <> 0
dim n as cfa, b as cf
dim as short c, c2 = 1 shl 14
dim prc as cnv = cnvmax
if qirr(g) then prc = g.s
if @g <> @a then a = g

if k = 0 then
   cvcf a, 1: exit sub
elseif k < 0 then
   k = -k: Cinv a
end if
if k < 0 then
   print " illegal argument Cpwr"
   empty(a): exit sub
end if

n.ini
b = a
c = abs(k)
do until c and c2
   c2 shr= 1: loop '                   highest set bit(k)
c2 shr= 1
while c2 '                             L=>R square-and-multiply
   n.x = a: n.y = a
   tset n, 1,0,0,0, 0,0,0,1 '          a ^ 2
   dobf(a, n)
   least prc, a

   if fl then outcf "a^2", a

   if c and c2 then
      n.x = a: n.y = b
      tset n, 1,0,0,0, 0,0,0,1 '       a * b
      dobf(a, n)
      least prc, a

      if fl then outcf "a*b", a
   end if
   c2 shr= 1
wend
least a, prc
clrs
end sub

#else
'Cf_r:= Cf_a ^ k, integer k
'depth-first
sub Cpwr (byref r as cf, byref a as cf, byval k as short)
dim jx as const integer = 14
dim as integer b(jx), i, j
dim as cfa m(jx-1), n(jx-1)
dim as cf c = a
dim as long d, d0
dim prc as cnv = cnvmax
if qirr(a) then prc = a.s

if k = 0 then
   cvcf r, 1: exit sub
elseif k < 0 then
   k = -k: Cinv c
end if
if k < 0 then
   print " illegal argument Cpwr"
   empty(r): exit sub
end if

j =-1
do
   j += 1: b(j) = k and 1 '            bitstring for k
   k shr= 1
loop until j = jx or k = 0
if k > 0 then
   print " overflow Cpwr"
   empty(r): exit sub
end if
j -= 1

clrs
for i = 0 to jx - 1
   m(i).ini: n(i).ini
   tset m(i), 1,0,0,0, 0,0,0,1
   tset n(i), 1,0,0,0, 0,0,0,1
next i

r.ini: c.r = 0
do
   d = popd(c): d0 = d
   for i = j to 0 step -1 '            L=>R bitscan(k)
      push m(i).x, d
      push m(i).y, d
      d = nextbfd(m(i)) '              square
      if b(i) then
         push n(i).x, d
         push n(i).y, d0
         d = nextbfd(n(i)) '           multiply
      end if
   next i
   push r, d
loop until tm(r)
least r, prc
clrs
end sub

#endif

'Cf_r:= Cf_a ^ (1 / 3)
sub Cubrt (byref r as cf, byref a as cf)
dim as integer fl = (getsw and 2) <> 0
dim as integer sg, sw, t = 1
dim as cfa m, n, s
dim prc as cnv = cnvmax
if qirr(a) then prc = a.s
if @a <> @r then r = a

if is0(r) then exit sub
with r
   if .q(0) = Inf then exit sub
   sg = (.q(0) < 0)
   if sg then Cneg r
   sw = (.q(0) = 0)
   if sw then Cinv r
   '
   n.y = r '                           radicand R > 1
  .q(0) = int((.q(0) - (.q(1) = 1)) ^ (1 / 3))
  .q(1) = Inf '                        initial guess X
end with
s.ini: m.ini: n.ini

if fl then
   outcf "R", n.y
   outcf "X0", r
end if

while (t shr 1) < Cfx
   s.x = r: s.y = r '                  feed back
   m.x = r: r.ini
   tset s, 1,0,0,0, 0,0,0,1
   tset n, 0,0,1,0, 0,2,0,0
   tset m, 0,2,2,0, 0,0,0,3
   do '                                Newton iteration
      push n.x, nextbfd(s)
      push m.y, nextbfd(n) '           y:= R/ 2X^2
      push r, nextbfd(m) '             mean 2(X + y)/ 3
   loop until tm(r) or (r.w = t + 3)
   push r, Inf

   if fl then
      print "Cubrt ="; t
      outcf " X", r
   end if
   t shl= 1
wend
if sw then Cinv r
if sg then Cneg r
least r, prc
clrs
end sub

'continued fraction for e
function eul (byref t as integer) as long
if t = 0 then t = 1: return 2 '        first digit
dim as integer q = t \ 3
eul = iif(t - q * 3 = 2, t - q, 1)
t += 1
end function

' *****************************************************************************
'building functions
sub building
dim as cfa m, n, s, p, q
dim as cf a, b, c, e, f, g, h
dim as long x, y
dim as integer t
s.ini: n.ini: p.ini
q.ini: m.ini

print : ? "building" : ?

'' basic arithmetic
   ''linear
   '113x / 355
   tset n, 0,113,0,0, 0,0,0,355
   cvcf n.x, "3,7"
    outcf " x ", n.x
   dolf(c, n)
    outcf "x * 113/355 =", c

   '(7x - 22) / 7
   tset n, 0,7,0,-22, 0,0,0,7
   cvcf n.x, "3,7,16"
    outcf " x ", n.x
   dolf(c, n)
    outcf "x - 22/7 =", c

   ''bi-linear
   print
   tset n, 0,1,0,0, 0,0,1,0
   cvcf n.x, "0,3,7"
   cvcf n.y, 113/355
    outcf " x ", n.x
    outcf " y ", n.y
   dobf(c, n)
    outcf "x/y =", c

   tset n, 0,1,1,0, 0,0,0,1
   cvcf n.x, "-4,1,6,16"
   cvcf n.y, 22 / 7
    outcf " x ", n.x
    outcf " y ", n.y
   dobf(c, n)
    outcf "x+y =", c

   print
   'Kornerup & Matula p.1114
   tset s, 0,1,1,0, 0,0,0,1 '          sum
   cvcf n.x, 8/5
   cvcf n.y, 3/2
   tset n, 1,0,0,0, 0,0,0,1 '          product 1
   cvcf p.x, 7/9
   cvcf p.y, 1/5
   tset p, 1,0,0,0, 0,0,0,1 '          product 2
   c.ini
   do
      push s.x, nextbfd(n)
      push s.y, nextbfd(p)
      push c, nextbfd(s)
   loop until tm(c)
    outcf "(8/5*3/2)+(7/9*1/5)", c

   print
   tset s, 0,1,1,0,  0,0,0,1 '         sum
   tset n, 0,1,-1,0, 0,0,0,1 '         difference
   tset p, 1,0,0,0,  0,0,0,1 '         product 1
   tset q, 0,1,0,0,  0,0,1,0 '         quotient
   tset m, 1,0,0,0,  0,0,0,1 '         product 2
   t = 0
   a.ini: b.ini: c.ini
   e.ini: f.ini: g.ini: h.ini
   do
      x = 1 '                          fi
      y = eul(t) '                     e
      push f, x : push e, y
      push s.x, x: push s.y, y
      push n.x, x: push n.y, y
      push p.x, x: push p.y, y
      push q.x, x: push q.y, y
      x = nextbfd(s) '                 fi + e
      y = nextbfd(n) '                 fi - e
      push m.x, x: push m.y, y
      push a, x: push b, y
      push c, nextbfd(p)
      push g, nextbfd(q)
      push h, nextbfd(m)
   loop until tm(h)
    outcf "x", f: outcf "y", e
    outcf "x+y", a
    outcf "x-y", b
    outcf "x*y", c
    outcf "x/y", g
    outcf "(x+y)*(x-y)", h

   push f, Inf
   push e, Inf

   clrs
   print
   Cpwr c, f, 3
    outcf "y = fi^3 ", c
   Cubrt a, c
    outcf "x = y^1/3 ", a

   Cpwr c, c, 3
    outcf "y = fi^9 ", c
   Cubrt a, c: Cubrt a, a
    outcf "x = y^1/9 ", a

   Cpwr a, c, 3
    outcf "y = fi^27 ", a
   for t = 1 to 3: Cubrt a, a: next
    outcf "x = y^1/27 ", a

   print
   Cubrt c, e
    outcf "y = e^1/3 ", c
   Cpwr a, c, 3
    outcf "x = y^3 ", a

   Cubrt c, c
    outcf "y = e^1/9 ", c
   Cpwr a, c, 9
    outcf "x = y^9 ", a

   Cubrt c, c
    outcf "y = e^1/27 ", c
   Cpwr a, c, 27
    outcf "x = y^27 ", a

   print
   cvcf a, 2
   Cubrt c, a
    outcf "y = 2^ 1/3 ", c
   Cpwr a, c, 3
    outcf "x = y^3 ", a

   print
'' Newton stress-tests on rationals
   cvcf a, "0,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,2,1,1,1,1,2,1,2"
    outcf "x ", a

   Cpwr a, a, 27
    outcf "y = x^27", a
   for t = 1 to 3: Cubrt a, a: next
    outcf "x = y^1/27", a

   for t = 1 to 5: Cubrt a, a: next
    outcf "z = x^1/243", a
   Cpwr a, a, 243
    outcf "x = z^243", a

   print
   cvcf a, "-2,1,7,9,1,6171,1,15,2,1,3,10000"
    outcf "x ", a
   Cpwr a, a, 81
    outcf "y = x^81", a
   for t = 1 to 4: Cubrt a, a: next
    outcf "x = y^1/81", a

   for t = 1 to 5: Cubrt a, a: next
    outcf "z = x^1/243", a
   Cpwr a, a, 243
    outcf "x = z^243", a

''Pythagorean tuning
   print
   cvcf a, 3/2
   Cpwr n.x, a, 12 '                   12 fifths
   cvcf a, 2/1
   Cpwr n.y, a, 7 '                    7 octaves
   n.ini
   tset n, 0,1,0,0, 0,0,1,0
   dobf(c, n)
    outcf "ditonic comma", c

   clrs
end sub

' *****************************************************************************
'testing sub cf2q
sub cf2qtest
dim as cfa n
dim as cf a, b, c

print : ? "cf2qtest" : ?

   cvcf c, 3.14159265346744
   cf2q "~pi", c, Cfx, 0

   cvcf c, 1.5707963269607
   cf2q "~pi/2", c, Cfx, 0

   print
''near integer ratios in the solar system
   n.ini
''Jupiter / Saturn
   cvcf n.x, 29.45779 '                Saturn's orbit in years
   cvcf n.y, 11.86224 '                orbit of Jupiter
   tset n, 1,0,0,0, 0,1,-1,0
   dobf(c, n)
   cf2q "mean synodic period Jupiter/Saturn", c, 2, 0

   tset n, 0,1,0,0, 0,1,-1,0
   dobf(c, n)
   cf2q "triangle of Jupiter revolutions/conjunctions", c, 2, 0

''Venus / Earth
   cvcf n.x, 365.2564 '                Earth's orbit in days
   cvcf n.y, 224.7008 '                orbit of Venus
   tset n, 1,0,0,0, 0,1,-1,0
   dobf(b, n)
   cf2q "synodic period Venus/Sun in days", b, 4, 0

   tset n, 0,0,1,0, 0,1,-1,0
   dobf(c, n)
   cf2q "pentagram of years/inferior conjunctions", c, 3, 0

'' and Mars
   cvcf n.y, 686.9800 '                orbit of Mars
   tset n, 1,0,0,0, 0,-1,1,0
   dobf(c, n)
   cf2q "synodic period Sun/Mars", c, 4, 0

   n.x = b: n.y = c
   tset n, 0,0,1,0, 0,1,0,0
   dobf(c, n)
   cf2q "resonance Venus/Earth/Mars", c, 2, 0

   cvcf n.y, 365.2564
   tset n, 0,4,0,0, 0,0,1,0
   dobf(c, n)
   cf2q "years/resonance cycles", c, 3, 0

   print
''lunisolar calendar
   cvcf n.x, 365.2564 '                sidereal year
   cvcf n.y, 29.53059 '                synodic month
   tset n, 0,1,0,0, 0,0,1,0
   dobf(c, n)
   cf2q "Metonic cycle: lunar months/years", c, 4, 0

''periodicities related to the lunar orbit
   cvcf n.x, 27.21222
   tset n, 0,0,1,0, 0,1,0,0
   dobf(c, n)
   cf2q "Saros: draconic/ synodic months", c, 5, 0

   cvcf n.x, 27.55455
   tset n, 0,0,1,0, 0,1,0,0
   dobf(c, n)
   cf2q "anomalistic/ synodic months", c, 3, 0
   c.q(3) -= 1
   cf2q "Saros fit:", c, 3, 0

   cvcf n.x, 27.32166
   cvcf n.y, 27.21222
   tset n, 0,1,0,0, 0,0,1,0
   dobf(c, n)
   cf2q "draconic/ sidereal months", c, 1, 0
   c.q(1) -= 7
   cf2q "Saros fit:", c, 1, 0

   clrs
end sub

' *****************************************************************************
'testing sub cvcf
sub cvcftest
dim as cf a

print "cvcftest" : ?

   cvcf a, 343305231/121166534
    outcf "r2cf(x)", a

   cvcf a, 0.11000000011
    outcf "r2cf(x)", a
   cvcf a, 100000.10000011
    outcf "r2cf(x)", a
   cvcf a, 0.110001000000000000000001
    outcf "r2cf(x)", a

   print
   cvcf a, 38654705592
    outcf "r2cf(x)", a
   cvcf a, 38654705593
    outcf "r2cf(x+1)", a
   cvcf a, -3377777776.9077
    outcf "r2cf(x)", a
   cvcf a, -33777777769077, 10000
    outcf "q2cf(x)", a

   print
   cvcf a, 3377777777777700, 2600888888888929
    outcf "q2cf(x)", a
   cvcf a, "1.2987012987012487683198551285541"
    outcf "cvcf(x)", a
   cvcf a, 1.2987012987012487683198551285541
    outcf "r2cf(x)", a

   cvcf a, 0.0000000123
    outcf "r2cf(x)", a
   cvcf a, 1.0000000123
    outcf "r2cf(x)", a

   print
   'F45 / F46
   cvcf a, 1134903170, 1836311903
    rawcf "q2cf(x)", a, 0
   Cidem a
   print "cnv: ";a.s.t1;" * 2^";a.s.k

   'F91 / F92
   cvcf a, 4660046610375530309, 7540113804746346429
    rawcf "q2cf(x)", a, 0
   Cidem a
   print "cnv: ";a.s.t1;" * 2^";a.s.k

   cvcf a, "0.618033988749894848204586834365638117717304598812"
    outcf "cvcf(x)", a

   'F92 / F93 too big
   cvcf a, 7540113804746346429, 12200160415121876738
    outcf "q2cf(x)", a

   cvcf a, "3.141592e-10"
    outcf "cvcf(x)", a
   cvcf a, "3.-141592"
    outcf "cvcf(x)", a

   print
   cvcf a, "8,3,5,2,2,0,-6"
    outcf "readcf(x)", a
   cvcf a, "8,3,5,2,-4"
    outcf "readcf(x)", a
   cvcf a, "8,3,5,2,3,0,-7"
    outcf "readcf(x)", a
   Cidem a
    outcf "idem (x) ", a

   clrs
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

cvcftest
cf2qtest
building

print : ? "timer:"; csng(timer - tim); "s"
system
