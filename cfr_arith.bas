' *****************************************************************************
'Subject: R. W. Gosper's continued fraction arithmetic.
'Refs.  : HAKMEM, A.I. Lab Memo #239, M.I.T., 1972
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#define __LIB

#include once "cfr_lib.bi"


' *****************************************************************************
'global variables
dim shared as integer Verb = 0 '       verbose mode
dim shared as integer Vwx = Cfx - 1 '  view-port < Cfx

dim shared as integer trlf '           truncate linear form registers
dim shared as integer trlon = -1 '     do not truncate in outcf()
trlf = 4.5 * Cfx / mp_bits_per_limb

'global mpz pointers
dim shared as mpz_ptr rat, reg
'initialize all mpz integers with space for h bits
dim as ulong h = Cfx * 4
dim as integer i, j

redim pw2(-1 to mant)
pw2(-1) = 0.5 '                        initialize
for i = 0 to mant
   pw2(i) = pw2(i - 1) * 2
next i

rat = allocate(Rats * 2 * len(mpz_t))
tmp = @rat[0]
for i = 1 to Rats
   for j = 1 to 2
      mpz_init2(tmp, h) '              mpz pointer pairs
      tmp += 1
   next j
next i

reg = allocate(Regs * 8 * len(mpz_t))
tmp = @reg[0]
for i = 1 to Regs
   for j = 1 to 8
      mpz_init2(tmp, h) '              mpz register list
      tmp += 1
   next j
next i

tmp = allocate(len(mpz_t))
mpz_init2(tmp, h)
off2 = 0: off8 = 0

'constructor methods
' *****************************************************************************
constructor cf export
redim q(Cfx)
end constructor

constructor cfa export '               default
redim x.q(Cfx), y.q(Cfx)
end constructor

constructor cfa (byval fl as integer) export
'base register set without Cfs x and y
redim x.q(0), y.q(0)
end constructor

constructor poly export
redim a(dmx)
end constructor

' *****************************************************************************
'switch to verbose mode, set
'bit 1 (1) full Cf, ignoring precision
'bit 2 (2) print iteration steps
'bit 3 (4) linear form throughput streams
'bit 4 (8) bi-linear throughput streams
'bit 5 (16) print radix-2 exponential format
sub setsw (byval sw as integer) export
   Verb = sw
end sub
'return the current mode
function getsw as integer export
getsw = Verb
end function

'set view-port
sub setvw (byval sw as integer) export
   Vwx = sw
end sub
'return the current size
function getvw as integer export
getvw = Vwx
end function

'clear mpz stack; use often to reclaim register space
sub clrs export
   off2 = 0: off8 = 0
end sub

' *****************************************************************************
'assign mpz variable pair
sub pair.ini export
dim as mpz_ptr t = @rat[off2 * 2]
if off2 + 1 > Rats then
   print "out of mpz pairs:";
   print off2 + 1 : end
end if

a = t: b = t + 1
id = off2
off2 += 1
end sub

'assign x and y to pair
sub tset overload (byref p as pair, byval x as long, byval y as long) export
   mpz_set_si p.a, x
   mpz_set_ui p.b, abs(y)
end sub

'reset cf indices, initialize recurrence
sub cf.ini export
   r = 0: w = 0
   s.t0 = 1: s.t1 = 0: s.k = 0
end sub

'assign mpz registers to Cfa structure
sub cfa.ini export
dim as mpz_ptr t = @reg[off8 * 8]
dim as integer i
if off8 + 1 > Regs then
   print "out of mpz registers:";
   print off8 + 1 : end
end if

for i = 0 to 3
   u(i) = t: t += 1
   v(i) = t: t += 1
next i
id = off8
off8 += 1
end sub

'assign (a b c d) / (e f g h) to tensor and initialize Cf x and y
sub tset (byref n as cfa,_
 byval a as long, byval b as long, byval c as long, byval d as long,_
 byval e as long, byval f as long, byval g as long, byval h as long) export
with n
   mpz_set_si .u(0), a: mpz_set_si .v(0), e
   mpz_set_si .u(1), b: mpz_set_si .v(1), f
   mpz_set_si .u(2), c: mpz_set_si .v(2), g
   mpz_set_si .u(3), d: mpz_set_si .v(3), h
  .x.ini : .y.ini
end with
end sub

'assign mpz coefficients to polynomial
sub polr_t.ini (byref g as poly) export
dim as mpz_ptr tp = @reg[off8 * 8]
dim as integer i, k
n = g.n
if n < 1 or n > dmx then
   print "polr_t: degree out of bounds"
   t = -1: exit sub
end if
k = 1 + n \ 8
if off8 + k > Regs then
   print "polr_t: out of mpz registers:";
   print off8 + k : end
end if

for i = 0 to n
   this.h(i) = tp: tp += 1
   mpz_set_si this.h(i), g.a(i) '     (coefficients will explode)
next i
off8 += k
t = 0: t0 = 0
end sub

'mpz-matrix operations
' *****************************************************************************
'pre-multiply trafo * matrix, columns 0 & 2 or 1 & 3
sub do_absorb (byval d as long, byref n as cfa, byval i as integer)
dim as integer j = i + 2
with n
   if d > Inf then '                   absorb
      '.u(j) += d * .u(i)
      mpz_set_si  tmp, d
      mpz_addmul .u(j), .u(i), tmp
      'exchange
      swap .u(i), .u(j)
      'repeat with .v
      mpz_addmul .v(j), .v(i), tmp
      swap .v(i), .v(j)

   else
      mpz_set .u(j), .u(i) '           done: copy
      mpz_set .v(j), .v(i)
   end if
end with
end sub

sub do_emit (byval d as long, byref n as cfa, byval i as integer)
dim as integer j = i + 2
with n
   '.u(i) -= d * .u(j)
   mpz_set_si  tmp, d
   mpz_submul .u(i), .u(j), tmp '      emit
   'exchange
   swap .u(j), .u(i)
   'repeat with .v
   mpz_submul .v(i), .v(j), tmp
   swap .v(j), .v(i)
end with
end sub

sub do_genr (byref n as cfa, byval sw as integer) export
with n '                               trafo |u(0) u(2)|*|u(1) v(1)|
   '.u(3) *= .u(2)                           |v(0)  0  | |u(3) v(3)|
   mpz_mul .u(3), .u(3), .u(2)
   '.u(3) += .u(0) * .u(1)             absorb from sqrtcf, gencf,
   mpz_addmul .u(3), .u(1), .u(0) '        acotcf, besscf, exppcf
   'repeat with .v()
   mpz_mul .v(3), .v(3), .u(2)
   mpz_addmul .v(3), .v(1), .u(0)
   if sw then '.u(1) *= .v(0)          tanpcf
      mpz_mul .u(1), .u(1), .v(0)
      mpz_mul .v(1), .v(1), .v(0)
   end if
   'exchange
   swap .u(1), .u(3)
   swap .v(1), .v(3)
end with
end sub

'transpose matrices
sub trans (byref n as cfa)
   swap n.u(2), n.u(1)
   swap n.v(2), n.v(1)
end sub

'turn coefficient cube
sub turn (byref n as cfa)
   swap n.u(2), n.v(0)
   swap n.u(3), n.v(1)
end sub

' *****************************************************************************
'trial quotient registers n.u(i) / n.v(i)
function zquot (byref n as cfa, byval i as integer) as long
dim as double q, v
dim as clong r, s
zquot = Inf

v = mpz_get_d_2exp(@s, n.v(i))
if v then
   q = mpz_get_d_2exp(@r, n.u(i)) / v

   r -= s
   if abs(r) > mexp then
      ''print "overflow in zquot, 2^"; r
      ''print "id";n.id;", i";i
      exit function
   end if
   select case r
   case is < -1
      q = int(q * 0.25) '             -1 or 0
   case is < cfb
      q = int(q * pw2(r))
   case else
      q = int(q * pw2(cfb)) '          too big
   end select

   if abs(q) > Maxd then '             piecewise transmission
      q = iif(q > 0, Maxd, 1 - Maxd)
   end if

zquot = clng(q)
end if
end function

'next convergent
sub nxtcnv (byref q as cnv, byval d as long)
with q
   .s = .t0 + .t1 * d
   .t0 = .t1: .t1 = .s
   while abs(.t1) >= 1 '               normalize
      .t1 *= 0.5
      .t0 *= 0.5
      .k += 1
   wend
end with
end sub

'push digit in Cf queue
sub push (byref c as cf, byval d as long) export
dim as integer i, v = c.w

with c
   if v = Cfx then '                   shift down fifo queue,
      if .q(Cfx) <> Inf and .r > 1 then ' transient Cfs only
        .r -= 1: v -= .r
         for i = 1 to v
           .q(i) = .q(i + .r)
         next i
        .r = 1: .w = v
      end if
   end if

   if v < Cfx then
      if d > Inf then
         if .r = 0 then nxtcnv .s, d ' update convergent size
        .w += 1 '                      advance write index
      end if
     .q(.w) = Hold
     .q(v) = d

   else
     .q(Cfx) = Inf
   end if
end with
end sub

'pop Cf digit
function popd (byref c as cf) as long export
dim as long d = c.q(c.r)
if d > Inf then c.r += 1 '             advance read index
popd = d
end function

' *****************************************************************************
'print digit <-in
sub strin (byref s as string, byval d as long)
select case d
   case Hold : s += "Hold"
   case Inf  : s += "Inf"
   case Maxd : s += "Maxd"
   case else : s += str(d)
end select
print s;
end sub

'print bilinear fractional form
sub outfrm (byref n as cfa)
dim as integer i, fl = (Verb and 16) <> 0
dim as zstring ptr s
dim as string fs, gs = ""
dim as double h
dim as clong r

fs = "["+str(n.id)+"] "
if fl then
   for i = 0 to 3
      h = mpz_get_d_2exp(@r, n.u(i))
      fs += " "+str(h)+"b+"+str(r)+" "
      h = mpz_get_d_2exp(@r, n.v(i))
      gs += "  "+str(h)+"b+"+str(r)
   next i
else
   for i = 0 to 3
      s = allocate(mpz_sizeinbase(n.u(i), 10) + 2)
      mpz_get_str  s, 10, n.u(i)
      fs += " "+*s+" "
      s = allocate(mpz_sizeinbase(n.v(i), 10) + 2)
      mpz_get_str  s, 10, n.v(i)
      gs += "  "+*s
   next i
end if
print fs;"/";gs;"  "
end sub

'bi-homographic routines
' *****************************************************************************
'update tensor for input x (sw = 0), or y (sw = 1)
function absorb overload (byref n as cfa, byval sw as integer) as integer
dim as integer fl = (Verb and 8) <> 0
dim as long d
absorb = 0

   if sw then
      d = popd(n.y) '                  Cf digit
      if fl then strin " <-y  ", d : ?
   else
      d = popd(n.x)
      if fl then strin " <-x  ", d : ?
   end if
   if d = Hold then exit function '    wait for next digit

   if sw then trans n '                transpose for y
   do_absorb d, n, 0 '                 trafo * tensor
   do_absorb d, n, 1
   if sw then trans n '                restore
absorb = -1
end function

'update tensor for output digit d
sub emit01 (byref n as cfa, byval d as long)
if d > Inf then
   turn n '                            tilted cube
   do_emit d, n, 0 '                   trafo * tensor
   do_emit d, n, 1
   turn n '                            restore
end if
end sub

'next bi-linear form digit z(x, y)
function nextbfd (byref n as cfa) as long export
dim as integer i, sw, fl = (Verb and 8) <> 0
dim as long d(3), df(2)
nextbfd = Hold
do
   sw = -1 '                           flat z-range
   for i = 0 to 3
      d(i) = zquot(n, i) '             four corner values z
      if i then sw and= d(i - 1) = d(i)
   next i

   if fl then
      outfrm n
      for i = 0 to 3 : strin "  ", d(i) : next
      if sw then print " ->out"
   end if

   if sw then '                        happens eventually
      emit01 n, d(0)
      nextbfd = d(0): exit do
   end if

   for i = 1 to 2 '                    choose Cf x or y
      df(i) = abs(d(0) - d(i)) '      (this test limits Maxd:
      df(0) = abs(d(3 - i) - d(3)) '   Maxd - Inf < 2^31)
      if df(i) < df(0) then df(i) = df(0)
   next i
   sw = iif(df(1) > df(2), 1, 0)

loop while absorb(n, sw) '             digit stream <-in
end function

'homographic routines separate for efficiency
' *****************************************************************************
'update matrix for input x
function absorb (byref n as cfa) as integer
dim as integer fl = (Verb and 4) <> 0
dim as long d
absorb = 0
   d = popd(n.x) '                     Cf_n.x digit
   if fl then strin " <-x  ", d : ?
   if d = Hold then exit function '    wait for next digit
   do_absorb d, n, 1 '                 trafo * matrix columns 1 & 3
absorb = -1
end function

'update matrix for output digit d
sub emit_1 (byref n as cfa, byval d as long)
if d > Inf then
   swap n.u(3), n.v(1) '               turn
   do_emit d, n, 1 '                   trafo * matrix columns 1 & 3
   swap n.u(3), n.v(1)
end if
end sub

'truncate (b d) / (f h) registers
sub truncate (byref n as cfa)
dim as integer i, fl = (Verb and 4) <> 0
dim as long d = mpz_size(n.v(3))

if d > trlf then
   d -= trlf
   d *= mp_bits_per_limb
   for i = 1 to 3 step 2
      mpz_tdiv_q_2exp  n.u(i), n.u(i), d
      mpz_tdiv_q_2exp  n.v(i), n.v(i), d
   next i
   if fl then print "red M (";str(d);")"
end if
end sub

'linear form digit z(x)
function linfd (byref n as cfa) as long export
dim as integer i, sw, fl = (Verb and 4) <> 0
dim as long d(1)
linfd = Hold

   for i = 1 to 3 step 2
      d(i shr 1) = zquot(n, i) '       both z values
   next i
   sw = (d(0) = d(1))

   if fl then
      outfrm n
      for i = 0 to 1 : strin "  ", d(i) : next
      if sw then print " ->out"
   end if

   if sw then '                        flat z-range
      emit_1 n, d(0)
      if trlon then truncate n '       limit register size
      linfd = d(0)
   end if
end function

'the above wrapped in...
'next linear form digit z(x)
function nextlfd (byref n as cfa) as long export
dim as long d
nextlfd = Hold
do
   d = linfd(n) '                      digit z(x)
   if d > Hold then exit do
loop while absorb(n) '                 Cf_n.x stream <-in
nextlfd = d
end function

'precision routines
' *****************************************************************************
'compare convergents
function cmpcnv (byref a as cnv, byref b as const cnv) as integer
dim as integer dk = a.k - b.k, sw = 0
dim as double df

if dk <> 0 then
   sw = iif(dk < 0,-1, 1)
else
   df = abs(a.t1 / b.t1)
   sw += df < 1 '  -1
   sw -= df > 1 '   1
end if
cmpcnv = sw
end function

'return least precision in prc
sub least overload (byref prc as cnv, byref c as cf) export
if c.w > (Cfx - 2) then
   if cmpcnv(prc, c.s) = 1 then prc = c.s
end if
end sub

'set a.s to least precision
sub least (byref a as cf, byref prc as const cnv) export
if qirr(a) then
   if cmpcnv(a.s, prc) = 1 then a.s = prc
else
   a.s = cnvmax
end if
end sub

sub least (byref a as cf, byref c as cf) export
least a, c.s
end sub

' *****************************************************************************
'compare Cfs, mismatch index r,
'return a < b: -1, a = b: 0, a > b: 1
function comp (byref r as integer,_
 byref a as cf, byref b as cf) as integer export
dim as integer i, s
   for i = 0 to Cfx - 1 '              first discrepant terms
      if a.q(i) - b.q(i) then exit for
      if a.q(i) <= Inf then exit for ' both finite
   next i
   r = a.q(i): if r <= Inf then r = Maxd
   s = b.q(i): if s <= Inf then s = Maxd
   if i and 1 then swap r, s '         invert decision
comp = sgn(r - s): r = i + 1
end function

'negation and reciprocation
'Cf_a:= -Cf_a
sub Cneg (byref a as cf) export
dim s as cnv = a.s
dim as cfa n
   n.ini: n.x = a
   tset n, 0,-1,0,0, 0,0,0,1
   dolf(a, n)
a.s = s
off8 -= 1
end sub

'Cf_a:= 1 / Cf_a
sub Cinv (byref a as cf) export
dim as integer i, sw
with a
   sw = (.q(0) > Inf) and (.q(0) < 0)
   if sw then Cneg a

   if .q(0) = 0 then
     .w -= 1
      for i = 0 to .w '                remove initial 0
        .q(i) = .q(i + 1)
      next i

   else
      if .w < Cfx then .w += 1
      for i = .w to 1 step -1
        .q(i) = .q(i - 1)
      next i
     .q(i) = 0 '                       prefix 0
     .q(.w) = Inf
   end if

   if sw then Cneg a
end with
end sub

'identity
sub Cidem (byref a as cf) export
dim s as cnv = a.s
dim as cfa n
   n.ini: n.x = a
   tset n, 0,1,0,0, 0,0,0,1
   dolf(a, n)
a.s = s
off8 -= 1
end sub

'Cf output
' *****************************************************************************
'check convergent size
function bigcnv (byref q as cnv) as integer
dim as double t
bigcnv = -1
if q.k <= mant then
   t = q.t1 * pw2(q.k)
bigcnv = (t * t > Flint)
end if
end function

'print long float convergent p.t1 / q.t1
sub prntcnv (byref p as cnv, byref q as cnv, byval sw as integer)
dim as integer dk = p.k - q.k
dim as double df

if sw then
   if abs(dk) <= mant then '           fp quotient
      df = p.t1 / q.t1
      if dk >= 0 then
         print df * pw2(dk);
      else
         print df / pw2(abs(dk));
      end if
   end if
else
   if p.k < mant and q.k < mant then ' tiny convergent
      df = q.t1 * pw2(q.k)
      if int(df + 0.5) <> 1 then
         print ", ";p.t1 * pw2(p.k);"/";df;
      end if
   else
      print ", ";p.t1;"b+";str(p.k);"/";
      print str(q.t1);"b+";str(q.k);
   end if
end if
end sub

'print partial quotient
#macro printqi
   if i then
      print ","; .q(i);
   else
      print str(.q(i));
   end if
#endmacro

'no truncate, clear or normalize
sub rawcf (byref s as string, byref c as cf, byval sw as integer) export
dim i as integer
dim as cnv p = (0,0,1,0)
dim as cnv q = (0,1,0,0)

with c
   for i = 0 to .w - 1
      nxtcnv p, .q(i) '                Cf recurrence
      nxtcnv q, .q(i)
   next i

   print s; " [";
   if sw and (.w > 24) then '          Cf head and tail only
      for i = 0 to 9
         printqi
      next i
      print " ...  ";
      for i = .w - 10 to .w - 1
         printqi
      next i

   else
      for i = 0 to .w - 1
         printqi
      next i
   end if
   print "] "; i

   prntcnv p, q, 1 '                   fp quotient
   prntcnv p, q, 0 '                   convergent
   print
end with
end sub

'print Cf_c and decimal expansion
sub outcf (byref s as string, byref c as cf) export
if (Verb and 1) <> 0 then
   rawcf s, c, 0
   exit sub
end if
dim as integer i, sw
dim as cfa n
dim as long d
dim as double ip
dim as cnv p = (0,0,1,0)
dim as cnv q = (0,1,0,0)
dim prc as cnv = c.s

if c.r > 0 then
   print "illegal argument outcf: transient Cf"
   exit sub
end if

n.ini
d = iif(c.q(0) < 0,-1, 1)
tset n, 0,d,0,0, 0,0,0,1

n.x = c
with n.x
   if .w = Cfx then .w -= 1 '          truncate
   if .w > Vwx then .w = Vwx '         view-port
  .q(.w) = Inf

   sw = 0
   for i = 1 to .w - 1
      sw or= .q(i) <= 0
   next i
   if sw then Cidem n.x '        clear zeros and negative terms

   for i = 0 to .w - 1
      d = .q(i)
      nxtcnv p, d '                    Cf recurrence
      nxtcnv q, d
      if cmpcnv(q, prc) = 1 then
        .w = i: exit for
      end if
   next i

   if .w > 1 then
      if .q(.w - 1) = 1 then '         normalize
        .q(.w - 2) += 1: .w -= 1
      end if
   end if
  .q(.w) = Inf

   print s; " [";
   for i = 0 to .w - 1
      printqi
   next i
   print "] "; i
end with

if q.t1 then
   if bigcnv(q) then
      print iif(c.q(0) < 0,"-"," ");

      d = 0: i = p.k - q.k
      if i > 0 then '                  big integer part ?
         ip = abs(p.t1 / q.t1)
         d = int(ip * (pw2(i) / Maxd))
      end if
      ip = 0: trlon = 0
      for i = 0 to d * 2
         ip += nextlfd(n)
      next i
      print str(ip) + "."; '           integer part

      do
         for i = 1 to 3 step 2 '       decimal expansion
            'n.v(i) *= 10
            mpz_mul_ui  n.v(i), n.v(i), 10
            'exchange
            swap n.u(i), n.v(i)
         next i
         sw = n.x.q(n.x.r) = Inf
         d = nextlfd(n)
         if d > Inf then print str(d);
      loop until sw or (d = Inf)
      print : trlon = -1

   else
      prntcnv p, q, 1 '                fp quotient
      prntcnv p, q, 0 '                tiny convergent
      print
   end if

else
   print " Inf"
end if
off8 -= 1
end sub

'print the first t (small) convergents to Cf_c,
'set flag to print intermediate convergents
sub cf2q (byref g as string, byref c as cf, byval t as integer,_
 byval fl as integer) export
const as double Lint = 2.0^63 - 1 '    max signed longint
dim as double h
dim as string pq = "1/ 0"
dim as integer r, s, sw, i = 0, j = 1
dim as longint p(1) = {0, 1}
dim as longint q(1) = {1, 0}
dim as long d

print g
for r = 0 to t
   d = c.q(r)
   if d > Inf then
      print d; chr(9); pq

      if fl then '                     intermediate convergents
         sw = (d + 1) shr 1
         for s = sw to d - 1
            if s > sw then print ",";
            print " "; p(i) + s * p(j);
            print "/"; q(i) + s * q(j);
         next s
         if d > 1 and r < t then print
      end if

      h = p(i) + d * cdbl(p(j)) '      fp check
      sw = abs(h) > Lint
      h = q(i) + d * cdbl(q(j))
      sw or= abs(h) > Lint

      p(i) += d * p(j) '               Cf recurrence
      q(i) += d * q(j)
      pq = str(p(i)) + "/ " + str(q(i))

      if sw then
         pq = "break": exit for '      longint overflow
      end if

      swap i, j
   else

      exit for
   end if
next r
print chr(9); pq

if q(j) then print p(j) / q(j)
end sub

'input conversion
' *****************************************************************************
'input continued fraction "q0,q1,q2,..."
sub readcf (byref c as cf, byval w as zstring ptr) export
dim as string s = trim(*w) + ","
dim as integer k = len(s), j, i = 1, v = 0
do
   j = instr(i, s, ",")
   if v < Cfx then '                   read partial quotient
      c.q(v) = valint(mid(s, i, j - i))
   else
      print "Cf truncated: readcf"
      exit do
   end if
   v += 1: i = j + 1
loop until i > k
c.q(v) = Inf
c.r = 0: c.w = v
c.s = cnvmax
end sub

#macro roundcf
if sw then
   v = v0
   if v = 0 then
      print "overflow: cvcf"
   elseif v > 1 then
      if c.q(v - 1) = 1 then '         rounding
         c.q(v - 2) += 1: v -= 1
      end if
   end if
end if
#endmacro

'convert real x to Cf_c
sub cvcf overload (byref c as cf, byval x as double) export
const as double meps = 2.0^-52 '       machine epsilon
dim as cnv q = (0,1,0,0)
dim as double d
dim as integer sw = 0, t = 0, v0 = 0, v = 0

do while v < Cfx
   t = 0: v0 = v
   while abs(x) > Maxd
      t += 1
      sw = t > mexp or v > (Cfx - 2)
      if sw then exit do
      d = iif(x > 0, Maxd, 1 - Maxd)
      c.q(v) = d: v += 1
      c.q(v) = 0: v += 1
      nxtcnv q, d
      nxtcnv q, 0
      x -= d
   wend
   d = int(x) '                        d = floor(x)

   nxtcnv q, d '                       convergent
   if bigcnv(q) then
      if d = 1 then c.q(v - 1) += 1 '  rounding
      exit do
   end if

   c.q(v) = d: v += 1

   x -= d
   if x < meps then exit do
   x = 1 / x
loop
roundcf
c.q(v) = Inf
c.r = 0: c.w = v
c.s = cnvmax
end sub

'convert mpz-rational n.u(i) / n.v(i) to Cf_c
sub mpq2cf (byref c as cf, byref n as cfa, byval i as integer) export
dim as integer sw = 0, t = 0, v0 = 0, v = 0
dim as long d

do while v < Cfx
   if mpz_sgn(n.v(i)) = 0 then exit do

   d = zquot(n, i) '                   Euclidean division

   sw = (d = Inf)
   if d = Maxd then '                  piecewise transmission
      t += 1: sw or= t > mexp
   elseif d > 0 then
      t = 0: v0 = v + 1 '              reset
   end if
   if sw then exit do

   c.q(v) = d: v += 1 '                partial quotient

   if d > 0 then
      mpz_submul_ui  n.u(i), n.v(i), d
   elseif d < 0 then
      mpz_set_si  tmp, d
      mpz_submul  n.u(i), n.v(i), tmp
   end if
   swap n.v(i), n.u(i)
loop
sw or= (v = Cfx)
roundcf
c.q(v) = Inf
c.r = 0: c.w = v
c.s = cnvmax
end sub

sub cvcf (byref a as cf, byval w as zstring ptr) export
dim n as cfa
dim as string g, gd
dim as integer i, j, k, sw

if instr(*w, ",") then '               "q0,q1,q2,..."
   readcf a, w
   exit sub
end if

n.ini
g = trim(*w)
i = instr(g, "/")
j = instr(g, ".")
if i then '                            rational
   gd = mid(g, i + 1) '                denominator,
   g = left(g, i - 1) '                numerator
   sw = mpz_set_str(n.u(0), strptr(g), 10)
   sw or= mpz_set_str(n.v(0), strptr(gd), 10)

elseif j then '                        decimal fraction
   k = len(g) - j
   mid(g, j, 1) = " "  '               erase '.'
   sw = mpz_set_str(n.u(0), strptr(g), 10)
   mpz_ui_pow_ui  n.v(0), 10, k

else '                                 integer
   sw = mpz_set_str(n.u(0), strptr(g), 10)
   mpz_set_ui  n.v(0), 1
end if

if sw then
   print " illegal argument cvcf: "; *w
   empty(a)
else
   mpq2cf a, n, 0
end if
off8 -= 1
end sub

'convert longint quotient x / y to Cf_c
sub cvcf (byref c as cf, byval x as longint, byval y as longint) export
cvcf c, str(x) + "/" + str(y)
end sub

'input polynomial coefficients "a_n,..., a_0"
sub readcoef (byref f as poly, byval w as zstring ptr) export
dim as string s = trim(*w) + ","
dim as integer k = len(s), j, i = 1, v =-1
do
   v += 1: j = instr(i, s, ",")
   if v <= dmx then '                  read coefficient
      f.a(v) = valint(mid(s, i, j - i))
   else
      print "poly truncated: readcoef"
      exit do
   end if
   i = j + 1
loop until i > k

i = (v - 1) shr 1
for j = 0 to i
   swap f.a(j), f.a(v - j)
next j
f.n = v
end sub

' *****************************************************************************
'convert Cf_a to mpz-rational n.u(1) / n.v(1)
sub cf2mpq (byref n as cfa, byref a as cf) export
dim d as long
n.x = a
tset n, 0,a.q(0),0,1, 0,1,0,0
n.x.r = 1
do
   d = popd(n.x) '                     Cf recurrence
   do_absorb d, n, 1
loop until d = Inf

with n
   if mpz_sgn(.v(1)) < 0 then
      mpz_abs .v(1), .v(1)
      mpz_neg .u(1), .u(1)
   end if
   if trlon then truncate n '          limit register size
   mpz_set_ui .u(3), 0
   mpz_set_ui .v(3), 1
end with
end sub

' *****************************************************************************
'convergence check, set print frequency m to some 2^k - 1
sub cnvchck (byref c as cf, byval t as integer, byval m as integer) export
if t < 1 then exit sub
if (c.w and m) = 0 then
   print using " w,t  & / &  #.###"; c.w; t; c.w / t
end if
end sub

sub rate (byval t as integer) export
if t > 0 then
   print using "Cfx &, t &, q ##.###"; Cfx; t; Cfx / t
end if
end sub

