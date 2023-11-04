' *****************************************************************************
'Subject: Include file for continued fraction arithmetic.
'Author : Djoser.j.Spacher
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#ifdef __FB_LINUX__
  #include "gmp.bi"
   type mpz_t as __mpz_struct
#else
  #include "gmp_short.bi"
#endif


'you can change this and recompile:

const Cfx = 101 '                      fixed Cf length + 1

const Rax = Cfx shr 1 '             Cf is reckoned quasi-irrational if longer

#define trees '                      define to use module cf_trees.bas

#ifdef trees
const Regs = 500
#else
const Regs = 31 '                      #state registers
#endif
const Rats = 10 '                      #mpz quotient pairs


'do not change
' *****************************************************************************
const cfb = 32 '                       partial quotient bit size
const mexp = cfb + 3 '                 declare pq infinite

const Maxd = &h3ffffffe '              maximum long pq digit 2^30 - 2
const  Inf = -Maxd - 1 '               infinity  1 - 2^30
const Hold = -Maxd - 2 '               wait state -2^30

const dmx = 20         '               highest poly degree
const mant = 53        '               float mantissa bits
const as double Flint = 2.0^53 - 1 '   max integral float
const Ulng = (clngint(1) shl 32) - 1 ' max unsigned long

'checks
' *****************************************************************************
if Cfx < 46 then
'fit the largest long Fibonacci quotient F45 / F46
   print "Cfx too small:"; Cfx : end
'set view-port size to display less than 46 digits
end if

'is Cf_a zero ?
#define is0(a) (a.q(0) = 0 and a.q(1) = Inf)
'is a less than or equal to 0 ?
#define leq0(a) ((a.q(0) < 0) or is0(a))
'terminate loop ?
#define tm(a) (a.q(a.w) = Inf)
'irrational Cf ?
#define qirr(a) (a.w > Rax)

'set a []
#define empty(a) a.q(0) = Inf: a.w = 0

'simple loops
' *****************************************************************************
#macro dobf(t, s)
   t.ini
   do : push t, nextbfd(s) '          (axy + bx + cy + d) / (exy + fx + gy + h)
   loop until tm(t)
#endmacro
#macro dolf(t, s)
   t.ini
   do : push t, nextlfd(s) '          (bx + d) / (fx + h)
   loop until tm(t)
#endmacro

'custom data types
' *****************************************************************************
type cnv '                            Scaled convergent:
   as double s, t0, t1 '               terms
   as integer k '                      base 2 exponent
end type

'maximum convergent
dim shared cnvmax as const cnv = (0, 0, 1, 1 shl 30)

type pair
declare sub ini ()
   as mpz_ptr a, b '                   mpz parameter pair
   as integer id '                     pair id
end type

type cf '                             Continued fraction:
declare constructor ()
declare sub ini ()
   as long q(any) '                    partial quotients
   as integer r, w '                   read and write pts
   as cnv s '                          convergent size
end type

type cfa '                            Cf arithmetic structure:
declare constructor ()
declare constructor (byval fl as integer)
declare sub ini ()
   as mpz_ptr u(3), v(3) '             mpz 2 x 2 matrices
   as cf x, y '                        continued fractions
   as integer id '                     set id
end type

type poly '                           Polynomial:
declare constructor ()
   as long a(any) '                    coefficients
   as integer n '                      degree
end type

'function types, lazy evaluation
' *****************************************************************************
'root finding with Lagrange's 1767 Cf method
type polr_t
declare constructor ()
declare sub ini (byref g as poly)
declare function f (byref c as cf) as long
   as mpz_ptr h(any)
   as integer n, t0, t
end type

#ifndef __LIB

'solve ax^2 + bx - c = 0
type quad_t
declare sub ini (byval a as long, byval b as long, byval c as long)
declare function f () as long
   as long d, rD, p(1), q(1)
   as integer i, j, t
end type

'pseudo-random partial quotients
type rnd_t
declare sub ini ()
declare sub seed (byval s as ushort)
declare function f () as long
   as quad_t cgenr
end type

'cf input
' *****************************************************************************
type sqrt_2
declare sub ini ()
declare function f () as long
   as integer t
end type

type sqrt_t
declare sub ini (byref a as cf)
declare function f () as long
   as cfa m=cfa(0), n
   as sqrt_2 sqrt2
   as pair q
   as long d
   as integer k, fl, sw, t
end type

type log_t
declare sub ini (byref a as cf)
declare function f () as long
   as cfa m=cfa(0), n=cfa(0), s
   as pair q
   as long d
   as integer k, fl = 0, t2, t
end type

type logr_t extends log_t
declare constructor ()
end type

type exp_t
declare sub ini (byref a as cf)
declare function f () as long
   as cfa m=cfa(0), n(any)
   as pair q
   as long d
   as integer k, t
end type

type tanh_t
declare sub ini (byref a as cf)
declare function f () as long
   as cfa m=cfa(0), n(any), s
   as pair q
   as long d
   as integer i = 1, fl = 0, k, t
end type

type sinh_t extends tanh_t
declare constructor ()
end type

type cosh_t extends tanh_t
declare constructor ()
end type

type tan_t extends tanh_t
declare constructor ()
end type

type sin_t extends sinh_t
declare constructor ()
end type

type cos_t extends cosh_t
declare constructor ()
end type

type pi_t
declare sub ini (byval k as long)
declare function f () as long
   as cfa m=cfa(0), n=cfa(0), s
   as long d
   as integer t2, t
end type

type atan_t
declare sub ini (byref a as cf)
declare function f () as long
   as cfa m=cfa(0), n
   as pi_t cpi
   as pair q
   as long d
   as integer k, sw, t
end type

'eager evaluation
' *****************************************************************************
type sqrtN_t
declare sub f (byref r as cf, byref a as cf)
   as cfa m, n
end type

type asin_t
declare sub f (byref r as cf, byref a as cf)
   csqrt as sqrtN_t
   catan as atan_t
   cpi as pi_t
   n as cfa
   fl as integer = 0
end type

type acos_t extends asin_t
declare constructor ()
end type

type atanh_t
declare sub f (byref r as cf, byref a as cf)
   clog as logr_t
   n as cfa
end type

type asinh_t
declare sub f (byref r as cf, byref a as cf)
   csqrt as sqrt_t
   clog as log_t
   n as cfa
   fl as integer = 0
end type

type acosh_t extends asinh_t
declare constructor ()
end type

type gd_t
declare sub f (byref r as cf, byref a as cf)
   catan as atan_t
   cexp as exp_t
   cpi as pi_t
   n as cfa
end type

type gdinv_t
declare sub f (byref r as cf, byref a as cf)
   clog as log_t
   csin as sin_t
   ccos as cos_t
   n as cfa
end type

type pow_t
declare sub f (byref r as cf, byref a as cf, byref e as cf)
   clog as log_t
   cexp as exp_t
   n as cfa
end type

type agm_t
declare sub f (byref r as cf, byref g as cf, byref h as cf)
   csqrt as sqrtN_t
   as cfa m, n
   b as cf
end type

'lazy evaluation with rational input
' *****************************************************************************
type sin3_t
declare sub ini (byval k as long)
declare function f () as long
   as cfa h=cfa(0), m, n, s
   as pair q
   as long d
   as integer t
end type

type tanhr_t
declare sub ini (byval k as short)
declare function f () as long
   as quad_t cquad
   as cfa m=cfa(0), n
   as pair q
   as integer t
end type

type exppi_t
declare sub ini (byval x as long, byval y as long)
declare function f () as long
   as cfa n=cfa(0)
   as long a2, b
   as integer t
end type

type tanpi_t
declare sub ini (byval x as long, byval y as long)
declare function f () as long
   as cfa n=cfa(0)
   as long a, b
   as integer t
end type

type sinpi_t
declare sub ini (byval x as long, byval y as long)
declare function f () as long
   as cfa m=cfa(0), n
   as long a, b, d
   as integer fl = 1, t
end type

type cospi_t extends sinpi_t
declare constructor ()
end type

type besseli_t
declare sub ini (byval v as integer, byval x as long, byval y as long)
declare function f () as long
   as cfa n=cfa(0)
   as long a2, b, u
   as integer i = 1, t
end type

type besselj_t extends besseli_t
declare constructor ()
end type


'input conversion
' *****************************************************************************
'convert real x to Cf_c
declare sub cvcf overload (byref c as cf, byval x as double)
'convert longint quotient x / y to Cf_c
declare sub cvcf (byref c as cf, byval x as longint, byval y as longint)
'input decimal fraction "a.b0b1b2..." or rational "x / y"
declare sub cvcf (byref c as cf, byval w as zstring ptr)
'input continued fraction "q0,q1,q2,..."
declare sub readcf (byref c as cf, byval w as zstring ptr)

'input polynomial coefficients "a_n,..., a_0"
declare sub readcoef (byref f as poly, byval w as zstring ptr)

'Cf output
' *****************************************************************************
'print Cf_c and decimal expansion
declare sub outcf (byref s as string, byref c as cf)
'print the first t (small) convergents to Cf_c,
'set flag to print intermediate convergents
declare sub cf2q (byref g as string, byref c as cf,_
 byval t as integer, byval fl as integer)
'print raw Cf, set switch for head and tail only
declare sub rawcf (byref s as string, byref c as cf, byval sw as integer)

' *****************************************************************************
'assign x and y to pair
declare sub tset overload (byref p as pair, byval x as long, byval y as long)
'assign (a b c d) / (e f g h) to tensor and initialize Cf x and y
declare sub tset (byref n as cfa,_
 byval a as long, byval b as long, byval c as long, byval d as long,_
 byval e as long, byval f as long, byval g as long, byval h as long)

'push digit in Cf_c queue
declare sub push (byref c as cf, byval d as long)
'pop Cf digit
declare function popd (byref c as cf) as long
'next bi-linear form digit z(x, y)
declare function nextbfd (byref n as cfa) as long
'next linear form digit z(x)
declare function nextlfd (byref n as cfa) as long

'return least precision in prc
declare sub least overload (byref prc as cnv, byref c as cf)
'set a.s to least precision
declare sub least (byref a as cf, byref prc as const cnv)
declare sub least (byref a as cf, byref c as cf)

'compare Cf's, mismatch pointer t,
'return a < b: -1, a = b: 0, a > b: 1
declare function comp (byref t as integer,_
 byref a as cf, byref b as cf) as integer
'Cf_a:= -Cf_a
declare sub Cneg (byref a as cf)
'Cf_a:= 1 / Cf_a
declare sub Cinv (byref a as cf)
'identity
declare sub Cidem (byref a as cf)

'other functions
' *****************************************************************************
'convert Cf_a to mpz-rational n.u(1) / n.v(1)
declare sub cf2mpq (byref n as cfa, byref a as cf)
'convert mpz-rational n.u(i) / n.v(i) to Cf_a
declare sub mpq2cf (byref a as cf, byref n as cfa, byval i as integer)

'convergence check, set print frequency m to some 2^k - 1
declare sub cnvchck (byref c as cf, byval t as integer, byval m as integer)
'convergence speed
declare sub rate (byval t as integer)

'generalized Cf for transcendental functions,
'set switch to evaluate tanh (sw = 0) or atan (sw = 1)
declare function gencf (byref n as cfa, byref q as pair,_
 byref t as integer, byval sw as integer) as long
'absorb trafo
declare sub do_genr (byref n as cfa, byval sw as integer)
'linear form digit z(x)
declare function linfd (byref n as cfa) as long

'set or get
' *****************************************************************************
'switch to verbose mode, set
'bit 1 (1) full Cf, ignoring precision
'bit 2 (2) print iteration steps
'bit 3 (4) linear form throughput streams
'bit 4 (8) bi-linear throughput streams
'bit 5 (16) print binary exponential format
declare sub setsw (byval sw as integer)
'return the current mode
declare function getsw () as integer
'set view-port
declare sub setvw (byval sw as integer)
'return the current size
declare function getvw () as integer

'clear mpz stack; use often to reclaim register space
declare sub clrs ()

#endif

'variables shared between both library modules
' *****************************************************************************
common shared pw2() as double '        powers of two
common shared tmp as mpz_ptr '         global mpz pointer
'offset for assigning the next 2 or 8-variable mpz set
common shared as integer off2, off8

