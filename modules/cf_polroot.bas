' *****************************************************************************
'Subject: Continued fraction expansion of real polynomial roots.
'Code   : FreeBasic 1.06.0 with GMP library 6.1.1

#include "cfr_lib.bi"
#inclib "cfr_math"

' *****************************************************************************
'polynomial roots
sub polroots
dim as cf a, b, c
dim as poly f
dim as polr_t cpolroot
dim as quad_t cquad

print "polynomial roots" : ?

   readcoef f, "8,0,-6,1" '            trisection
   c.ini: cvcf a, "-1" '               cos((8*pi/3)/ 3)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a1)=0", c
   c.ini: cvcf a, "0,5" '              cos((4*pi/3)/ 3)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a2)=0", c
   c.ini: cvcf a, "0,1" '              cos((2*pi/3)/ 3)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a3)=0", c

   print
   readcoef f, "1,0,-2,-5" '           Newton 1669
   c.ini: cvcf a, "2"
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   'Brent, vd Poorten, te Riele (1996): Table 2
   readcoef f, "1,0,-8,-10"
   c.ini: cvcf a, "3"
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   readcoef f, "1, 8,14,-8,-23" '      sqrt(2) + sqrt(3) - 2
   cpolroot.ini f
   c.ini: empty(a)
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   readcoef f, "1,-8,-12,184,-178,-664,580,744,-71"
   cpolroot.ini f
   c.ini: empty(a)
   do : push c, cpolroot.f(a) '        sqrt(5) - sqrt(3) - sqrt(2) + 1
   loop until tm(c)
    outcf "f(x)=0", c

   readcoef f, "1,0,0,0,-20,0,0,0,-666,0,0,0,-3860,0,0,0,1"
   c.ini: readcf a, "0"
   cpolroot.ini f
   do : push c, cpolroot.f(a) '        3^(1/4) - 2^(1/4)
   loop until tm(c)
    outcf "f(x)=0", c

   print
   'Cantor, Galyean, Zimmer (1972)
   readcoef f, "1,0,0,0,0,0,-7,3"
   c.ini: readcf a, "-2"
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a1)=0", c
   c.ini: readcf a, "0,2"
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a2)=0", c
   c.ini: readcf a, "1"
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a3)=0", c

   print
   readcoef f, "1,0,-10667,0,0,0,-3367,0,15782,0,0,0,0,0,0,-10541"
   cpolroot.ini f
   c.ini: readcf a, "-104"
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a1)=0", c
   cpolroot.ini f
   c.ini: cvcf a, "-2,1"
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a2)=0", c
   cpolroot.ini f
   c.ini: readcf a, "103"
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(a3)=0", c

   print
   c.ini: cquad.ini 5, 18, 15 '        solve quadratic
   do : push c, cquad.f
   loop until tm(c)
    outcf "5xx + 18x - 15 = 0", c

   c.ini: cquad.ini 349,-500, 0 '      rational to Cf
   do : push c, cquad.f
   loop until tm(c)
    outcf "349x = 500", c

   readcoef f, "343,0,-7,-120"
   c.ini: empty(a)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   readcoef f, "3,-2,3,-2"
   c.ini: empty(a)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   readcoef f, "121166534,-343305231"
   c.ini: empty(a)
   cpolroot.ini f
   do : push c, cpolroot.f(a)
   loop until tm(c)
    outcf "f(x)=0", c

   print
   'big Disc 9007199254711125, small regulator
   c.ini: cquad.ini 31,94906257,13201299
   do : push c, cquad.f
   loop until tm(c)
    outcf "pos.root", c

   c.ini: cquad.ini 31,-94906257,13201299
   do : push c, cquad.f
   loop until tm(c)
   Cneg c
    outcf "neg.root", c

   c.ini: cquad.ini 31,94906247,28508759
   do : push c, cquad.f
   loop until tm(c)
    outcf "pos.root", c

   c.ini: cquad.ini 31,-94906247,28508759
   do : push c, cquad.f
   loop until tm(c)
   Cneg c
    outcf "neg.root", c

   clrs
end sub

' *****************************************************************************
' main
dim as double tim = timer

cls

setvw Cfx '                            view-port size
setsw 0 '                              verbosity switch

polroots

print : ? "timer:"; csng(timer - tim); "s"
system
