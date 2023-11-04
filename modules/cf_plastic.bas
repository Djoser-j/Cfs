' *****************************************************************************
'Subject: The smallest Pisot-Vijayaraghavan number
'Code   : FreeBasic 1.10.0 with GMP library 6.2.1

#include "cfr_lib.bi"
#inclib "cfr_math"

dim as integer fl = (getsw and 2) <> 0
dim as integer t = 1
dim as cfa m, n
dim as cf r
m.ini: n.ini

print : ? "solve x^3 - x = 1" : ?

cvcf r, 4/3
while (t shr 1) < Cfx
   n.x = r: n.y = r
   m.x = r: r.ini '                    feed back
   tset n, 1,0,0,0, 0,0,0,1
   tset m, 2,0,0,1, 0,0,3,-1
   do '                                Newton iteration
      push m.y, nextbfd(n) '           y:= X^2
      push r, nextbfd(m) '            (2Xy + 1) / (3y - 1)
   loop until tm(r) or (r.w = t + 3)
   push r, Inf

   if fl then
      print "Pisot_t ="; t
      outcf " X", r
   end if
   t shl= 1
wend

outcf "Le nombre plastique", r '       A072117

system
