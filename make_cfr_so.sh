#!/bin/bash
# ---------------------------------------------------------------
# Compile dynamic cfr library

# run this script in the Cfs base directory
# ---------------------------------------------------------------

# library names
lib1=cfr_arith
lib2=cfr_math

# library directory
ldir=./library
# path to the include file
tool=./modules

bins=$tool/bin

echo
echo  compiling $lib2 DLL

opts="-dll -O 2 -i $tool -w pedantic -x"

fbc $ldir/$lib1.bas $ldir/$lib2.bas $opts $bins/lib$lib2.so >compile.log
