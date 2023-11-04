#!/bin/bash
# ---------------------------------------------------------------
# Compile all cf demo programs

# run this script in the Cfs base directory
# ---------------------------------------------------------------

# source files directory
tool=./modules

bins=$tool/bin

ls $bins/lib*.so &>/dev/null || ./make_cfr_so.sh
echo
echo  compiling all cf demos...

opts="-O 2 -p $bins -s console -w pedantic -x"

for f in $tool/*.bas;
do
  fbc $f $opts "${f%.*}".exe >>compile.log
done

mv $tool/*.exe $bins >/dev/null

echo
read -p "press enter to continue"
