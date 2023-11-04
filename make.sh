#!/bin/bash
# ---------------------------------------------------------------
# Command line: ./make.sh <demo name> without .bas

# run this script in the Cfs base directory
# ---------------------------------------------------------------

# source files directory
tool=./modules

bins=$tool/bin

ls $bins/lib*.so &>/dev/null || ./make_cfr_so.sh

echo "Make demo: $1"
demo=$1

echo
if [ ! -f $tool/$demo.bas ]; then
  echo  "Not found:" $tool/$demo.bas
  read -p " press enter to quit"

else
  echo  "compiling demo" $demo.bas

  opts="-O 2 -p $bins -s console -w pedantic -x"

  fbc $tool/$demo.bas $opts $bins/$demo.exe >compile.log
fi
echo
