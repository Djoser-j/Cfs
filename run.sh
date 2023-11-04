#!/bin/bash
# ---------------------------------------------------------------
# Command line: ./run.sh <demo name> without .bas

# run this script in the Cfs base directory
# ---------------------------------------------------------------

# executables directory
tool=./modules/bin
# working directory
work=./workdir

echo "Run demo: $1"
demo=$1

echo
if [ ! -f $tool/$demo.exe ]; then
  echo  "not found:" $tool/$demo.exe
  read -p " press enter to quit"

else
  $tool/$demo.exe > $work/$demo.log
  more $work/$demo.log
fi
echo
