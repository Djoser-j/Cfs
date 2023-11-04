#!/bin/bash
# ---------------------------------------------------------------
# Run all cf demo programs

# run this script in the Cfs base directory
# ---------------------------------------------------------------

# executables directory
tool=./modules/bin
# working directory
work=./workdir

echo
if ! ls $tool/*.exe &>/dev/null; then
  echo " demo.exes not found,"
  read -p " press enter to quit"
  echo
  exit
fi

echo  running all cf demos...

rm $work/*.log &>/dev/null

for n in $tool/*.exe; do
  infile=$(basename -s .exe $n)
  $n > $work/$infile.log
done
echo
