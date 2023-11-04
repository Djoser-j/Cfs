#!/bin/bash
# Delete all compiles and logfiles
# run this script in the Cfs base directory

# executables directory
tool=./modules/bin
# working directory
work=./workdir

rm $tool/lib*.so >/dev/null
rm $tool/*.exe >/dev/null

rm $work/*.log >/dev/null

rm ./compare.txt >/dev/null
rm ./compile.log >/dev/null
