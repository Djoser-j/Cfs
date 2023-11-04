#!/bin/bash
# ---------------------------------------------------------------
# Install shared library

# run in the base directory after compiling the library
# ---------------------------------------------------------------

# shared library directory
tool=$PWD/modules/bin
# FreeBasic config file
file=/etc/ld.so.conf.d/freebasic.conf

if ! grep -qx "$tool" "$file"; then
  echo Add library path
  echo "$tool" | sudo tee -a "$file"
fi

echo Update cache
sudo ldconfig
echo
