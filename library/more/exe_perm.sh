#!/bin/bash
# Set execute permission on bash scripts

for n in ./*.sh; do
  chmod +x $n
done
