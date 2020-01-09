#!/bin/bash

if [ $# -lt 7 ]
then
  echo "Usage: $0 <inner radius> <outer radius> <separation> <x0> <y0> <z0> <colour>"
  exit
fi
lattice -L `echo "3*$2" | bc -l` -r `echo "0.5*$3" | bc -l` -t fcc -c $7 -f $3 |
  grep -v "#" |
  awk '($1**2 + $2**2 + $3**2 >= '$1'**2) && ($1**2 + $2**2 + $3**2 < '$2'**2) && (NF == 5){printf "%f %f %f %f %d\n", $1 + '$4', 1*$2 + '$5', 1*$3 + '$6', 1*$4, $5}'
