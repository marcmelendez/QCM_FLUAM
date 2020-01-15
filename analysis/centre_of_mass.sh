#!/bin/sh
if [ $# -lt 3 ]
then
  echo "Calculate the centre of mass trajectory of particles of a given type."
  echo "Usage: $0 <Filename> <Particle type> <time step between frames>"
  exit
fi

awk 'BEGIN{np = 0; x = 0; y = 0; z = 0; t = '$3';}
     ($5 == '$2'){x = x + $1; y = y + $2; z = z + $3; np = np + 1}
     (NF == 1 && NR > 1){print t, x/np, y/np, z/np; np = 0; x = 0; y = 0; z = 0; t = t + '$3'}' $1
