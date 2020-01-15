#!/bin/bash

# Calculate the shear stress on the QCM wall due to a parallel stationary wall
# as a function of the distance y between the walls

# Stokes flow parameters
rho=1 # Fluid density
omega=1 # Angular frequency of oscillation
delta=1.414213562 # Penetration depth of oscillating flow
u0=1 # Amplitude of velocity

echo "# Ly  # shear stress (real part) # shear stress (imaginary part)"
for y in `seq 0.05 0.05 6`
do
  echo "$rho $delta $omega $u0 $y" | bc -l complex.bc shear_stress.bc
done
