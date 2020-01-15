#!/bin/sh
if [ $# -lt 1 ]
then
  echo "Usage: $0 <simulation directory>"
  exit
fi

PARAM_WALLCOLOUR=255 # Particle type (colour) of the QCM resonator
PARAM_MASS=0.52 # Particle mass (calibrated from fit for a given resolution)
PARAM_KSPRING=0.2 # Strength of the fixed bonds
PARAM_DELTA=12 # Penetration of the oscillating flow

for dir in $@
do

SIMDIR=$dir

### Get parameters from simulation files ###
Lx=$(grep "celldimension" $SIMDIR/data.main | awk '{print $2}') # Box size (x direction)
Ly=$(grep "celldimension" $SIMDIR/data.main | awk '{print $3}') # Box size (y direction)
AREA=$(echo "$Lx * $Ly" | bc -l) # QCM surface area
NSTEPS=$(grep "numsteps[^R]" $SIMDIR/data.main | awk '{print $2}') # Number of simulation steps
TIMESTEP=$(grep "dt" $SIMDIR/data.main | awk '{print $2}') # Simulation time step
SAMPLING=$(grep "samplefreq" $SIMDIR/data.main | awk '{print $2}') # Sampling step
TMAX=$(echo "$NSTEPS*$TIMESTEP" | bc -l) # Simulation time
LAG=$(echo "$TIMESTEP*$SAMPLING" | bc -l) # Time lag between frames
FREQUENCY=$(grep "stokesFlowParams" $SIMDIR/data.main | awk '{print $2}') # Forcing frequency
F0=$(grep "stokesFlowParams" $SIMDIR/data.main | awk '{print $3}') # Forcing amplitude
omega=$(echo "pi = 4*a(1); 2*pi*$FREQUENCY" | bc -l) # Forcing angular frequency
eta=$(grep "viscosity" $SIMDIR/data.main | awk '{print $2}') # Shear viscosity
delta=$PARAM_DELTA # Penetration of the oscillating flow
FIT_T0=$(echo "0.2*$TMAX" | bc -l) # Time at which to start fit

# Count the number of particles in the QCM resonator
WN=$(awk 'BEGIN{n = 0} ($4 == '$PARAM_WALLCOLOUR'){n = n + 1} END{print n}' $SIMDIR/initial_configuration*.dat)

### x coordinate of the centre-of-mass motion of the wall ###
if [ ! -e $SIMDIR/wall.dat ]
then
  PARTICLE_FILE=$(ls $SIMDIR/*particles | rev | cut -d '/' -f 1 | rev)
  ./centre_of_mass.sh $SIMDIR/$PARTICLE_FILE $PARAM_WALLCOLOUR $LAG | awk '{print $1, $2}' > $SIMDIR/wall.dat
fi

### Fit centre of mass trajectory to function x(t) ###
cat << EOF > fit.gnuplot.tmp # Create gnuplot fitting script
# Generic sinusoidal function with the right frequency
x(t) = xr*cos($omega * t) + xi*sin($omega * t) + B
# Fit parameters
fit [$FIT_T0:] x(x) "$SIMDIR/wall.dat" via xr, xi, B
# Output results
set print "wall.fit.tmp"
print xr, xi, B
# plot "$SIMDIR/wall.dat", x(x) lw 2
# pause -1
EOF

# Run gnuplot fit
gnuplot fit.gnuplot.tmp 2> /dev/null

# x(t) phasor amplitude (real and imaginary parts)
xr=$(awk '{printf "%2.8f", $1}' wall.fit.tmp)
xi=$(awk '{printf "%2.8f", $2}' wall.fit.tmp)

# clean up
rm fit.gnuplot.tmp wall.fit.tmp

### Calculate the impedance of the analyte ###

# Complex amplitude of the x(t) phasor squared
x2=$(echo "${xr}^2 + ${xi}^2" | bc -l)

# Force density on the QCM due to forcing
Fr=$(echo "0" | bc -l)
Fi=$(echo "$F0 * $WN / $AREA" | bc -l)

# Shear stress caused by the fluid on the QCM resonator
shearr=$(echo "2 * $eta * $omega * ( $xr - $xi ) / $delta" | bc -l)
sheari=$(echo "2 * $eta * $omega * ( $xr + $xi ) / $delta" | bc -l)

# Fluid inertia
inertiar=$(echo "-$PARAM_MASS * ${omega}^2 * $xr * $WN / $AREA" | bc -l)
inertiai=$(echo "-$PARAM_MASS * ${omega}^2 * $xi * $WN / $AREA" | bc -l)

# Spring force
springr=$(echo "-$PARAM_KSPRING * $xr * $WN / $AREA" | bc -l)
springi=$(echo "-$PARAM_KSPRING * $xi * $WN / $AREA" | bc -l)

# Stress due to the analyte (indirect measurement)
sigmar=$(echo "$inertiar -($Fr + $shearr + ${springr})" | bc -l)
sigmai=$(echo "$inertiai -($Fi + $sheari + ${springi})" | bc -l)

# Complex impedance of the analyte (through indirect measurement: calibrate mass before using xw)
Zr=$(echo "scale = 6; (${sigmar}*${xi} - ${sigmai}*${xr})/(${omega}*${x2})" | bc -l)
Zi=$(echo "scale = 6; (${sigmar}*${xr} + ${sigmai}*${xi})/(${omega}*${x2})" | bc -l)
AR=$(echo "scale = 6; -2 * $Zr / $Zi" | bc -l)

### Fit velocity gradient to numerical data  ###
# Stress due to analyte (numerical measurement)
if [ -e $SIMDIR/*gradv ]
then

  # Get data filename
  GRADV_FILE=$(ls $SIMDIR/*gradv | rev | cut -d '/' -f 1 | rev)

# Create gnuplot fitting script
cat << EOF > fit.gnuplot.tmp
# Generic sinusoidal function with the right frequency
v(t) = vr*cos($omega * t) + vi*sin($omega * t) + C
gradv(t) = gradvr*cos($omega * t) + gradvi*sin($omega * t) + B
# Fit parameters
fit [$FIT_T0:] gradv(x) "$SIMDIR/$GRADV_FILE" u ($TIMESTEP*\$1):2 via gradvr, gradvi, B
fit [$FIT_T0:] v(x) "$SIMDIR/$GRADV_FILE" u ($TIMESTEP*\$1):3 via vr, vi, C
# Output results
set print "gradv.fit.tmp"
print gradvr, gradvi, B, vr, vi, C
EOF

  # Run gnuplot fit
  gnuplot fit.gnuplot.tmp 2> /dev/null

  # Velocity phasor amplitude (real and imaginary parts)
  vr=$(awk '{printf "%2.8f", $4}' gradv.fit.tmp)
  vi=$(awk '{printf "%2.8f", $5}' gradv.fit.tmp)
  v2=$(echo "$vr^2 + $vi^2" | bc -l)

  # gradv phasor amplitude (real and imaginary parts)
  sigmar_2=$(awk '{printf "%2.8f", '$eta'*($1 - ('$xr' - '$xi')*'$omega'/'$delta')}' gradv.fit.tmp)
  sigmai_2=$(awk '{printf "%2.8f", '$eta'*($2 - ('$xr' + '$xi')*'$omega'/'$delta')}' gradv.fit.tmp)

  # Impedance and acoustic ratio (second method)
  Zr_2=$(echo "scale = 6; (${sigmar_2}*${xi} - ${sigmai_2}*${xr})/(${omega}*${x2})" | bc -l)
  Zi_2=$(echo "scale = 6; (${sigmar_2}*${xr} + ${sigmai_2}*${xi})/(${omega}*${x2})" | bc -l)
  AR_2=$(echo "scale = 6; -2 * $Zr_2 / $Zi_2" | bc -l)

  # Impedance and acoustic ratio (third method)
  Zr_3=$(echo "scale = 6; (${sigmar_2}*${vr} + ${sigmai_2}*${vi})/$v2" | bc -l)
  Zi_3=$(echo "scale = 6; (${sigmai_2}*${vr} - ${sigmar_2}*${vi})/$v2" | bc -l)
  AR_3=$(echo "scale = 6; -2 * $Zr_3 / $Zi_3" | bc -l)

  # Clean up
rm fit.gnuplot.tmp gradv.fit.tmp
else
  Zr_2="N/A"
  Zi_2="N/A"
  AR_2="N/A"
  Zr_3="N/A"
  Zi_3="N/A"
  AR_3="N/A"
fi

# echo "x0: $xr + $xi i"
# echo "F0: $Fr + $Fi i"
# echo "Shear: $shearr + $sheari i"
# echo "Fluid inertia: $inertiar + $inertiai i"
# echo "Spring: $springr + $springi i"
# echo "Stress: $sigmar + $sigmai i, $sigmar_2 + $sigmai_2 i"

# Output results
#   Directory - Method 1 -  - Method 2 -        - Method 3 -
echo "$SIMDIR $Zr $Zi $AR   $Zr_2 $Zi_2 $AR_2   $Zr_3 $Zi_3 $AR_3"

done
