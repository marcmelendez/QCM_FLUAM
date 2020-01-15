#!/bin/bash

### Parameters ###
PARAM_SIMNAME="QCM_sim" # Simulation name
PARAM_NODE="0-9" # GPU node
PARAM_BDNODE="0-9" # CPU node
PARAM_CASE="run6" # Case name

# Experimental parameters
PARAM_LDNAEXP=50 # Length of DNA strand (nm)
PARAM_LPEREXP=50  # Persistence length of DNA strand (nm)
PARAM_DELTAEXP=95 # Penetration length of the oscillating fluid (nm)
PARAM_RADIUSEXP=50 # Liposome radius (nm)
PARAM_LIPOSOMETHICKNESSEXP=10 # Thickness of the liposome membrane (nm)

# Constants
PARAM_DELTA=12 # Penetration length in simulation units (fixes the simulation resolution)
PARAM_CELLSIZE=0.5 # Size of simulation cell in FLUAM.
# PARAM_MF=0.5 # Fluid mass in blob

# Files
PARAM_RUNSCRIPT="run.sh" # Cluster run script for FLUAM
PARAM_BDRUNSCRIPT="runBD.sh" # Cluster run script for Brownian Dynamics
PARAM_PARAMETERFILE="data.main" # Fluam parameter file
PARAM_WALLFILE="wall.tmp" # Wall particle xy coordinates
PARAM_PARTICLEFILE="initial_configuration.dat" # Initial configuration file name
PARAM_BONDFILE="bonds.fluam.dat" # Bond file name
PARAM_ANGULARSPRINGSFILE="bonds3.fluam.dat"
PARAM_BDINITFILE="BD_initial_configuration.dat" # Initial configuration for Brownian dynamics file
PARAM_BDBONDFILE="bonds.nanometre.dat" # Bond file for nanometre BD simulation
PARAM_BDCONFIGURATIONS="BD_configurations.dat" # Brownian dynamics output configurations

# Simulation
PARAM_TIMESTEP=0.002 # Integration time step
PARAM_NPERIODS=20 # Number of QCM oscillation periods
PARAM_NDATAPER=200 # Number of points per period to save

# System
PARAM_LX=32 # Simulation box length
PARAM_LY=32 # Simulation box width
PARAM_LZ=128 # Simulation box height

# Particles
PARAM_BOTTOMWALLC=255 # 30   # Colour for bottom wall
PARAM_TOPWALLC=250 # 30     # Colour for top wall
PARAM_LJEPSILON=0 # Lennard-Jones epsilon parameter
PARAM_EXCESSMASS=0 # Excess mass of particles

# Walls
PARAM_WALLBOND=10.0 # Wall bond parameter (holding the wall at the right z position --no xy force--).
PARAM_WALLAMP=1 # Amplitude of QCM oscillation
PARAM_WALLFORCETMAX=10000000 # Duration of forcing
PARAM_WALLFREQ=0.0025 # Frequency of forcing
PARAM_FBOND=0.2 # Fixed bond parameter for walls (holding each wall particle to its initial position).
                # Set FBOND to 0 (integer zero) to ignore fixed bonds

# Polymers
PARAM_PBOND=100.0 # Polymer bond parameter
PARAM_CHAINC=123123 # 3    # Polymer chain colour

# Liposome
PARAM_NLIPOSOMES=0 # Number of liposomes to simulate
PARAM_LIPOSOMEC=321123 # 5 # Liposome colour
PARAM_LBOND=100.0 # Liposome bond parameter

# Fluid
PARAM_RHO=1 # Fluid density
PARAM_KT=1.0 # Thermal energy
PARAM_THERMALFLUCTUATIONS=N # Set to "N" to ignore thermal fluctuations

# Initial configurations and sampling by Brownian dynamics
PARAM_NEWCONFIGURATION=Y # Generate walls, liposome and bonds
PARAM_RUNBD=N # Run Brownian dynamics simulation (set to Y to generate initial configurations by Brownian dynamics)
PARAM_BDTRELAX=100000 # Relaxation time
PARAM_BDTMAX=200000  # Final time
PARAM_BDDT=0.001   # Integration time step
PARAM_BDSAMPLES=1  # Number of initial configurations
PARAM_BDKT=$PARAM_KT # Thermal energy for the Brownian dynamics simulation

# Shear tensor for Brownian Dynamics simulations
read -r -d '' SHEAR_TENSOR <<'EOF'
  0  0  0
  0  0  0
  0  0  0
EOF

### Insert parameters from scripts here ###

# Script parameters

### Derived parameters ###

# Number of cells in each direction
PARAM_NX=`echo "scale = 0; $PARAM_LX/$PARAM_CELLSIZE" | bc -l`
PARAM_NY=`echo "scale = 0; $PARAM_LY/$PARAM_CELLSIZE" | bc -l`
PARAM_NZ=`echo "scale = 0; $PARAM_LZ/$PARAM_CELLSIZE" | bc -l`
PARAM_SIGMA=$PARAM_CELLSIZE # Excluded volume parameter
PARAM_RADIUS=`echo "0.5*$PARAM_SIGMA" | bc -l` # Particle radius
# PARAM_MASS=`echo "$PARAM_MF + $PARAM_EXCESSMASS" | bc -l` # Particle mass
PARAM_MASS=`echo "pi = 4*a(1); 4*pi*$PARAM_RADIUS^3*$PARAM_RHO/3 + $PARAM_EXCESSMASS" | bc -l` # Particle mass
PARAM_LJCUTOFF=`echo "$PARAM_SIGMA * 1.122462" | bc -l` # Lennard-Jones cut off
# PARAM_WALLPSEPARATION=`echo "4*$PARAM_RADIUS" | bc -l` # Separation between closest particles in the walls
# PARAM_WALLPSEPARATION=`echo "3*$PARAM_RADIUS" | bc -l` # Separation between closest particles in the walls (overlapping blobs)
PARAM_WALLPSEPARATION=`echo "2*$PARAM_RADIUS" | bc -l` # Separation between closest particles in the walls (overlapping blobs)
PARAM_WALLWIDTH=`echo "2*$PARAM_RADIUS" | bc -l`
PARAM_RCUT=`echo "1.5*$PARAM_WALLPSEPARATION" | bc -l` # Cut-off for elastic network model
PARAM_NBEADS=`echo "scale = 0; $PARAM_LDNAEXP*$PARAM_DELTA/($PARAM_DELTAEXP*$PARAM_SIGMA)" | bc -l` # Number of beads in a polymer chain
PARAM_OUTERRADIUS=`echo "$PARAM_RADIUSEXP*$PARAM_DELTA/$PARAM_DELTAEXP" | bc -l` # Liposome outer radius
PARAM_LIPOSOMETHICKNESS=`echo "$PARAM_LIPOSOMETHICKNESSEXP*$PARAM_DELTA/$PARAM_DELTAEXP" | bc -l` # Thickness of the liposome membrane
PARAM_INNERRADIUS=`echo "$PARAM_OUTERRADIUS - $PARAM_LIPOSOMETHICKNESS" | bc -l` # Liposome inner radius
PARAM_ANGULARK=`echo "$PARAM_LPEREXP*$PARAM_DELTA*$PARAM_SIGMA/($PARAM_DELTAEXP*$PARAM_KT)" | bc -l` # Angular spring constant
PARAM_ETA=`echo "pi = 4*a(1); pi*$PARAM_WALLFREQ*$PARAM_DELTA^2*$PARAM_RHO" | bc -l` # Shear viscosity
PARAM_NSTEPS=`echo "nsteps = $PARAM_NPERIODS/($PARAM_WALLFREQ*$PARAM_TIMESTEP); scale = 0; nsteps/1" | bc -l` # Number of time steps to perform
PARAM_SAVEFREQ=`echo "savefreq = 1/($PARAM_NDATAPER*$PARAM_WALLFREQ*$PARAM_TIMESTEP); scale = 0; savefreq/1" | bc -l`
PARAM_SAMPLEFREQ=`echo $PARAM_SAVEFREQ`
PARAM_FRAMESTEP=`echo "$PARAM_SAMPLEFREQ * $PARAM_TIMESTEP" | bc -l`

if [ $PARAM_THERMALFLUCTUATIONS == "N" ]
then
  PARAM_KT=0.0 # Ignore thermal fluctuations
fi

####################################################################

# Adjust harmonic bonds to compensate for bonds in angular bonds
BDBOND=$PARAM_PBOND

# Number of angular springs
NANGULAR_SPRINGS=$(($PARAM_NLIPOSOMES*($PARAM_NBEADS - 2)))

# Move to simulation directory
echo "Creating simulation dir $PARAM_SIMNAME..."
mkdir -p $PARAM_SIMNAME
cd $PARAM_SIMNAME

# Save all parameters to a parameter file
echo "Simulation box: $PARAM_LX x $PARAM_LY x $PARAM_LZ ($PARAM_NX x $PARAM_NY x $PARAM_NZ cells)" > parameters.$PARAM_CASE.dat
grep "PARAM_.*[a-Z]=" ../$0 | grep -v "echo" >> parameters.$PARAM_CASE.dat

echo "PARAM_MASS=$PARAM_MASS" >> parameters.$PARAM_CASE.dat
echo "PARAM_SIGMA=$PARAM_SIGMA" >> parameters.$PARAM_CASE.dat
echo "PARAM_LJCUTOFF=$PARAM_LJCUTOFF" >> parameters.$PARAM_CASE.dat
echo "PARAM_NBEADS=$PARAM_NBEADS" >> parameters.$PARAM_CASE.dat
echo "PARAM_OUTERRADIUS=$PARAM_OUTERRADIUS" >> parameters.$PARAM_CASE.dat
echo "PARAM_INNERRADIUS=$PARAM_INNERRADIUS" >> parameters.$PARAM_CASE.dat
echo "PARAM_ANGULARK=$PARAM_ANGULARK" >> parameters.$PARAM_CASE.dat
echo "PARAM_ETA=$PARAM_ETA" >> parameters.$PARAM_CASE.dat
echo "PARAM_NSTEPS=$PARAM_NSTEPS" >> parameters.$PARAM_CASE.dat
echo "PARAM_SAVEFREQ=$PARAM_SAVEFREQ" >> parameters.$PARAM_CASE.dat
echo "PARAM_SAMPLEFREQ=$PARAM_SAMPLEFREQ" >> parameters.$PARAM_CASE.dat
echo "PARAM_LIPOSOMETHICKNESS=$PARAM_LIPOSOMETHICKNESS" >> parameters.$PARAM_CASE.dat

# Save launch script to simulation directory
cp ../$0 .

# Write shear tensor to file
echo $SHEAR_TENSOR > shear_tensor.dat

# Create a wall file
echo "Creating wall file..."
../lattice -t hcp -f $PARAM_WALLPSEPARATION -Lx $PARAM_LX -Ly $PARAM_LY -Lz $PARAM_WALLWIDTH -r $PARAM_SIGMA | \
  grep -v "#" | \
  awk '{if(NF > 2) print $1, $2, $3}' > $PARAM_WALLFILE

WN=`cat $PARAM_WALLFILE | wc -l` # Number of particles in the wall file

echo "$WN particles per wall."

echo "Creating liposome file..."
../make_ball.sh $PARAM_INNERRADIUS $PARAM_OUTERRADIUS $PARAM_WALLPSEPARATION 0 0 0 255 | awk '{if(NF > 2) print $1, $2, $3}' > liposome.tmp

LIPOSOMEN=`cat liposome.tmp | wc -l` # Number of particles in the liposome file
echo "$LIPOSOMEN particles per liposome."

# Amplitude of the forcing
MDELTABLOB=`echo "$PARAM_RHO*$PARAM_DELTA*$PARAM_LX*$PARAM_LY/$WN" | bc -l`
OMEGA=`echo "pi = 4*a(1); 2*pi*$PARAM_WALLFREQ" | bc -l`
WALLF0=`echo "pi = 4*a(1); sqrt((($PARAM_MASS + $MDELTABLOB)*$OMEGA^2 - $PARAM_FBOND)^2 + $MDELTABLOB^2*$OMEGA^4)*$PARAM_WALLAMP" | bc -l` # Magnitude of wall forcing

echo "Writing data.main..."
# Create simulation configuration file
cat << EOF > $PARAM_PARAMETERFILE
#### GPU parameters ####
# Set GPU device
#setDevice 3

# Memory allocation limits
maxNumberPartInCellNonBonded    16
maxNumberPartInCell     	16

#### Numerical integration ####
# Choose scheme
quasiNeutrallyBuoyant
#stokesLimit

#Dimensions of the simulation box
celldimension  $PARAM_LX $PARAM_LY $PARAM_LZ

#Number of fluid cells in the directions x, y and z
cells   $PARAM_NX $PARAM_NY $PARAM_NZ

# Time steps
numstepsRelaxation	0
numsteps            $PARAM_NSTEPS
samplefreq          $PARAM_SAMPLEFREQ
savefreq            $PARAM_SAVEFREQ
dt			$PARAM_TIMESTEP

#### Fluid properties ####
initfluid 1
densfluid               $PARAM_RHO
shearviscosity          $PARAM_ETA
temperature             $PARAM_KT

saveFluid 0
saveVTK 0

#### Particle properties ####
loadparticles		1
particles		1
colors 	 	1

# Excess mass of particles
mass $PARAM_EXCESSMASS

# Interaction forces
computeNonBondedForces 0
cutoff                  $PARAM_LJCUTOFF

# Harmonic bonds
bondedForces               $PARAM_BONDFILE
threeBondedForces          $PARAM_ANGULARSPRINGSFILE

# Initial configuration file
coordinates		$PARAM_PARTICLEFILE

#Horizontal wall parameters   Last particle    K spring    z coordinate
wall1Params                    $WN              $PARAM_WALLBOND       $(($PARAM_LZ/2))
wall2Params                    $WN              $PARAM_WALLBOND       $(($PARAM_LZ/2))
wall3Params                    $((2*$WN))              $PARAM_WALLBOND       0.0

#Stokes flow parameters        Frequency  Amplitude  Duration
stokesFlowParams               $PARAM_WALLFREQ        $WALLF0        $PARAM_WALLFORCETMAX

outputname		./$PARAM_CASE
EOF

# Create run script
cat << EOF > $PARAM_RUNSCRIPT
  mkfifo $PARAM_CASE.particles
  awk 'BEGIN{np = 0; x = 0; y = 0; z = 0; t = $PARAM_FRAMESTEP;} \
       (\$5 == $PARAM_BOTTOMWALLC){x = x + \$1; y = y + \$2; z = z + \$3; np = np + 1} \
       (\$5 == $PARAM_LIPOSOMEC || \$5 == $PARAM_CHAINC){printf "%s\n", \$0 > "$PARAM_CASE.liposome"} \
       (NF == 1 && NR > 1){print t, x/np, y/np, z/np; np = 0; x = 0; y = 0; z = 0; t = t + $PARAM_FRAMESTEP; printf "\n" > "$PARAM_CASE.liposome" }' \
       < $PARAM_CASE.particles > wall.dat &
  ./fluam > /dev/null
  rm $PARAM_CASE.particles
EOF

# Get fluam executable program
cp ../fluam .

###### BEGIN create new configurations ######

if [ $PARAM_NEWCONFIGURATION == "Y" ]
then
# Total number of particles
echo $((2*$WN + $PARAM_NLIPOSOMES*($LIPOSOMEN + $PARAM_NBEADS))) > $PARAM_PARTICLEFILE

# Top wall
awk "{print \$1, \$2, \$3 + 0.5*$PARAM_LZ + 0.25*$PARAM_WALLPSEPARATION, $PARAM_TOPWALLC}" $PARAM_WALLFILE >> $PARAM_PARTICLEFILE

# Centre wall
awk "{print \$1, \$2, \$3 + 0.25*$PARAM_WALLPSEPARATION, $PARAM_BOTTOMWALLC}" $PARAM_WALLFILE >> $PARAM_PARTICLEFILE

# Elastic network
echo "Calculating bonds..."
../elastic_network $((2*$WN)) $PARAM_RCUT $PARAM_WALLBOND $PARAM_LX $PARAM_LY -1 $PARAM_PARTICLEFILE | \
  awk '{if(NF > 1) print $1, $2, $3, $4 "\n" $2, $1, $3, $4; else print $0}' > bonds.tmp

echo "Liposome bonds..."
../elastic_network $LIPOSOMEN $PARAM_RCUT $PARAM_LBOND -1 -1 -1 liposome.tmp | \
  awk '{if(NF > 1) print $1, $2, $3, $4 "\n" $2, $1, $3, $4; else print $0}' > liposomebonds.tmp

NBONDS=$((`cat bonds.tmp | wc -l` + $PARAM_NLIPOSOMES*$((`cat liposomebonds.tmp | wc -l` + 2*$PARAM_NBEADS + 2))))

# Number of bonds in bond files
echo $NBONDS > $PARAM_BONDFILE
echo $NANGULAR_SPRINGS > $PARAM_ANGULARSPRINGSFILE

# Strand beads positions
cat <<EOF > strand.bc
/* Read x and y coordinates of strand */
x = read(); y = read();
for(i = 0; i < $PARAM_NBEADS; i++)
  print x, "\t", y, "\t", $PARAM_SIGMA*(1.5 + i), "\t", $PARAM_CHAINC, "\n";
quit;
EOF

# Strand bonds
cat <<EOF > strand_bonds.bc
offset = read();
for(i = 0; i < $PARAM_NBEADS - 1; i++) {
  print offset + i, "\t", offset + i + 1, "\t", $PARAM_PBOND, "\t", $PARAM_SIGMA, "\n";
  print offset + i + 1, "\t", offset + i, "\t", $PARAM_PBOND, "\t", $PARAM_SIGMA, "\n";
}
quit;
EOF

# Angular bonds
cat <<EOF > angular_bonds.bc
offset = read();
for(i = 0; i < $PARAM_NBEADS - 2; i++) {
  print offset + i, "\t", offset + i + 1, "\t", offset + i + 2, "\t", $PARAM_ANGULARK, "\t", $PARAM_SIGMA, "\n";
}
quit;
EOF

# Find the southernmost particle
cat <<EOF > liposome_link.bc
linkx = linky = linkz = 0;
for(i = 0; i < $LIPOSOMEN; i++) {
  x = read(); y = read(); z = read()
  if(z < linkz) {
    linkx = x;
    linky = y;
    linkz = z;
    linkn = i;
  }
  if(z == linkz) {
    if(x^2 + y^2 < linkx^2 + linky^2) {
      linkx = x;
      linky = y;
      linkn = i;
    }
  }
}

print linkn;
quit;
EOF

LIPOSOME_LINK=`bc -l liposome_link.bc < liposome.tmp`

# Polymers and liposome
echo "Building liposomes and strands..."

# Create list of random positions for the liposomes
../no_overlap $PARAM_NLIPOSOMES $PARAM_OUTERRADIUS $PARAM_LX $PARAM_LY | awk '{print $1, $2, 1}' > random_positions.tmp

# Clear bondinfo.dat
rm -f bondinfo.dat

# Clear BD initial configurations
rm -f $PARAM_BDINITFILE

POLYMER_NUM=0
OFFSET=$((2*$WN))
# for PARTICLE_NUM in `seq 1 $WN | sort -R --random-source=/dev/urandom | head -n $PARAM_NLIPOSOMES`
for PARTICLE_NUM in `../find_closest_particle $PARAM_NLIPOSOMES random_positions.tmp $PARAM_WALLFILE`
do
    POLYMER_NUM=$(($POLYMER_NUM + 1))
    PARTICLE_XY=`awk 'NR == '$PARTICLE_NUM' {print $1, $2; exit}' $PARAM_WALLFILE` # Get xy coordinates
    PARTICLE_X=`echo $PARTICLE_XY | awk '{print $1}'` # x coordinate
    PARTICLE_Y=`echo $PARTICLE_XY | awk '{print $2}'` # y coordinate

    echo "$PARTICLE_X $PARTICLE_Y" | bc -l strand.bc >> $PARAM_PARTICLEFILE
    echo "$PARTICLE_X $PARTICLE_Y" | bc -l strand.bc | awk '{print $1, $2, $3, '$PARAM_RADIUS', '$PARAM_MASS', $4}' >> $PARAM_BDINITFILE

    echo "Wall-strand link at ($PARTICLE_X, $PARTICLE_Y): wall particle $(($PARTICLE_NUM + $WN - 1)) , strand particle $OFFSET" >> bondinfo.dat
    # Link strand to bottom wall
    echo "$(($PARTICLE_NUM + $WN - 1)) $OFFSET $PARAM_PBOND $PARAM_SIGMA" >> bonds.tmp
    echo "$OFFSET $((PARTICLE_NUM + $WN - 1)) $PARAM_PBOND $PARAM_SIGMA" >> bonds.tmp

    # Strand bonds
    echo $OFFSET | bc -l strand_bonds.bc >> bonds.tmp

    # Angular bonds
    echo $OFFSET | bc -l angular_bonds.bc >> $PARAM_ANGULARSPRINGSFILE

    # Increase offset
    OFFSET=$(($OFFSET + $PARAM_NBEADS))

    # Add liposome to particle file
    ZCOORD=`echo "$PARAM_SIGMA*($PARAM_NBEADS + 1) + $PARAM_OUTERRADIUS" | bc -l`
    awk '{if(NF > 1) print $1 + '$PARTICLE_X', $2 + '$PARTICLE_Y', $3 + '$ZCOORD', '$PARAM_LIPOSOMEC'}' liposome.tmp >> $PARAM_PARTICLEFILE

    # Add liposome to BD initial configuration file
    LIPOSOMER=`echo "$PARAM_OUTERRADIUS + $PARAM_RADIUS" | bc -l`
    echo -e "$PARTICLE_X\t$PARTICLE_Y\t$ZCOORD\t$LIPOSOMER\t$PARAM_MASS\t$PARAM_LIPOSOMEC" >> $PARAM_BDINITFILE

    # Link liposome to strand
#    echo "$(($OFFSET - 1)) $(($OFFSET + $LIPOSOME_LINK + 2)) $PARAM_LBOND $PARAM_SIGMA" >> bonds.tmp
#    echo "$(($OFFSET + $LIPOSOME_LINK + 2)) $(($OFFSET - 1)) $PARAM_LBOND $PARAM_SIGMA" >> bonds.tmp
    echo "$(($OFFSET - 1)) $(($OFFSET + $LIPOSOME_LINK)) $PARAM_LBOND $PARAM_SIGMA" >> bonds.tmp
    echo "$(($OFFSET + $LIPOSOME_LINK)) $(($OFFSET - 1)) $PARAM_LBOND $PARAM_SIGMA" >> bonds.tmp

#    echo "Strand-liposome link: strand particle $(($OFFSET - 1)), liposome particle $(($OFFSET + $LIPOSOME_LINK + 2))" >> bondinfo.dat
    echo "Strand-liposome link: strand particle $(($OFFSET - 1)), liposome particle $(($OFFSET + $LIPOSOME_LINK))" >> bondinfo.dat

    # Add liposome bonds to bondfile
    awk '{if(NF > 1) print '$OFFSET' + $1, '$OFFSET' + $2, $3, $4}' liposomebonds.tmp >> bonds.tmp

    # Increase offset
    OFFSET=$(($OFFSET + $LIPOSOMEN))
done

# Sort bonds for fluam
sort -k1 -k2 -n bonds.tmp >> $PARAM_BONDFILE

# Add fixed bonds to bondfile
if [ $PARAM_FBOND != "0" ]
then
  echo $((2*$WN)) >> $PARAM_BONDFILE
  awk '(NR > 1 && NR <= 2*'$WN' + 1){print NR - 2, '$PARAM_FBOND', 0, $1, $2, $3}' $PARAM_PARTICLEFILE >> $PARAM_BONDFILE
fi

# Create files for nanometre Brownian dynamics simulation
if [ $PARAM_RUNBD == "Y"  ]
then
# Initial configuration file for nanometre
# awk '(NR > 1 && $4 == '$PARAM_CHAINC'){print $1, $2, $3, '$PARAM_RADIUS', '$PARAM_MASS', $4}' $PARAM_PARTICLEFILE > $PARAM_BDINITFILE
# cat liposome_centres.tmp >> $PARAM_BDINITFILE

# Fixed bonds
awk '((NR % ('$PARAM_NBEADS' + 1)) == 1){print NR - 1, $1, $2, 0.25, '$PARAM_PBOND', '$PARAM_SIGMA'}' $PARAM_BDINITFILE| head -n $PARAM_NLIPOSOMES > $PARAM_BDBONDFILE

# Harmonic bonds
BDNHARMONIC=$(($PARAM_NLIPOSOMES*$PARAM_NBEADS))
OFFSET=0
for i in `seq 1 $PARAM_NLIPOSOMES`
do
  echo $OFFSET | bc -l strand_bonds.bc | awk '(NR % 2 == 1){print}' >> $PARAM_BDBONDFILE
  VSKIP=`echo $PARAM_OUTERRADIUS + $PARAM_RADIUS | bc -l`
#  echo -e "$(($i*$PARAM_NBEADS - 1))\t$(($PARAM_NLIPOSOMES*$PARAM_NBEADS + $i - 1)) $PARAM_PBOND $VSKIP" >> $PARAM_BDBONDFILE
  echo -e "$(($i*($PARAM_NBEADS + 1) - 2))\t$(($i*($PARAM_NBEADS + 1) - 1)) $PARAM_PBOND $VSKIP" >> $PARAM_BDBONDFILE
  OFFSET=$(($OFFSET + $PARAM_NBEADS + 1))
done

# Add angular springs to bond file
# awk '(NR>1){print $1 - 2*'$WN', $2 - 2*'$WN', $3 - 2*'$WN', $4, 0}' $PARAM_ANGULARSPRINGSFILE >> $PARAM_BDBONDFILE
OFFSET=0
for i in `seq 1 $PARAM_NLIPOSOMES`
do
  echo $OFFSET | bc -l angular_bonds.bc | awk '{print $1, $2, $3, $4, 0}' >> $PARAM_BDBONDFILE
  OFFSET=$(($OFFSET + $PARAM_NBEADS))
done

### Parameters for the Brownian dynamics initial configuration simulation ###
BDNPARTICLES=$(($PARAM_NLIPOSOMES*($PARAM_NBEADS + 1))) # Number of particles
BDLJCUTOFF=$PARAM_LJCUTOFF
BDSAMPLING=`echo "scale = 4; ($PARAM_BDTMAX - $PARAM_BDTRELAX)/$PARAM_BDSAMPLES" | bc -l` # Sampling frequency
BDCELLSIZE=`echo "0.25*$PARAM_LX" | bc -l` # Cell size for linked lists
BDNLRADIUS=`echo "0.5*$BDCELLSIZE" | bc -l` # Neighbour list radius
BDGAMMA=`echo "scale = 6; pi = 4*a(1); 6*pi*$PARAM_ETA*$PARAM_RADIUS" | bc -l`

echo "Writing Brownian dynamics nanometre.parameters.dat..."
# Create initial configurations parameter file
cat << EOF > nanometre.parameters.dat
# Free comment
comment Brownian dynamics simulation

# Number of particles
np $BDNPARTICLES

# System dimensionality
dim 3

### Interactions ###

# Lennard-Jones
# force_cutoff 2.5 # Standard Lennard-Jones cutoff (in units of sigma)
force_cutoff $BDLJCUTOFF # Weeks-Chandler-Anderson potential
epsilon $PARAM_LJEPSILON # Lennard-Jones epsilon

### Time ###
t0 0              # Initial time
trelax $PARAM_BDTRELAX    # Relaxation time
tmax $PARAM_BDTMAX        # Final time
dt $PARAM_BDDT            # Integration time step
sampling_freq $BDSAMPLING # Sampling frequency

### Initial conditions ###
initial_conditions $PARAM_BDINITFILE

### Shear tensor ###
shear_tensor shear_tensor.dat

### Space ###
L 4 4 8 # Simulation box size in cells
cellsize $BDCELLSIZE
neighbour_list_radius $BDNLRADIUS

### Bonds ###
bonds_filename $PARAM_BDBONDFILE
nfixed $PARAM_NLIPOSOMES
nbonds $BDNHARMONIC
nangular_springs $NANGULAR_SPRINGS

### Langevin thermostat ###
kT $PARAM_BDKT # Thermal energy
gamma $BDGAMMA # Friction coefficient
EOF

cp ../Brownian_dynamics .
fi

fi ###### END of create new configurations ######

cat << EOF > launch.sh
#!/bin/sh

for dir in \$(ls -d ${PARAM_CASE}-*)
do
  cd \$dir
  echo "Launching simulation \$dir..."
  ssubgpu $PARAM_NODE \$dir $PARAM_RUNSCRIPT \$dir.out
  cd ..
done
EOF

chmod a+x launch.sh

cat << EOF > rotate_liposome.bc
  scale = 5;

  /* Read x, y and z position of the liposome */
  liposome_r[0] = read();
  liposome_r[1] = read();
  liposome_r[2] = read();

  /* Read x, y and z position of the strand end */
  strand_r[0] = read();
  strand_r[1] = read();
  strand_r[2] = read();

  /* Read number of particles in the liposome file */
  liposomen = read();
  /* Read positions of the particles in the liposome */
  for(i = 0; i < liposomen; i++) {
    q[3*i]     = read();
    q[3*i + 1] = read();
    q[3*i + 2] = read();
  }

  /* Distance from liposome centre to strand */
  for(i = 0; i < 3; i++) dist2 += (strand_r[i] - liposome_r[i])^2;
  dist = sqrt(dist2);

  /* Calculate the cosine of the angle of rotation */
  costheta = (liposome_r[2] - strand_r[2])/dist;
  sintheta = sqrt(1 - costheta^2);

  /* Calculate the axis of rotation */
  u[0] = strand_r[1] - liposome_r[1];
  u[1] = liposome_r[0] - strand_r[0];
  u[2] = 0;

  /* Norm of the axis of rotation */
  for(i = 0; i < 3; i++) norm2 += u[i]^2;
  norm = sqrt(norm2);

  if(norm > 0) {
    /* Normalise axis vector */
    for(i = 0; i < 3; i++) u[i] /= norm;
  }

  /* Rotated orientation vector */
  newdir[0] = -u[1]*sintheta;
  newdir[1] = u[0]*sintheta;
  newdir[2] = -costheta;

  /* Orientation check (newdir should be parallel to strand_r - liposome_r) */
  for(i = 0; i < 3; i++) orientation_check += newdir[i]*(strand_r[i] - liposome_r[i]);

  /* If orientation is wrong, flip the sign of sintheta */
  if(orientation_check < 0.98*dist) sintheta *= -1;

  /* Calculate the rotation matrix (for u[2] = 0) */
  if(norm > 0) {
    r[0] = costheta + u[0]^2*(1 - costheta); r[1] = u[0]*u[1]*(1 - costheta);         r[2] = u[1]*sintheta;
    r[3] = r[1];                             r[4] = costheta + u[1]^2*(1 - costheta); r[5] = -u[0]*sintheta;
    r[6] = -r[2];                            r[7] = -r[5];                            r[8] = costheta;
  } else {
    r[0] = r[4] = r[8] = 1;
  }

  /* Rotate positions */
  for(i = 0; i < liposomen; i++) {
    newq[0] = newq[1] = newq[2] = 0;
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        newq[j] += r[3*j + k]*q[3*i + k];
      }
    }

    /* Output rotated position */
    print newq[0] + liposome_r[0], "\t", newq[1] + liposome_r[1], "\t", newq[2] + liposome_r[2], "\t", $PARAM_LIPOSOMEC, "\n";
  }

  quit;
EOF

cat << EOF > $PARAM_BDRUNSCRIPT
if [ $PARAM_RUNBD == "Y" ]
then
export OMP_NUM_THREADS=1
./Brownian_dynamics nanometre.parameters.dat > $PARAM_BDCONFIGURATIONS 2> /dev/null
fi
rm -Rf ${PARAM_CASE}-*

if [ $PARAM_NLIPOSOMES -lt 1 ]
then
  mkdir ${PARAM_CASE}-0
  head -n $((2*$WN + 1)) $PARAM_PARTICLEFILE > ${PARAM_CASE}-0/$PARAM_PARTICLEFILE
  cp $PARAM_PARAMETERFILE ${PARAM_CASE}-0/data.main
  cp $PARAM_BONDFILE ${PARAM_CASE}-0
  cp $PARAM_ANGULARSPRINGSFILE ${PARAM_CASE}-0
  cp $PARAM_RUNSCRIPT ${PARAM_CASE}-0
  cp fluam ${PARAM_CASE}-0
  # Ignore velocity files
  rm -f $PARAM_SIMNAME.velocityParticles*
  ln -s /dev/null ${PARAM_CASE}-0/$PARAM_CASE.velocityParticles
  ln -s /dev/null ${PARAM_CASE}-0/$PARAM_CASE.velocityParticlesI
else
i=0
grep -v "#" $PARAM_BDCONFIGURATIONS | while read x y z r col
do
  # Get the initial x and y coordinates at the bottom of the first DNA strand
  STRANDX=\$(head -n 1 $PARAM_BDINITFILE | awk '{print \$1}')
  STRANDY=\$(head -n 1 $PARAM_BDINITFILE | awk '{print \$2}')

  # Get the index of the wall particle connected to the DNA strand
  WALLPARTICLEN=\$(awk '{print \$8}' bondinfo.dat | head -n 1)

  # Get the initial x and y coordinates of the wall particle connected to the DNA strand
  WALLPARTICLEX=\$(awk '(NR == '\$WALLPARTICLEN' - 2){print \$1; exit}' $PARAM_PARTICLEFILE)
  WALLPARTICLEY=\$(awk '(NR == '\$WALLPARTICLEN' - 2){print \$2; exit}' $PARAM_PARTICLEFILE)

  if [ "\$x" == "" ]; then
    i=\$((\$i + 1))
  else
    if [ ! -d ${PARAM_CASE}-\$i ]; then
      mkdir ${PARAM_CASE}-\$i
      head -n $((2*$WN + 1)) $PARAM_PARTICLEFILE > ${PARAM_CASE}-\$i/$PARAM_PARTICLEFILE
      cp $PARAM_PARAMETERFILE ${PARAM_CASE}-\$i/data.main
      cp $PARAM_BONDFILE ${PARAM_CASE}-\$i
      cp $PARAM_ANGULARSPRINGSFILE ${PARAM_CASE}-\$i
      cp $PARAM_RUNSCRIPT ${PARAM_CASE}-\$i
      cp fluam ${PARAM_CASE}-\$i
      # Ignore velocity files
      rm -f $PARAM_SIMNAME.velocityParticles*
      ln -s /dev/null ${PARAM_CASE}-\$i/$PARAM_CASE.velocityParticles
      ln -s /dev/null ${PARAM_CASE}-\$i/$PARAM_CASE.velocityParticlesI
    fi
    if [ "\$col" == "$PARAM_LIPOSOMEC" ]; then
      echo -e "\$x\t\$y\t\$z" > liposome_info.tmp
      echo -e "\$PREVX\t\$PREVY\t\$PREVZ" >> liposome_info.tmp
      echo $LIPOSOMEN >> liposome_info.tmp
       cat liposome_info.tmp liposome.tmp | bc -l rotate_liposome.bc >>  ${PARAM_CASE}-\$i/$PARAM_PARTICLEFILE
    elif [ "\$col" == "$PARAM_CHAINC" ]; then
#      echo -e "\$x\t\$y\t\$z\t\$col" | awk '{print \$1 - '\$WALLPARTICLEX', \$2 - '\$WALLPARTICLEY', \$3, \$4}' >> ${PARAM_CASE}-\$i/$PARAM_PARTICLEFILE
      echo -e "\$x\t\$y\t\$z\t\$col" >> ${PARAM_CASE}-\$i/$PARAM_PARTICLEFILE
      PREVX=\$x
      PREVY=\$y
      PREVZ=\$z
    else
      echo -e "\$x\t\$y\t\$z\t\$col" >> ${PARAM_CASE}-\$i/$PARAM_PARTICLEFILE
    fi
  fi
done
fi

rm -f liposome.tmp
rm -f liposome_info.tmp
rm -f rotate_liposome.bc
EOF

# Submit job
if [ $PARAM_RUNBD == "Y" ]
then
echo "Submitting job $PARAM_SIMNAME"
ssubcpu $PARAM_BDNODE $PARAM_SIMNAME $PARAM_BDRUNSCRIPT $PARAM_SIMNAME.out
else
echo "Creating $PARAM_CASE directories..."
bash ./$PARAM_BDRUNSCRIPT
fi

rm -f angular_bonds.bc
rm -f bonds.tmp
rm -f liposomebonds.tmp
rm -f liposome_link.bc
rm -f strand.bc
rm -f strand_bonds.bc
rm -f wall.tmp
rm -r random_positions.tmp
