// Filename: main.cpp
//
// Copyright (c) 2010-2017, Florencio Balboa Usabiaga
//
// This file is part of Fluam
//
// Fluam is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fluam is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Fluam. If not, see <http://www.gnu.org/licenses/>.


#define GLOBAL 1
#define TEMPERATURE 1

#include "header.h"
#include "headerOtherFluidVariables.h"
#include "temperature.h"
#include "binary_src_backup.h"

int main(int argc, char* argv[]){

  std::cerr<<"COMPILED AT "<<__TIME__<<" on "<<__DATE__<<endl;
  for(int i = 1; i<argc; i++){
    if(std::string(argv[i]).compare("--dump-src") == 0){
      std::ofstream srcout("src.tar.xz", std::ios::binary);
      std::cerr<<"DUMPING ORIGINAL SOURCE AND EXITING"<<endl;
      srcout.write((const char*)src_tar_xz, src_tar_xz_len);
      exit(0);
    }
  }
    
  
  //Read a file with data for start the simulation
  if(!loadDataMain(argc,argv)) return 0;
  
  //Save GPU information
  if(!cudaDevice()) return 0;
  
  //Initialize random number generator
  if(!initializeRandomNumber()) return 0;

  //Choose scheme
  if(setContinuousGradient==1){
    if(!schemeContinuousGradient()) return 0;
  }
  else if(setGiantFluctuations==1){
    if(!schemeGiantFluctuations()) return 0;
  }
  else if(setBinaryMixtureWall==1){
    if(!schemeBinaryMixtureWall()) return 0;
  }
  else if(setBinaryMixture==1){
    if(!schemeBinaryMixture()) return 0;
  }
  else if(particlesWall){
    if(!schemeParticlesWall()) return 0;
  }
  else if((setboundary || setparticles) && (thermostat)){
    if(!schemeCompressibleParticles()) return 0;
  }
  else if(freeEnergyCompressibleParticles){
    if(!schemeFreeEnergyCompressibleParticles()) return 0;
  }
  else if(semiImplicitCompressibleParticles){
    if(!schemeSemiImplicitCompressibleParticles()) return 0;
  }
  else if(momentumCoupling){
    if(!schemeMomentumCoupling()) return 0;
  }
  else if(thermostat){
    if(!schemeThermostat()) return 0;
  }
  else if((setboundary || setparticles) && incompressible){
    if(!schemeIncompressibleBoundary()) return 0;
  }
  else if(incompressible){
    if(!schemeIncompressible()) return 0;
  }
  else if(incompressibleBinaryMixture){
    if(!schemeIncompressibleBinaryMixture()) return 0;
  }
  else if(incompressibleBinaryMixtureMidPoint){
    if(!schemeIncompressibleBinaryMixtureMidPoint()) return 0;
  }
  else if(quasiNeutrallyBuoyant){
    if(!schemeQuasiNeutrallyBuoyant()) return 0;
  }
  else if(quasiNeutrallyBuoyant2D){
    if(!schemeQuasiNeutrallyBuoyant2D()) return 0;
  }
  else if(quasiNeutrallyBuoyant4pt2D){
    if(!schemeQuasiNeutrallyBuoyant4pt2D()) return 0;
  }
  else if(quasi2D or stokesLimit2D){
    if(!schemeQuasi2D()) return 0;
  }
  else if(stokesLimitFirstOrder==1){
    if(!schemeStokesLimitFirstOrder()) return 0;
  }
  else if(stokesLimit==1 && bigSystem==1){
    if(!schemeStokesLimitBigSystem()) return 0;
  }
  else if(stokesLimit==1){
    if(!schemeStokesLimit()) return 0;
  }
  else if(setboundary || setparticles){
    if(!schemeBoundary()) return 0;
  }
  else if(setGhost==1){
    if(!schemeRK3Ghost()) return 0;
  }
  else{
    if(!schemeRK3()) return 0;
  }

  cout << "END" << endl;
  return 0;
}
