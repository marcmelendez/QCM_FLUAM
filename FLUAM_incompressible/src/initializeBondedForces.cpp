// Filename: initializeBondedForces.cu
//
// Copyright (c) 2010-2016, Florencio Balboa Usabiaga
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



#include <cstring>
#include <math.h>
#include "header.h"
#include "particles.h"
#include "fluid.h"
#include "parameters.h"
#include "cells.h"
#include<vector>


#define fori(x,y) for(int i=x; i<y;i++)
#define forj(x,y) for(int j=x; j<y;j++)



bool initializeBondedForces(){
  
  int index1, index2;
  int trashInt;
  double trashDouble;

  if(bondedForcesVersion){
    initializeBondedForcesOldVersion();
    return 1;
  }
  

  //OPEN FILE
  ifstream file(bondedForcesFile.c_str());

  //Number of particle-particle bonds
  file >> nbondsParticleParticle;
  
  //Allocate memory
  bondsParticleParticle = new int [np];
  bondsParticleParticleOffset = new int [np];

  //Initially no particle has bonds
  for(int i=0;i<np;i++)
    bondsParticleParticle[i] = 0;

  //Information bonds particle-particle
  for(int i=0;i<nbondsParticleParticle;i++){  
    file >> index1 >> index2 >> trashDouble >> trashDouble;   
    bondsParticleParticle[index1]++;
    bondsParticleParticle[index2]++;
  }
  




  //Number of particle-fixedPoints bonds
  file >> nbondsParticleFixedPoint;

  //Allocate memory
  bondsParticleFixedPoint = new int [np];
  bondsParticleFixedPointOffset = new int [np];

  //Initially no particle has bonds
  for(int i=0;i<np;i++)
    bondsParticleFixedPoint[i] = 0;

  //Information bonds particle-fixedPoint
  for(int i=0;i<nbondsParticleFixedPoint;i++){  
    file >> index1 >> trashDouble >> trashDouble >> trashDouble >> trashDouble >> trashDouble;
    bondsParticleFixedPoint[index1]++;
  }
  
  //Important, lear how to rewind a file
  //CLOSE FILE 
  file.close();
  
  //Allocate memory
  bondsIndexParticleParticle = new int [nbondsParticleParticle * 2];
  kSpringParticleParticle = new double [nbondsParticleParticle * 2];
  r0ParticleParticle = new double [nbondsParticleParticle * 2];

  //bondsIndexParticleFixedPoint = new int [nbondsParticleFixedPoint];  
  kSpringParticleFixedPoint = new double [nbondsParticleFixedPoint];
  r0ParticleFixedPoint = new double [nbondsParticleFixedPoint];
  rxFixedPoint = new double [nbondsParticleFixedPoint];
  ryFixedPoint = new double [nbondsParticleFixedPoint];
  rzFixedPoint = new double [nbondsParticleFixedPoint];

  //Compute offset
  bondsParticleParticleOffset[0] = 0;
  bondsParticleFixedPointOffset[0] = 0;
  for(int i=1; i<np; i++){
    bondsParticleParticleOffset[i] = bondsParticleParticleOffset[i-1] + bondsParticleParticle[i-1];
    bondsParticleFixedPointOffset[i] = bondsParticleFixedPointOffset[i-1] + bondsParticleFixedPoint[i-1];
  }
  
  //Create tmp offset
  int *tmpOffset = new int [np];
  for(int i=0;i<np;i++){
    tmpOffset[i] = 0;
  }

  //OPEN THE FILE AGAIN
  file.open(bondedForcesFile.c_str());

  //Number of particle-particle bonds
  file >> nbondsParticleParticle;

  //Information bonds particle-particle
  double k, r0;
  for(int i=0;i<nbondsParticleParticle;i++){  
    file >> index1 >> index2 >> k >> r0;
    
    // Data for particle index1
    bondsIndexParticleParticle[bondsParticleParticleOffset[index1]+tmpOffset[index1]] = index2;
    kSpringParticleParticle[   bondsParticleParticleOffset[index1]+tmpOffset[index1]] = k;
    r0ParticleParticle[        bondsParticleParticleOffset[index1]+tmpOffset[index1]] = r0;

    // Data for particle index2
    bondsIndexParticleParticle[bondsParticleParticleOffset[index2]+tmpOffset[index2]] = index1;
    kSpringParticleParticle[   bondsParticleParticleOffset[index2]+tmpOffset[index2]] = k;
    r0ParticleParticle[        bondsParticleParticleOffset[index2]+tmpOffset[index2]] = r0;

    // Increase tmpOffset
    tmpOffset[index1]++;
    tmpOffset[index2]++;
  }

  // Reset tmp offset to zero for particle-fixed point interactions
  for(int i=0;i<np;i++){
    tmpOffset[i] = 0;
  }

  //Number of particle-fixedPoints bonds
  file >> nbondsParticleFixedPoint;

  //Information bonds particle-fixedPoint
  for(int i=0;i<nbondsParticleFixedPoint;i++){  
    file >> index1 >> kSpringParticleFixedPoint[bondsParticleFixedPointOffset[index1]+tmpOffset[index1]]
	 >> r0ParticleFixedPoint[               bondsParticleFixedPointOffset[index1]+tmpOffset[index1]]
	 >> rxFixedPoint[                       bondsParticleFixedPointOffset[index1]+tmpOffset[index1]]
	 >> ryFixedPoint[                       bondsParticleFixedPointOffset[index1]+tmpOffset[index1]]
	 >> rzFixedPoint[                       bondsParticleFixedPointOffset[index1]+tmpOffset[index1]];

    // Increase tmpOffset
    tmpOffset[index1]++;
  }

  //CLOSE FILE
  file.close();

  // Free tmpOffset
  delete[] tmpOffset;
  
  cout << "INITALIZE BONDED FORCES :       DONE " << endl;

  return 1;
}

  //!*R THREE BONDED INIT STARTS
bool initializeThreeBondedForces(){
  ifstream file(threeBondedForcesFile.c_str());  
  //ONLY WRITE EACH BOND ONCE!
  /*
          O 4
          |
    0     |     2     3
    O-----O-----O-----O
         1|
          |
	  O 5


	  BondList.dat would be:
	  3
	  0 1 2 k r0
	  1 2 3 k r0
	  4 1 5 k r0
	  
    The order of the bonds in the file does not matter
   */
  
  uint nbonds;
  file>>nbonds;
  NbondsThreeParticle = nbonds;
  threebondList = new int[3*nbonds]; //Three particles per bond, just one entry per spring
  threeNbonds = new int[np]; //Number of bonds in which each particle is involved
  threeCumulativeIndex = new int[np]; //Sum of threeNbonds up to i-1
  vector<vector<int> > isinbondshelper(np); //vector equivalent of threeisinbonds, for dynamic alloc
  //Is in matrix form for convinience, later we place it as a 1D array in threeisinbonds

  threekSprings = new double[nbonds];
  threer0Springs = new double[nbonds];

  
  for(int i=0; i<np; i++) threeNbonds[i] = 0; //Initialize bond counter

  for(int i=0; i<nbonds; i++){
    for(int j=0; j<3;j++){
      int particle;
      file>>particle;
      threebondList[3*i + j]= particle;
      threeNbonds[ particle ]++; //Sum the bondcounter for the particle in this spring
      isinbondshelper[particle].push_back(i); //particle is involved in spring i 
    }
    file>> threekSprings[i] >> threer0Springs[i];
  }

  //Now initialize the list of bonds index for each particle from the helper vector :
  int Nbondsi;

  int isinbondsDIM = 0;
  fori(0, np){
    threeCumulativeIndex[i] = isinbondsDIM;
    isinbondsDIM += threeNbonds[i];
  }
  threeisinbonds = new int[isinbondsDIM];


  for(int i=0;i<np;i++){//Transform 2D array in 1D array
    Nbondsi = isinbondshelper[i].size();
    if(Nbondsi!= threeNbonds[i]) cout<<"ERROR in create bonded forces!!"<<endl;
    forj(0,Nbondsi) threeisinbonds[threeCumulativeIndex[i]+j] = isinbondshelper[i][j];
  }


  /*
  fori(0, nbonds*3) cout<<threebondList[i]<<" ";
  cout<<endl;
  
  fori(0, np) cout<<threeNbonds[i]<<" ";
  cout<<endl;

  fori(0, isinbondsDIM) cout<<threeisinbonds[i]<<" ";
  cout<<endl;
  
  fori(0, np) cout<<threeCumulativeIndex[i]<<" ";
  cout<<endl;
  */
 

  //Now all CPU info is stored
  cout<<"THREE BONDED SPRINGS DONE!!!"<<endl;
  return true;
}














bool freeBondedForces(){

  delete[] bondsParticleParticle;
  delete[] bondsParticleFixedPoint;

  delete[] bondsParticleParticleOffset;
  delete[] bondsParticleFixedPointOffset;


  delete[] bondsIndexParticleParticle;
  delete[] kSpringParticleParticle;
  delete[] r0ParticleParticle;
  

  //delete[] bondsIndexParticleFixedPoint;
  delete[] kSpringParticleFixedPoint;
  delete[] r0ParticleFixedPoint;
  delete[] rxFixedPoint;
  delete[] ryFixedPoint;
  delete[] rzFixedPoint;
  
}
