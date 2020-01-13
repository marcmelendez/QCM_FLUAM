// Filename: forceBondedGPU.cu
//
// Copyright (c) 2010-2015, Florencio Balboa Usabiaga
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

__device__ void forceBondedThreeParticleGPU(const int i,
					    double& fx, //Pass by reference
					    double& fy,
					    double& fz,
					    const double rx,
					    const double ry,
					    const double rz,
					    const threeParticleBondsVariables* tPBV){	


  double ri[3], rj[3], rl[3];
  double rij[3], ril[3];
  double rij2, ril2, rij1, ril1;
  double a1, a2, a3;
  double r, r0;
  double kSpring;
  double ampli;
  int p1, p2, p3;

  //Particle-Particle Force
  int nBonds = tPBV->Nbonds[i];
  int offset = tPBV->cumulative_index[i];
  int bond;
  for(int j=0;j<nBonds;j++){
    bond = tPBV->isinbonds[offset+j];
    
    p1 = tPBV->bondList[3*bond];
    p2 = tPBV->bondList[3*bond+1];
    p3 = tPBV->bondList[3*bond+2];

    //Particle bonded coordinates
    rj[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p1);
    rj[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p1);
    rj[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p1);
    //Particle bonded coordinates
    ri[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p2);
    ri[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p2);
    ri[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p2);
    //Particle bonded coordinates
    rl[0] = fetch_double(texrxboundaryGPU,nboundaryGPU+p3);
    rl[1] = fetch_double(texryboundaryGPU,nboundaryGPU+p3);
    rl[2] = fetch_double(texrzboundaryGPU,nboundaryGPU+p3);
    
    /*
    printf("i= %d     current bond: %d %d %d\n   pos\n p1: %.3f %.3f %.3f\n p2: %.3f %.3f %.3f\n p3: %.3f %.3f %.3f\n     ",
	   i, p1, p2, p3, rj[0], rj[1], rj[2], ri[0], ri[1], ri[2], rl[0], rl[1], rl[2]);
    */
    rij[0] = ri[0]-rj[0];    rij[1] = ri[1]-rj[1];    rij[2] = ri[2]-rj[2];
    rij2 = rij[0]*rij[0] + rij[1]*rij[1] + rij[2]*rij[2];
    rij1 = sqrt(rij2);

    ril[0] = ri[0]-rl[0];    ril[1] = ri[1]-rl[1];    ril[2] = ri[2]-rl[2];
    ril2 = ril[0]*ril[0] + ril[1]*ril[1] + ril[2]*ril[2];
    ril1 = sqrt(ril2);

    a1 = rij[0]*ril[0] + rij[1]*ril[1] + rij[2]*ril[2];
    a2 = rij1*ril1;
    a3 = a1/a2;             //a3 = cos (teta) = rij*ril / mod(rij)*mod(ril)                           


    //Spring constant
    kSpring = tPBV->kSprings[bond];
    //Equilibrium distance 
    r0 = tPBV->r0Springs[bond];  
    /*    
    if(a3<=-0.9999){
      ampli = -kSpring*sqrt(2.0/(1-a3));
    }
    */
     ampli = kSpring * (acos(a3)-3.1415)/sqrt(1-a3*a3);
    

    //p1 is j, p2 is i, p3 is l
    
    if(i==p1){
      fx += ampli * (a3*rij[0]/rij2 - ril[0]/a2);
      fy += ampli * (a3*rij[1]/rij2 - ril[1]/a2);
      fz += ampli * (a3*rij[2]/rij2 - ril[2]/a2);
      
      /* Harmonic bond (deactivated)
      r = rij1;
      fx += -kSpring * (1 - r0/r) * (rj[0] - ri[0]);
      fy += -kSpring * (1 - r0/r) * (rj[1] - ri[1]);
      fz += -kSpring * (1 - r0/r) * (rj[2] - ri[2]);
      */
    }
    else if(i==p2){
      
      fx += ampli * (-a3*(rij[0]/rij2 + ril[0]/ril2) + (1/a2)*(rij[0] + ril[0]));
      fy += ampli * (-a3*(rij[1]/rij2 + ril[1]/ril2) + (1/a2)*(rij[1] + ril[1]));
      fz += ampli * (-a3*(rij[2]/rij2 + ril[2]/ril2) + (1/a2)*(rij[2] + ril[2]));
      
      /* Harmonic bonds (deactivated)
      //First spring
      r = rij1;

      fx += -kSpring * (1 - r0/r) * (ri[0] - rj[0]);
      fy += -kSpring * (1 - r0/r) * (ri[1] - rj[1]);
      fz += -kSpring * (1 - r0/r) * (ri[2] - rj[2]);
      
      //Second spring
      r = ril1;

      fx += -kSpring * (1 - r0/r) * (ri[0] - rl[0]);
      fy += -kSpring * (1 - r0/r) * (ri[1] - rl[1]);
      fz += -kSpring * (1 - r0/r) * (ri[2] - rl[2]);
      */
    }
    else if(i==p3){
      fx += ampli * (a3*ril[0]/ril2 - rij[0]/a2);
      fy += ampli * (a3*ril[1]/ril2 - rij[1]/a2);
      fz += ampli * (a3*ril[2]/ril2 - rij[2]/a2);
      
      /* Harmonic bond (deactivated)
      r = ril1;

      fx += -kSpring * (1 - r0/r) * (rl[0] - ri[0]);
      fy += -kSpring * (1 - r0/r) * (rl[1] - ri[1]);
      fz += -kSpring * (1 - r0/r) * (rl[2] - ri[2]);
      */
    }

    //printf("i=%d\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",i, a1, a2, a3, ampli, fx, fy, fz);
    
  }
}




__device__ void forceBondedParticleParticleGPU(const int i,
					       double& fx, //Pass by reference
					       double& fy,
					       double& fz,
					       const double rx,
					       const double ry,
					       const double rz,
					       const bondedForcesVariables* bFV){	

  double x, y, z;
  double r, r0;
  double kSpring;
  int index;


  //Particle-Particle Force
  int nBonds = bFV->bondsParticleParticleGPU[i];
  int offset = bFV->bondsParticleParticleOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){

    index = bFV->bondsIndexParticleParticleGPU[offset + j];

    //if(i==0) index=1;
    //if(i==1) index=0;


    //Particle bonded coordinates
    x = fetch_double(texrxboundaryGPU,nboundaryGPU+index);
    y = fetch_double(texryboundaryGPU,nboundaryGPU+index);
    z = fetch_double(texrzboundaryGPU,nboundaryGPU+index);
    
    //Equilibrium distance 
    r0 = bFV->r0ParticleParticleGPU[offset+j];
    
    //Spring constant
    kSpring = bFV->kSpringParticleParticleGPU[offset+j];

    double rxx = rx-x;
    rxx -= floorf(rxx/lxGPU+0.5)*lxGPU;
    double ryy = ry-y;
    ryy -= floorf(ryy/lyGPU+0.5)*lyGPU;
    double rzz = rz-z;
    rzz -= floorf(rzz/lzGPU+0.5)*lzGPU;
    
    if(r0==0){
      fx += -kSpring * rxx;
      fy += -kSpring * ryy;
      fz += -kSpring * rzz;
    }
    else{     //If r0!=0 calculate particle particle distance
      r = sqrt(rxx*rxx+ryy*ryy+rzz*rzz);      
      if(r>0){//If r=0 -> f=0
	fx += -kSpring * (1 - r0/r) * rxx;
	fy += -kSpring * (1 - r0/r) * ryy;
	fz += -kSpring * (1 - r0/r) * rzz;
      }
    }
  }











  //Particle-FixedPoint Force
  nBonds = bFV->bondsParticleFixedPointGPU[i];
  offset = bFV->bondsParticleFixedPointOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){
    

    //Fixed point coordinates
    x = bFV->rxFixedPointGPU[offset+j];
    y = bFV->ryFixedPointGPU[offset+j];
    z = bFV->rzFixedPointGPU[offset+j];
    
    //Equilibrium distance 
    r0 = bFV->r0ParticleFixedPointGPU[offset+j];

    //Spring constant
    kSpring = bFV->kSpringParticleFixedPointGPU[offset+j];
    
    if(r0==0){
      fx += -kSpring * (rx - x);
      fy += -kSpring * (ry - y);
      fz += -kSpring * (rz - z);
    }  
    else{     //If r0!=0 calculate particle particle distance
      r = sqrt( (x-rx)*(x-rx) + (y-ry)*(y-ry) + (z-rz)*(z-rz) );
      if(r>0){//If r=0 -> f=0
	fx += -kSpring * (1 - r0/r) * (rx - x);
	fy += -kSpring * (1 - r0/r) * (ry - y);
	fz += -kSpring * (1 - r0/r) * (rz - z);
      }
    }
  }

    




  return ;
}


__device__ void forceBondedParticleParticleGPU_changed(const int i,
					       double& fx, //Pass by reference
					       double& fy,
					       double& fz,
					       double& Fatt,
                                               const double rx,
					       const double ry,
					       const double rz,
					       const bondedForcesVariables* bFV,
					       double xWall){	

  double x, y, z;
  double r, r0;
  double kSpring;
  int index;


  //Particle-Particle Force
  int nBonds = bFV->bondsParticleParticleGPU[i];
  int offset = bFV->bondsParticleParticleOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){

    index = bFV->bondsIndexParticleParticleGPU[offset + j];

    //if(i==0) index=1;
    //if(i==1) index=0;


    //Particle bonded coordinates
    x = fetch_double(texrxboundaryGPU,nboundaryGPU+index);
    y = fetch_double(texryboundaryGPU,nboundaryGPU+index);
    z = fetch_double(texrzboundaryGPU,nboundaryGPU+index);
    
    //Equilibrium distance 
    r0 = bFV->r0ParticleParticleGPU[offset+j];
    
    //Spring constant
    kSpring = bFV->kSpringParticleParticleGPU[offset+j];

    double rxx = rx-x;
    rxx -= floorf(rxx/lxGPU+0.5)*lxGPU;
    double ryy = ry-y;
    ryy -= floorf(ryy/lyGPU+0.5)*lyGPU;
    double rzz = rz-z;
    rzz -= floorf(rzz/lzGPU+0.5)*lzGPU;
    
    if(r0==0){
      fx += -kSpring * rxx;
      fy += -kSpring * ryy;
      fz += -kSpring * rzz;
    }
    else{     //If r0!=0 calculate particle particle distance
      r = sqrt(rxx*rxx+ryy*ryy+rzz*rzz);      
      if(r>0){//If r=0 -> f=0
	fx += -kSpring * (1 - r0/r) * rxx;
	fy += -kSpring * (1 - r0/r) * ryy;
	fz += -kSpring * (1 - r0/r) * rzz;
      }
    }
  }











  //Particle-FixedPoint Force
  nBonds = bFV->bondsParticleFixedPointGPU[i];
  offset = bFV->bondsParticleFixedPointOffsetGPU[i];
  
  
  for(int j=0;j<nBonds;j++){
    

    //Fixed point coordinates
    x = bFV->rxFixedPointGPU[offset+j];
    y = bFV->ryFixedPointGPU[offset+j];
    z = bFV->rzFixedPointGPU[offset+j];

    //*** Added by Adolfo ***
    // This is to fixed the fixed point springs to the wall
    x = x + xWall; 
    
    //Equilibrium distance 
    r0 = bFV->r0ParticleFixedPointGPU[offset+j];

    //Spring constant
    kSpring = bFV->kSpringParticleFixedPointGPU[offset+j];
    
    if(r0==0){
      fx += -kSpring * (rx - x);
      fy += -kSpring * (ry - y);
      fz += -kSpring * (rz - z);
      /**** Added by Adolfo ****/
      if (r0 != 0) // Only fixed points on the wall (the other ones has an eq dist = 0)
   	Fatt += kSpring * (rx - x);
	/********/
    }  
    else{     //If r0!=0 calculate particle particle distance
      r = sqrt( (x-rx)*(x-rx) + (y-ry)*(y-ry) + (z-rz)*(z-rz) );
      if(r>0){//If r=0 -> f=0
	fx += -kSpring * (1 - r0/r) * (rx - x);
	fy += -kSpring * (1 - r0/r) * (ry - y);
	fz += -kSpring * (1 - r0/r) * (rz - z);
	/****** Added by Adolfo ******/
	if (r0 != 0) // Only fixed points on the wall (the other ones has an eq dist = 0)
	  Fatt += kSpring * (1 - r0/r) * (rx - x);
	  /****************/   
      }
    }
  }

    

  /*** Added by Adolfo ***/
  //--- For checking  purposes 
  /*
  if (nBonds > 0)
     printf("i=%d\n %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",i, fx, Fatt, kSpring, r0, r, rx, ry, rz, x, y, z);
     */


  return ;
}
