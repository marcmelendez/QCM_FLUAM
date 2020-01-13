// Filename: quasiNeutrallyBuoyantFunctions.cu
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




//Raul added, a kernel to add a shear flow to the fluid. A sinusoidal force is added to each fluid cell. This shear flow is intended to be used for viscosity measurements.
__global__ void addShearFlowQuasiNeutrallyBuoyant(cufftDoubleComplex *vxZ,
					cufftDoubleComplex *vyZ,
					cufftDoubleComplex *vzZ,
					double viscosityMeasureAmplitude,
					double viscosityMeasureMode,
					int viscosityMeasurePlane,
					int viscosityMeasureDir){
  const int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=ncellsGPU) return;
  double posfluid;
  int ic;
  switch(viscosityMeasurePlane){
  case 0: ic = i%mxGPU;                    posfluid = (ic-0.5*mxGPU+0.5)/(double)mxGPU;    break;
  case 1: ic = (i%(mxGPU*mytGPU)) /mxGPU;  posfluid = (ic-0.5*mytGPU+0.5)/(double)mytGPU;  break;
  case 2: ic = i/(mxGPU*mytGPU);           posfluid = (ic-0.5*mzGPU+0.5)/(double)mzGPU;    break;
  }
  constexpr double pi2 = 6.28318530717958648;
  double force = viscosityMeasureAmplitude * sin(pi2*viscosityMeasureMode* posfluid)*dtGPU/densfluidGPU;
  switch(viscosityMeasureDir){
  case 0: vxZ[i].x -= force; break;
  case 1: vyZ[i].x -= force; break;
  case 2: vzZ[i].x -= force; break;
  }
}



//!*R New kernelSpread with Colors and three springs, tested with quasiNeutrally
//Fill "countparticlesincellX" lists
//and spread particle force F, The force comes from LJ function in nonBondedForce.cu
__global__ void kernelSpreadParticlesForceColors(const double* rxcellGPU,
						 const double* rycellGPU,
						 const double* rzcellGPU,
						 double* fxboundaryGPU,
						 double* fyboundaryGPU,
						 double* fzboundaryGPU,
						 particlesincell* pc,
						 int* errorKernel,
						 const bondedForcesVariables* bFV,
						 const particle_type *pt,
						 const threeParticleBondsVariables* tPBV){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;

  int type_i = pt->types[nboundaryGPU+i];
  int type_j;


  double fx = 0.;//pressurea0GPU;// * 0.267261241912424385 ;//0;
  double fy = 0.;//pressurea0GPU * 0.534522483824848769 ;
  double fz = 0.;//pressurea0GPU * 0.801783725737273154 ;
  double f;


  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  //INCLUDE EXTERNAL FORCES HERE

  //Example: harmonic potential
  // V(r) = (1/2) * k * ((x-x0)**2 + (y-y0)**2 + (z-z0)**2)
  //
  //with spring constant k=0.01
  //and x0=y0=z0=0
  //
  //fx = -0.01*rx;
  //fy = -0.01*ry;
  //fz = -0.01*rz;


  if(confinementZGPU){
    fz = confinementZKGPU*rz;
  }

  /* CODE MODIFIED by Marc: QCM vibrating wall */
  double time = nstepGPU*dtGPU;
  if(i < nwall2GPU){
  //  if(time < stokesFlowTimeGPU)
  //    fx = -stokesFlowAmpGPU*sin(6.2831853*time*stokesFlowFreqGPU);
    if(i < nwall1GPU)
      fz = -kwall1GPU*(rz - zwall1GPU);
    else
      fz = -kwall2GPU*(rz - zwall2GPU);
    }
  else if(i < nwall3GPU) {
    if(time < stokesFlowTimeGPU)
      fx = stokesFlowAmpGPU*sin(6.2831853*time*stokesFlowFreqGPU);
    fz = -kwall3GPU*(rz - zwall3GPU);
  }
  /* END OF MODIFICATION by Marc */

  if(particlesWallGPU){
    //INCLUDE WALL REPULSION HERE
    //We use a repulsive Lennard-Jones
    double sigmaWall = 2*dyGPU;
    double cutoffWall = 1.12246204830937302 * sigmaWall; // 2^{1/6} * sigmaWall

    //IMPORTANT, for particleWall lyGPU stores ly+2*dy
    if(ry<(-0.5*lyGPU+cutoffWall+dyGPU)){//Left wall

      double distance = (0.5*lyGPU-dyGPU) + ry; //distance >= 0
      fy += 48 * temperatureGPU * (pow((sigmaWall/distance),13) - 0.5*pow((sigmaWall/distance),7));
    }
    else if(ry>(0.5*lyGPU-cutoffWall-dyGPU)){//Right wall
      double distance = (0.5*lyGPU-dyGPU) - ry; //distance >= 0
      fy -= 48 * temperatureGPU * (pow((sigmaWall/distance),13) - 0.5*pow((sigmaWall/distance),7));
    }
  }





  //NEW bonded forces
  if(bondedForcesGPU){
    //call function for bonded forces particle-particle
    forceBondedParticleParticleGPU(i,
				   fx,
				   fy,
				   fz,
				   rx,
				   ry,
				   rz,
				   bFV);
  }
  if(threeBondedForcesGPU){
    forceBondedThreeParticleGPU(i,
				fx,
				fy,
				fz,
				rx,
				ry,
				rz,
				tPBV);
  }

    
  double rxij, ryij, rzij, r2;

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;
  
  int icel;
  double r, rp, rm;

  {
    double invdx = double(mxNeighborsGPU)/lxGPU;
    double invdy = double(myNeighborsGPU)/lyGPU;
    double invdz = double(mzNeighborsGPU)/lzGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdx + 0.5*mxNeighborsGPU) % mxNeighborsGPU;
    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdy + 0.5*myNeighborsGPU) % myNeighborsGPU;
    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdz + 0.5*mzNeighborsGPU) % mzNeighborsGPU;
    icel  = jx;
    icel += jy * mxNeighborsGPU;
    icel += jz * mxNeighborsGPU * myNeighborsGPU;    
  }
  
  int np;
  if(computeNonBondedForcesGPU){
    //Particles in Cell i
    np = tex1Dfetch(texCountParticlesInCellNonBonded,icel);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+icel);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }  
    //Particles in Cell vecino0
    vecino0 = tex1Dfetch(texneighbor0GPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino0);
    //printf("RRRRRRRRR %i %i %i\n",i,icel,vecino0);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino0);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecino1
    vecino1 = tex1Dfetch(texneighbor1GPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino1);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino1);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecino2
    vecino2 = tex1Dfetch(texneighbor2GPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino2);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino2);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecino3
    vecino3 = tex1Dfetch(texneighbor3GPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino3);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino3);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecino4
    vecino4 = tex1Dfetch(texneighbor4GPU, icel);
    //printf("VECINO %i %i \n",icel,vecino4);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino4);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino4);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecino5
    vecino5 = tex1Dfetch(texneighbor5GPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecino5);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecino5);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxpy
    vecinopxpy = tex1Dfetch(texneighborpxpyGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxpy);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxpy);    
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxmy
    vecinopxmy = tex1Dfetch(texneighborpxmyGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxmy);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxmy);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxpz
    vecinopxpz = tex1Dfetch(texneighborpxpzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxpz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxpz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxmz
    vecinopxmz = tex1Dfetch(texneighborpxmzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxmz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxmz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxpy
    vecinomxpy = tex1Dfetch(texneighbormxpyGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxpy);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxpy);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij; 
    }
    //Particles in Cell vecinomxmy
    vecinomxmy = tex1Dfetch(texneighbormxmyGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxmy);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxmy);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxpz
    vecinomxpz = tex1Dfetch(texneighbormxpzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxpz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxpz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxmz
    vecinomxmz = tex1Dfetch(texneighbormxmzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxmz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxmz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopypz
    vecinopypz = tex1Dfetch(texneighborpypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopymz
    vecinopymz = tex1Dfetch(texneighborpymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomypz
    vecinomypz = tex1Dfetch(texneighbormypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomymz
    vecinomymz = tex1Dfetch(texneighbormymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxpypz
    vecinopxpypz = tex1Dfetch(texneighborpxpypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxpypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxpypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxpymz
    vecinopxpymz = tex1Dfetch(texneighborpxpymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxpymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxpymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxmypz
    vecinopxmypz = tex1Dfetch(texneighborpxmypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxmypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxmypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinopxmymz
    vecinopxmymz = tex1Dfetch(texneighborpxmymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinopxmymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinopxmymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxpypz
    vecinomxpypz = tex1Dfetch(texneighbormxpypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxpypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxpypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxpymz
    vecinomxpymz = tex1Dfetch(texneighbormxpymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxpymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxpymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxmypz
    vecinomxmypz = tex1Dfetch(texneighbormxmypzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxmypz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxmypz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    //Particles in Cell vecinomxmymz
    vecinomxmymz = tex1Dfetch(texneighbormxmymzGPU, icel);
    np = tex1Dfetch(texCountParticlesInCellNonBonded,vecinomxmymz);
    for(int j=0;j<np;j++){
      int particle = tex1Dfetch(texPartInCellNonBonded,mNeighborsGPU*j+vecinomxmymz);
      rxij =  (rx - fetch_double(texrxboundaryGPU,particle));
      rxij =  (rxij - int(rxij*invlxGPU + 0.5*((rxij>0)-(rxij<0)))*lxGPU);
      ryij =  (ry - fetch_double(texryboundaryGPU,particle));
      ryij =  (ryij - int(ryij*invlyGPU + 0.5*((ryij>0)-(ryij<0)))*lyGPU);
      rzij =  (rz - fetch_double(texrzboundaryGPU,particle));
      rzij =  (rzij - int(rzij*invlzGPU + 0.5*((rzij>0)-(rzij<0)))*lzGPU);
      r2 = rxij*rxij + ryij*ryij + rzij*rzij;
      type_j = pt->types[particle];f =LJ(r2, pt->Aij_param, pt->Bij_param, type_i+ntypesGPU*type_j, nboundaryGPU+i, particle);
      fx += f * rxij;
      fy += f * ryij;
      fz += f * rzij;
    }
    
  }
  
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;
  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * mytGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }

  np = atomicAdd(&pc->countparticlesincellX[icelx],1);
  if(np >= maxNumberPartInCellGPU){
    errorKernel[0]=1;
    errorKernel[1]=maxNumberPartInCellGPU;
    return;
  }
  pc->partincellX[ncellstGPU*np+icelx] = nboundaryGPU+i;
  np = atomicAdd(&pc->countparticlesincellY[icely],1);
  if(np >= maxNumberPartInCellGPU){
    errorKernel[0]=1;
    errorKernel[2]=np;
    return;
  }
  pc->partincellY[ncellstGPU*np+icely] = nboundaryGPU+i;
  np = atomicAdd(&pc->countparticlesincellZ[icelz],1);
  if(np >= maxNumberPartInCellGPU){
    errorKernel[0]=1;
    errorKernel[3]=np;
    return;
  }
  pc->partincellZ[ncellstGPU*np+icelz] = nboundaryGPU+i;

  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;

  //FORCE IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);

  //double dlxS, dlxpS, dlxmS;
  //double dlyS, dlypS, dlymS;
  //double dlzS, dlzpS, dlzmS;

  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx  = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
  //dlxS = functionDeltaDerived(1.5*r);
  //dlxpS = functionDeltaDerived(1.5*rp);
  //dlxmS = functionDeltaDerived(1.5*rm);
                              
  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);
  //dlyS = functionDeltaDerived(1.5*r);
  //dlypS = functionDeltaDerived(1.5*rp);
  //dlymS = functionDeltaDerived(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  //dlzS = functionDeltaDerived(1.5*r);
  //dlzpS = functionDeltaDerived(1.5*rp);
  //dlzmS = functionDeltaDerived(1.5*rm);

    
  int offset = nboundaryGPU;
  fxboundaryGPU[offset+i] = -dlxp * dlyp * dlzp * fx;
  offset += nboundaryGPU+npGPU;//1
  fxboundaryGPU[offset+i] = -dlxp * dlyp * dlz  * fx;
  offset += nboundaryGPU+npGPU;//2
  fxboundaryGPU[offset+i] = -dlxp * dlyp * dlzm * fx;
  offset += nboundaryGPU+npGPU;//3
  fxboundaryGPU[offset+i] = -dlxp * dly  * dlzp * fx;
  offset += nboundaryGPU+npGPU;//4
  fxboundaryGPU[offset+i] = -dlxp * dly  * dlz  * fx;
  offset += nboundaryGPU+npGPU;//5
  fxboundaryGPU[offset+i] = -dlxp * dly  * dlzm * fx;
  offset += nboundaryGPU+npGPU;//6
  fxboundaryGPU[offset+i] = -dlxp * dlym * dlzp * fx;
  offset += nboundaryGPU+npGPU;//7
  fxboundaryGPU[offset+i] = -dlxp * dlym * dlz  * fx;
  offset += nboundaryGPU+npGPU;//8
  fxboundaryGPU[offset+i] = -dlxp * dlym * dlzm * fx;
  offset += nboundaryGPU+npGPU;//9
  fxboundaryGPU[offset+i] = -dlx  * dlyp * dlzp * fx;
  offset += nboundaryGPU+npGPU;//10
  fxboundaryGPU[offset+i] = -dlx  * dlyp * dlz  * fx;
  offset += nboundaryGPU+npGPU;//11
  fxboundaryGPU[offset+i] = -dlx  * dlyp * dlzm * fx;
  offset += nboundaryGPU+npGPU;//12
  fxboundaryGPU[offset+i] = -dlx  * dly  * dlzp * fx;
  offset += nboundaryGPU+npGPU;//13
  fxboundaryGPU[offset+i] = -dlx  * dly  * dlz  * fx;
  offset += nboundaryGPU+npGPU;//14
  fxboundaryGPU[offset+i] = -dlx  * dly  * dlzm * fx;
  offset += nboundaryGPU+npGPU;//15
  fxboundaryGPU[offset+i] = -dlx  * dlym * dlzp * fx;
  offset += nboundaryGPU+npGPU;//16
  fxboundaryGPU[offset+i] = -dlx  * dlym * dlz  * fx;
  offset += nboundaryGPU+npGPU;//17
  fxboundaryGPU[offset+i] = -dlx  * dlym * dlzm * fx;
  offset += nboundaryGPU+npGPU;//18
  fxboundaryGPU[offset+i] = -dlxm * dlyp * dlzp * fx;
  offset += nboundaryGPU+npGPU;//19
  fxboundaryGPU[offset+i] = -dlxm * dlyp * dlz  * fx;
  offset += nboundaryGPU+npGPU;//20
  fxboundaryGPU[offset+i] = -dlxm * dlyp * dlzm * fx;
  offset += nboundaryGPU+npGPU;//21
  fxboundaryGPU[offset+i] = -dlxm * dly  * dlzp * fx;
  offset += nboundaryGPU+npGPU;//22
  fxboundaryGPU[offset+i] = -dlxm * dly  * dlz  * fx;
  offset += nboundaryGPU+npGPU;//23
  fxboundaryGPU[offset+i] = -dlxm * dly  * dlzm * fx;
  offset += nboundaryGPU+npGPU;//24
  fxboundaryGPU[offset+i] = -dlxm * dlym * dlzp * fx;
  offset += nboundaryGPU+npGPU;//25
  fxboundaryGPU[offset+i] = -dlxm * dlym * dlz  * fx;
  offset += nboundaryGPU+npGPU;//26
  fxboundaryGPU[offset+i] = -dlxm * dlym * dlzm * fx;
  
  
  //FORCE IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
  //dlxS = functionDeltaDerived(1.5*r);
  //dlxpS = functionDeltaDerived(1.5*rp);
  //dlxmS = functionDeltaDerived(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);
  //dlyS = functionDeltaDerived(1.5*r);
  //dlypS = functionDeltaDerived(1.5*rp);
  //dlymS = functionDeltaDerived(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  //dlzS = functionDeltaDerived(1.5*r);
  //dlzpS = functionDeltaDerived(1.5*rp);
  //dlzmS = functionDeltaDerived(1.5*rm);

  offset = nboundaryGPU;
  fyboundaryGPU[offset+i] = -dlxp * dlyp * dlzp * fy;
  offset += nboundaryGPU+npGPU;//1
  fyboundaryGPU[offset+i] = -dlxp * dlyp * dlz  * fy;
  offset += nboundaryGPU+npGPU;//2
  fyboundaryGPU[offset+i] = -dlxp * dlyp * dlzm * fy;
  offset += nboundaryGPU+npGPU;//3
  fyboundaryGPU[offset+i] = -dlxp * dly  * dlzp * fy;
  offset += nboundaryGPU+npGPU;//4
  fyboundaryGPU[offset+i] = -dlxp * dly  * dlz  * fy;
  offset += nboundaryGPU+npGPU;//5
  fyboundaryGPU[offset+i] = -dlxp * dly  * dlzm * fy;
  offset += nboundaryGPU+npGPU;//6
  fyboundaryGPU[offset+i] = -dlxp * dlym * dlzp * fy;
  offset += nboundaryGPU+npGPU;//7
  fyboundaryGPU[offset+i] = -dlxp * dlym * dlz  * fy;
  offset += nboundaryGPU+npGPU;//8
  fyboundaryGPU[offset+i] = -dlxp * dlym * dlzm * fy;
  offset += nboundaryGPU+npGPU;//9
  fyboundaryGPU[offset+i] = -dlx  * dlyp * dlzp * fy;
  offset += nboundaryGPU+npGPU;//10
  fyboundaryGPU[offset+i] = -dlx  * dlyp * dlz  * fy;
  offset += nboundaryGPU+npGPU;//11
  fyboundaryGPU[offset+i] = -dlx  * dlyp * dlzm * fy;
  offset += nboundaryGPU+npGPU;//12
  fyboundaryGPU[offset+i] = -dlx  * dly  * dlzp * fy;
  offset += nboundaryGPU+npGPU;//13
  fyboundaryGPU[offset+i] = -dlx  * dly  * dlz  * fy;
  offset += nboundaryGPU+npGPU;//14
  fyboundaryGPU[offset+i] = -dlx  * dly  * dlzm * fy;
  offset += nboundaryGPU+npGPU;//15
  fyboundaryGPU[offset+i] = -dlx  * dlym * dlzp * fy;
  offset += nboundaryGPU+npGPU;//16
  fyboundaryGPU[offset+i] = -dlx  * dlym * dlz  * fy;
  offset += nboundaryGPU+npGPU;//17
  fyboundaryGPU[offset+i] = -dlx  * dlym * dlzm * fy;
  offset += nboundaryGPU+npGPU;//18
  fyboundaryGPU[offset+i] = -dlxm * dlyp * dlzp * fy;
  offset += nboundaryGPU+npGPU;//19
  fyboundaryGPU[offset+i] = -dlxm * dlyp * dlz  * fy;
  offset += nboundaryGPU+npGPU;//20
  fyboundaryGPU[offset+i] = -dlxm * dlyp * dlzm * fy;
  offset += nboundaryGPU+npGPU;//21
  fyboundaryGPU[offset+i] = -dlxm * dly  * dlzp * fy;
  offset += nboundaryGPU+npGPU;//22
  fyboundaryGPU[offset+i] = -dlxm * dly  * dlz  * fy;
  offset += nboundaryGPU+npGPU;//23
  fyboundaryGPU[offset+i] = -dlxm * dly  * dlzm * fy;
  offset += nboundaryGPU+npGPU;//24
  fyboundaryGPU[offset+i] = -dlxm * dlym * dlzp * fy;
  offset += nboundaryGPU+npGPU;//25
  fyboundaryGPU[offset+i] = -dlxm * dlym * dlz  * fy;
  offset += nboundaryGPU+npGPU;//26
  fyboundaryGPU[offset+i] = -dlxm * dlym * dlzm * fy;


  
  //FORCE IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
  //dlxS = functionDeltaDerived(1.5*r);
  //dlxpS = functionDeltaDerived(1.5*rp);
  //dlxmS = functionDeltaDerived(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);
  //dlyS = functionDeltaDerived(1.5*r);
  //dlypS = functionDeltaDerived(1.5*rp);
  //dlymS = functionDeltaDerived(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  //dlzS = functionDeltaDerived(1.5*r);
  //dlzpS = functionDeltaDerived(1.5*rp);
  //dlzmS = functionDeltaDerived(1.5*rm);


  offset = nboundaryGPU;
  fzboundaryGPU[offset+i] = -dlxp * dlyp * dlzp * fz;
  offset += nboundaryGPU+npGPU;//1
  fzboundaryGPU[offset+i] = -dlxp * dlyp * dlz  * fz;
  offset += nboundaryGPU+npGPU;//2
  fzboundaryGPU[offset+i] = -dlxp * dlyp * dlzm * fz;
  offset += nboundaryGPU+npGPU;//3
  fzboundaryGPU[offset+i] = -dlxp * dly  * dlzp * fz;
  offset += nboundaryGPU+npGPU;//4
  fzboundaryGPU[offset+i] = -dlxp * dly  * dlz  * fz;
  offset += nboundaryGPU+npGPU;//5
  fzboundaryGPU[offset+i] = -dlxp * dly  * dlzm * fz;
  offset += nboundaryGPU+npGPU;//6
  fzboundaryGPU[offset+i] = -dlxp * dlym * dlzp * fz;
  offset += nboundaryGPU+npGPU;//7
  fzboundaryGPU[offset+i] = -dlxp * dlym * dlz  * fz;
  offset += nboundaryGPU+npGPU;//8
  fzboundaryGPU[offset+i] = -dlxp * dlym * dlzm * fz;
  offset += nboundaryGPU+npGPU;//9
  fzboundaryGPU[offset+i] = -dlx  * dlyp * dlzp * fz;
  offset += nboundaryGPU+npGPU;//10
  fzboundaryGPU[offset+i] = -dlx  * dlyp * dlz  * fz;
  offset += nboundaryGPU+npGPU;//11
  fzboundaryGPU[offset+i] = -dlx  * dlyp * dlzm * fz;
  offset += nboundaryGPU+npGPU;//12
  fzboundaryGPU[offset+i] = -dlx  * dly  * dlzp * fz;
  offset += nboundaryGPU+npGPU;//13
  fzboundaryGPU[offset+i] = -dlx  * dly  * dlz  * fz;
  offset += nboundaryGPU+npGPU;//14
  fzboundaryGPU[offset+i] = -dlx  * dly  * dlzm * fz;
  offset += nboundaryGPU+npGPU;//15
  fzboundaryGPU[offset+i] = -dlx  * dlym * dlzp * fz;
  offset += nboundaryGPU+npGPU;//16
  fzboundaryGPU[offset+i] = -dlx  * dlym * dlz  * fz;
  offset += nboundaryGPU+npGPU;//17
  fzboundaryGPU[offset+i] = -dlx  * dlym * dlzm * fz;
  offset += nboundaryGPU+npGPU;//18
  fzboundaryGPU[offset+i] = -dlxm * dlyp * dlzp * fz;
  offset += nboundaryGPU+npGPU;//19
  fzboundaryGPU[offset+i] = -dlxm * dlyp * dlz  * fz;
  offset += nboundaryGPU+npGPU;//20
  fzboundaryGPU[offset+i] = -dlxm * dlyp * dlzm * fz;
  offset += nboundaryGPU+npGPU;//21
  fzboundaryGPU[offset+i] = -dlxm * dly  * dlzp * fz;
  offset += nboundaryGPU+npGPU;//22
  fzboundaryGPU[offset+i] = -dlxm * dly  * dlz  * fz;
  offset += nboundaryGPU+npGPU;//23
  fzboundaryGPU[offset+i] = -dlxm * dly  * dlzm * fz;
  offset += nboundaryGPU+npGPU;//24
  fzboundaryGPU[offset+i] = -dlxm * dlym * dlzp * fz;
  offset += nboundaryGPU+npGPU;//25
  fzboundaryGPU[offset+i] = -dlxm * dlym * dlz  * fz;
  offset += nboundaryGPU+npGPU;//26
  fzboundaryGPU[offset+i] = -dlxm * dlym * dlzm * fz;
  
}






__global__ void copyParticles(double *rx,
			      double *ry,
			      double *rz,
			      double *rxCopy,
			      double *ryCopy,
			      double *rzCopy){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  rxCopy[i] = rx[i];
  ryCopy[i] = ry[i];
  rzCopy[i] = rz[i];

}





__global__ void copyField(double *vx1,
			  double *vy1,
			  double *vz1,
			  double *vx2,
			  double *vy2,
			  double *vz2){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(ncellstGPU)) return;   

  vx2[i] = vx1[i];
  vy2[i] = vy1[i];
  vz2[i] = vz1[i];

}







__global__ void meanField(double *vx1,
			  double *vy1,
			  double *vz1,
			  cufftDoubleComplex *vx2,
			  cufftDoubleComplex *vy2,
			  cufftDoubleComplex *vz2){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(ncellsGPU)) return;   

  vx1[i] = 0.5 * ( vx1[i] + vx2[i].x / double(ncellsGPU) );
  vy1[i] = 0.5 * ( vy1[i] + vy2[i].x / double(ncellsGPU) );
  vz1[i] = 0.5 * ( vz1[i] + vz2[i].x / double(ncellsGPU) );

}





__global__ void updateParticlesTest(double *rxcellGPU,
				    double *rycellGPU,
				    double *rzcellGPU,
				    double *rxboundaryGPU,
				    double *ryboundaryGPU,
				    double *rzboundaryGPU,
				    double *rxboundaryPredictionGPU,
				    double *ryboundaryPredictionGPU,
				    double *rzboundaryPredictionGPU,
				    double *vxPredictionGPU,
				    double *vyPredictionGPU,
				    double *vzPredictionGPU){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  





  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  
  double v = dlxm * dlym * dlzm * vxPredictionGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxPredictionGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxPredictionGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxPredictionGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxPredictionGPU[vecino2] +
    dlxm * dly  * dlzp * vxPredictionGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxPredictionGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxPredictionGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxPredictionGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxPredictionGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxPredictionGPU[vecino1] +
    dlx  * dlym * dlzp * vxPredictionGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxPredictionGPU[vecino0] + 
    dlx  * dly  * dlz  * vxPredictionGPU[icelx] + 
    dlx  * dly  * dlzp * vxPredictionGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxPredictionGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxPredictionGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxPredictionGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxPredictionGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxPredictionGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxPredictionGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxPredictionGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxPredictionGPU[vecino3] +
    dlxp * dly  * dlzp * vxPredictionGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxPredictionGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxPredictionGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxPredictionGPU[vecinopxpypz];


  rxboundaryGPU[i] = rxboundaryPredictionGPU[i] + dtGPU * v;
  

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  v = dlxm * dlym * dlzm * vyPredictionGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyPredictionGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyPredictionGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyPredictionGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyPredictionGPU[vecino2] +
    dlxm * dly  * dlzp * vyPredictionGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyPredictionGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyPredictionGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyPredictionGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyPredictionGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyPredictionGPU[vecino1] +
    dlx  * dlym * dlzp * vyPredictionGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyPredictionGPU[vecino0] + 
    dlx  * dly  * dlz  * vyPredictionGPU[icely] + 
    dlx  * dly  * dlzp * vyPredictionGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyPredictionGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyPredictionGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyPredictionGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyPredictionGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyPredictionGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyPredictionGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyPredictionGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyPredictionGPU[vecino3] +
    dlxp * dly  * dlzp * vyPredictionGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyPredictionGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyPredictionGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyPredictionGPU[vecinopxpypz];

  
  ryboundaryGPU[i] = ryboundaryPredictionGPU[i] + dtGPU * v;
  

  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  v = dlxm * dlym * dlzm * vzPredictionGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzPredictionGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzPredictionGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzPredictionGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzPredictionGPU[vecino2] +
    dlxm * dly  * dlzp * vzPredictionGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzPredictionGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzPredictionGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzPredictionGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzPredictionGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzPredictionGPU[vecino1] +
    dlx  * dlym * dlzp * vzPredictionGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzPredictionGPU[vecino0] + 
    dlx  * dly  * dlz  * vzPredictionGPU[icelz] + 
    dlx  * dly  * dlzp * vzPredictionGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzPredictionGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzPredictionGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzPredictionGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzPredictionGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzPredictionGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzPredictionGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzPredictionGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzPredictionGPU[vecino3] +
    dlxp * dly  * dlzp * vzPredictionGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzPredictionGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzPredictionGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzPredictionGPU[vecinopxpypz];


  rzboundaryGPU[i] = rzboundaryPredictionGPU[i] + dtGPU * v;

}













__global__ void updateParticlesTest2(double *rxcellGPU,
				     double *rycellGPU,
				     double *rzcellGPU,
				     double *vxboundaryGPU,
				     double *vyboundaryGPU,
				     double *vzboundaryGPU,
				     double *vxGPU,
				     double *vyGPU,
				     double *vzGPU){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  





  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];


  vxboundaryGPU[i] = v;

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];

  
  vyboundaryGPU[i] = v;

  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];

  vzboundaryGPU[i] = v;


}



























__global__ void kernelCorrectionVQuasiNeutrallyBuoyantSemiImplicit_1_TEST3(double* rxcellGPU,
									   double* rycellGPU,
									   double* rzcellGPU,
									   double* vxGPU,
									   double* vyGPU,
									   double* vzGPU,
									   double* fxboundaryGPU,
									   double* fyboundaryGPU,
									   double* fzboundaryGPU,
									   double* vxboundaryGPU,
									   double* vyboundaryGPU,
									   double* vzboundaryGPU){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;
  int icelx, icely, icelz;

  double r, rp, rm;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }

  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;

  //CORRECTION IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);

  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx  = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
                              
  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  double v = dlxm * dlym * dlzm * fetch_double(texVxGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVxGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVxGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVxGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVxGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVxGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVxGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVxGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVxGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVxGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVxGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVxGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVxGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVxGPU,icelx) + 
    dlx  * dly  * dlzp * fetch_double(texVxGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVxGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVxGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVxGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVxGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVxGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVxGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVxGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVxGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVxGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVxGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVxGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVxGPU,vecinopxpypz);

  double u = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];



  double fact = massParticleGPU / (volumeGPU * densfluidGPU);
  //double deltaV = fact * (vxboundaryGPU[i] - v);
  double deltaV = fact * (u - v);

  int offset = nboundaryGPU;
  fxboundaryGPU[offset+i]        = dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fxboundaryGPU[offset+i] = dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fxboundaryGPU[offset+i] = dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fxboundaryGPU[offset+i] = dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fxboundaryGPU[offset+i] = dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fxboundaryGPU[offset+i] = dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fxboundaryGPU[offset+i] = dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fxboundaryGPU[offset+i] = dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fxboundaryGPU[offset+i] = dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fxboundaryGPU[offset+i] = dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fxboundaryGPU[offset+i] = dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fxboundaryGPU[offset+i] = dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fxboundaryGPU[offset+i] = dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fxboundaryGPU[offset+i] = dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fxboundaryGPU[offset+i] = dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fxboundaryGPU[offset+i] = dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fxboundaryGPU[offset+i] = dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fxboundaryGPU[offset+i] = dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fxboundaryGPU[offset+i] = dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fxboundaryGPU[offset+i] = dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fxboundaryGPU[offset+i] = dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fxboundaryGPU[offset+i] = dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fxboundaryGPU[offset+i] = dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fxboundaryGPU[offset+i] = dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fxboundaryGPU[offset+i] = dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fxboundaryGPU[offset+i] = dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fxboundaryGPU[offset+i] = dlxm * dlym * dlzm * deltaV;


  //CORRECTION IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  v = dlxm * dlym * dlzm * fetch_double(texVyGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVyGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVyGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVyGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVyGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVyGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVyGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVyGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVyGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVyGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVyGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVyGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVyGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVyGPU,icely) + 
    dlx  * dly  * dlzp * fetch_double(texVyGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVyGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVyGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVyGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVyGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVyGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVyGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVyGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVyGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVyGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVyGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVyGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVyGPU,vecinopxpypz);

  u = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];


  //deltaV = fact * v;
  //deltaV = fact * (vyboundaryGPU[i] - v);
  deltaV = fact * (u - v);

  offset = nboundaryGPU;
  fyboundaryGPU[offset+i]        = dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fyboundaryGPU[offset+i] = dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fyboundaryGPU[offset+i] = dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fyboundaryGPU[offset+i] = dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fyboundaryGPU[offset+i] = dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fyboundaryGPU[offset+i] = dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fyboundaryGPU[offset+i] = dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fyboundaryGPU[offset+i] = dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fyboundaryGPU[offset+i] = dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fyboundaryGPU[offset+i] = dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fyboundaryGPU[offset+i] = dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fyboundaryGPU[offset+i] = dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fyboundaryGPU[offset+i] = dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fyboundaryGPU[offset+i] = dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fyboundaryGPU[offset+i] = dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fyboundaryGPU[offset+i] = dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fyboundaryGPU[offset+i] = dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fyboundaryGPU[offset+i] = dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fyboundaryGPU[offset+i] = dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fyboundaryGPU[offset+i] = dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fyboundaryGPU[offset+i] = dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fyboundaryGPU[offset+i] = dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fyboundaryGPU[offset+i] = dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fyboundaryGPU[offset+i] = dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fyboundaryGPU[offset+i] = dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fyboundaryGPU[offset+i] = dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fyboundaryGPU[offset+i] = dlxm * dlym * dlzm * deltaV;


  //CORRECTION IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * fetch_double(texVzGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVzGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVzGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVzGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVzGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVzGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVzGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVzGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVzGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVzGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVzGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVzGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVzGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVzGPU,icelz) + 
    dlx  * dly  * dlzp * fetch_double(texVzGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVzGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVzGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVzGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVzGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVzGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVzGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVzGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVzGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVzGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVzGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVzGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVzGPU,vecinopxpypz);

  u = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];


  //deltaV = fact * v;
  //deltaV = fact * (vzboundaryGPU[i] - v);
  deltaV = fact * (u - v);

  offset = nboundaryGPU;
  fzboundaryGPU[offset+i]        = dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fzboundaryGPU[offset+i] = dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fzboundaryGPU[offset+i] = dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fzboundaryGPU[offset+i] = dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fzboundaryGPU[offset+i] = dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fzboundaryGPU[offset+i] = dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fzboundaryGPU[offset+i] = dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fzboundaryGPU[offset+i] = dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fzboundaryGPU[offset+i] = dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fzboundaryGPU[offset+i] = dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fzboundaryGPU[offset+i] = dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fzboundaryGPU[offset+i] = dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fzboundaryGPU[offset+i] = dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fzboundaryGPU[offset+i] = dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fzboundaryGPU[offset+i] = dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fzboundaryGPU[offset+i] = dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fzboundaryGPU[offset+i] = dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fzboundaryGPU[offset+i] = dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fzboundaryGPU[offset+i] = dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fzboundaryGPU[offset+i] = dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fzboundaryGPU[offset+i] = dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fzboundaryGPU[offset+i] = dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fzboundaryGPU[offset+i] = dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fzboundaryGPU[offset+i] = dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fzboundaryGPU[offset+i] = dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fzboundaryGPU[offset+i] = dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fzboundaryGPU[offset+i] = dlxm * dlym * dlzm * deltaV;


}



























__global__ void findNeighborParticlesQuasiNeutrallyBuoyantTEST3_1(particlesincell* pc, 
								  int* errorKernel,
								  double* rxcellGPU,
								  double* rycellGPU,
								  double* rzcellGPU,
								  double* rxboundaryGPU,  //q^{n-1/2}
								  double* ryboundaryGPU, 
								  double* rzboundaryGPU,
								  double* rxboundaryPredictionGPU, //q^{*,n+1/2}
								  double* ryboundaryPredictionGPU, 
								  double* rzboundaryPredictionGPU,
								  double* vxGPU, // tilde{v^n}
								  double* vyGPU, 
								  double* vzGPU){
								  
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];


  rxboundaryPredictionGPU[nboundaryGPU+i] = rxboundaryGPU[nboundaryGPU+i] + 
    dtGPU * v ;

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  
  ryboundaryPredictionGPU[nboundaryGPU+i] = ryboundaryGPU[nboundaryGPU+i] +
    dtGPU * v ;

  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];

  rzboundaryPredictionGPU[nboundaryGPU+i] = rzboundaryGPU[nboundaryGPU+i] + 
    dtGPU * v ;

}

















__global__ void findNeighborParticlesQuasiNeutrallyBuoyantTEST3_2(particlesincell* pc, 
								  int* errorKernel,
								  double* rxcellGPU,
								  double* rycellGPU,
								  double* rzcellGPU,
								  double* rxboundaryGPU,  //q^{n+1/2}
								  double* ryboundaryGPU, 
								  double* rzboundaryGPU,
								  double* rxboundaryPredictionGPU, //q^{*,n+1/2}
								  double* ryboundaryPredictionGPU, 
								  double* rzboundaryPredictionGPU,
								  double* vxGPU, //v^n
								  double* vyGPU, 
								  double* vzGPU){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = rxboundaryPredictionGPU[nboundaryGPU+i];
  double ry = ryboundaryPredictionGPU[nboundaryGPU+i];
  double rz = rzboundaryPredictionGPU[nboundaryGPU+i];
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];


  double rxNew = 0.5 * (rxboundaryGPU[nboundaryGPU+i] +
			rx + dtGPU * v);

  rxboundaryGPU[nboundaryGPU+i] = rxNew;

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  
 
  double ryNew = 0.5 * (ryboundaryGPU[nboundaryGPU+i] +
			ry + dtGPU * v);

  ryboundaryGPU[nboundaryGPU+i] = ryNew;

  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];


  double rzNew = 0.5 * (rzboundaryGPU[nboundaryGPU+i] +
			rz + dtGPU * v);

  rzboundaryGPU[nboundaryGPU+i] = rzNew;



  int icel;
  {
    double invdx = double(mxNeighborsGPU)/lxGPU;
    double invdy = double(myNeighborsGPU)/lyGPU;
    double invdz = double(mzNeighborsGPU)/lzGPU;
    rxNew = rxNew - (int(rxNew*invlxGPU + 0.5*((rxNew>0)-(rxNew<0)))) * lxGPU;
    int jx   = int(rxNew * invdx + 0.5*mxNeighborsGPU) % mxNeighborsGPU;
    ryNew = ryNew - (int(ryNew*invlyGPU + 0.5*((ryNew>0)-(ryNew<0)))) * lyGPU;
    int jy   = int(ryNew * invdy + 0.5*myNeighborsGPU) % myNeighborsGPU;
    rzNew = rzNew - (int(rzNew*invlzGPU + 0.5*((rzNew>0)-(rzNew<0)))) * lzGPU;
    int jz   = int(rzNew * invdz + 0.5*mzNeighborsGPU) % mzNeighborsGPU;
    icel  = jx;
    icel += jy * mxNeighborsGPU;
    icel += jz * mxNeighborsGPU * myNeighborsGPU;    
  }
  int np = atomicAdd(&pc->countPartInCellNonBonded[icel],1);
  if(np >= maxNumberPartInCellNonBondedGPU){
    errorKernel[0] = 1;
    errorKernel[4] = 1;
    return;
  }
  pc->partInCellNonBonded[mNeighborsGPU*np+icel] = i;



}

























__global__ void findNeighborParticlesQuasiNeutrallyBuoyant_1(particlesincell* pc, 
							     int* errorKernel,
							     const double* rxcellGPU,
							     const double* rycellGPU,
							     const double* rzcellGPU,
							     const double* rxboundaryGPU,  //q^{n}
							     const double* ryboundaryGPU, 
							     const double* rzboundaryGPU,
							     double* rxboundaryPredictionGPU, //q^{n+1/2}
							     double* ryboundaryPredictionGPU, 
							     double* rzboundaryPredictionGPU,
							     const double* vxGPU, //v^n
							     const double* vyGPU, 
							     const double* vzGPU){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * mytGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
  
  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];

  
  double rxNew = rx + 0.5 * dtGPU * v;
  rxboundaryPredictionGPU[nboundaryGPU+i] = rxNew;


  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  

  double ryNew = ry + 0.5 * dtGPU * v;

  ryboundaryPredictionGPU[nboundaryGPU+i] = ryNew;
 
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];



  double rzNew = rz + 0.5 * dtGPU * v;
  rzboundaryPredictionGPU[nboundaryGPU+i] = rzNew;



  int icel;
  {
    double invdx = double(mxNeighborsGPU)/lxGPU;
    double invdy = double(myNeighborsGPU)/lyGPU;
    double invdz = double(mzNeighborsGPU)/lzGPU;
    rxNew = rxNew - (int(rxNew*invlxGPU + 0.5*((rxNew>0)-(rxNew<0)))) * lxGPU;
    int jx   = int(rxNew * invdx + 0.5*mxNeighborsGPU) % mxNeighborsGPU;
    ryNew = ryNew - (int(ryNew*invlyGPU + 0.5*((ryNew>0)-(ryNew<0)))) * lyGPU;
    int jy   = int(ryNew * invdy + 0.5*myNeighborsGPU) % myNeighborsGPU;
    rzNew = rzNew - (int(rzNew*invlzGPU + 0.5*((rzNew>0)-(rzNew<0)))) * lzGPU;
    int jz   = int(rzNew * invdz + 0.5*mzNeighborsGPU) % mzNeighborsGPU;
    icel  = jx;
    icel += jy * mxNeighborsGPU;
    icel += jz * mxNeighborsGPU * myNeighborsGPU;    
  }
  int np = atomicAdd(&pc->countPartInCellNonBonded[icel],1);
  if(np >= maxNumberPartInCellNonBondedGPU){
    errorKernel[0] = 1;
    errorKernel[4] = 1;
    return;
  }
  pc->partInCellNonBonded[mNeighborsGPU*np+icel] = i;



}

















__global__ void kernelCorrectionVQuasiNeutrallyBuoyantSemiImplicitTEST4_2(double* rxcellGPU,
									  double* rycellGPU,
									  double* rzcellGPU,
									  cufftDoubleComplex* vxZ,
									  cufftDoubleComplex* vyZ,
									  cufftDoubleComplex* vzZ,
									  double* fxboundaryGPU,
									  double* fyboundaryGPU,
									  double* fzboundaryGPU,
									  double* vxboundaryGPU,
									  double* vyboundaryGPU,
									  double* vzboundaryGPU){
  
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;
  int icelx, icely, icelz;

  double r, rp, rm;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }

  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;

  //CORRECTION IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);

  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx  = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);
                              
  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  double v = dlxm * dlym * dlzm * vxZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vxZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vxZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vxZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vxZ[vecino2].x +
    dlxm * dly  * dlzp * vxZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vxZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vxZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vxZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vxZ[vecinomymz].x +
    dlx  * dlym * dlz  * vxZ[vecino1].x +
    dlx  * dlym * dlzp * vxZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vxZ[vecino0].x + 
    dlx  * dly  * dlz  * vxZ[icelx].x + 
    dlx  * dly  * dlzp * vxZ[vecino5].x + 
    dlx  * dlyp * dlzm * vxZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vxZ[vecino4].x + 
    dlx  * dlyp * dlzp * vxZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vxZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vxZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vxZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vxZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vxZ[vecino3].x +
    dlxp * dly  * dlzp * vxZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vxZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vxZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vxZ[vecinopxpypz].x;

  double fact = massParticleGPU / (volumeGPU * densfluidGPU);
  double deltaV = fact * v ;

  int offset = nboundaryGPU;
  fxboundaryGPU[offset+i]        -= dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fxboundaryGPU[offset+i] -= dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fxboundaryGPU[offset+i] -= dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fxboundaryGPU[offset+i] -= dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fxboundaryGPU[offset+i] -= dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fxboundaryGPU[offset+i] -= dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fxboundaryGPU[offset+i] -= dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fxboundaryGPU[offset+i] -= dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fxboundaryGPU[offset+i] -= dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fxboundaryGPU[offset+i] -= dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fxboundaryGPU[offset+i] -= dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fxboundaryGPU[offset+i] -= dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fxboundaryGPU[offset+i] -= dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fxboundaryGPU[offset+i] -= dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fxboundaryGPU[offset+i] -= dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fxboundaryGPU[offset+i] -= dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fxboundaryGPU[offset+i] -= dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fxboundaryGPU[offset+i] -= dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fxboundaryGPU[offset+i] -= dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fxboundaryGPU[offset+i] -= dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fxboundaryGPU[offset+i] -= dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fxboundaryGPU[offset+i] -= dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fxboundaryGPU[offset+i] -= dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fxboundaryGPU[offset+i] -= dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fxboundaryGPU[offset+i] -= dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fxboundaryGPU[offset+i] -= dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fxboundaryGPU[offset+i] -= dlxm * dlym * dlzm * deltaV;


  //CORRECTION IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  v = dlxm * dlym * dlzm * vyZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vyZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vyZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vyZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vyZ[vecino2].x +
    dlxm * dly  * dlzp * vyZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vyZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vyZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vyZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vyZ[vecinomymz].x +
    dlx  * dlym * dlz  * vyZ[vecino1].x +
    dlx  * dlym * dlzp * vyZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vyZ[vecino0].x + 
    dlx  * dly  * dlz  * vyZ[icely].x + 
    dlx  * dly  * dlzp * vyZ[vecino5].x + 
    dlx  * dlyp * dlzm * vyZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vyZ[vecino4].x + 
    dlx  * dlyp * dlzp * vyZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vyZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vyZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vyZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vyZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vyZ[vecino3].x +
    dlxp * dly  * dlzp * vyZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vyZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vyZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vyZ[vecinopxpypz].x;

  //deltaV = fact * v;
  deltaV = fact * v ;

  offset = nboundaryGPU;
  fyboundaryGPU[offset+i]        -= dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fyboundaryGPU[offset+i] -= dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fyboundaryGPU[offset+i] -= dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fyboundaryGPU[offset+i] -= dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fyboundaryGPU[offset+i] -= dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fyboundaryGPU[offset+i] -= dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fyboundaryGPU[offset+i] -= dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fyboundaryGPU[offset+i] -= dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fyboundaryGPU[offset+i] -= dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fyboundaryGPU[offset+i] -= dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fyboundaryGPU[offset+i] -= dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fyboundaryGPU[offset+i] -= dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fyboundaryGPU[offset+i] -= dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fyboundaryGPU[offset+i] -= dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fyboundaryGPU[offset+i] -= dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fyboundaryGPU[offset+i] -= dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fyboundaryGPU[offset+i] -= dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fyboundaryGPU[offset+i] -= dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fyboundaryGPU[offset+i] -= dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fyboundaryGPU[offset+i] -= dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fyboundaryGPU[offset+i] -= dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fyboundaryGPU[offset+i] -= dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fyboundaryGPU[offset+i] -= dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fyboundaryGPU[offset+i] -= dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fyboundaryGPU[offset+i] -= dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fyboundaryGPU[offset+i] -= dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fyboundaryGPU[offset+i] -= dlxm * dlym * dlzm * deltaV;


  //CORRECTION IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vzZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vzZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vzZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vzZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vzZ[vecino2].x +
    dlxm * dly  * dlzp * vzZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vzZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vzZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vzZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vzZ[vecinomymz].x +
    dlx  * dlym * dlz  * vzZ[vecino1].x +
    dlx  * dlym * dlzp * vzZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vzZ[vecino0].x + 
    dlx  * dly  * dlz  * vzZ[icelz].x + 
    dlx  * dly  * dlzp * vzZ[vecino5].x + 
    dlx  * dlyp * dlzm * vzZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vzZ[vecino4].x + 
    dlx  * dlyp * dlzp * vzZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vzZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vzZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vzZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vzZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vzZ[vecino3].x +
    dlxp * dly  * dlzp * vzZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vzZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vzZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vzZ[vecinopxpypz].x;

  //deltaV = fact * v;
  deltaV = fact * v ;

  offset = nboundaryGPU;
  fzboundaryGPU[offset+i]        -= dlxp * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//1
  fzboundaryGPU[offset+i] -= dlxp * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//2
  fzboundaryGPU[offset+i] -= dlxp * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//3
  fzboundaryGPU[offset+i] -= dlxp * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//4
  fzboundaryGPU[offset+i] -= dlxp * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//5
  fzboundaryGPU[offset+i] -= dlxp * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//6
  fzboundaryGPU[offset+i] -= dlxp * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//7
  fzboundaryGPU[offset+i] -= dlxp * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//8
  fzboundaryGPU[offset+i] -= dlxp * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//9
  fzboundaryGPU[offset+i] -= dlx  * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//10
  fzboundaryGPU[offset+i] -= dlx  * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//11
  fzboundaryGPU[offset+i] -= dlx  * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//12
  fzboundaryGPU[offset+i] -= dlx  * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//13
  fzboundaryGPU[offset+i] -= dlx  * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//14
  fzboundaryGPU[offset+i] -= dlx  * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//15
  fzboundaryGPU[offset+i] -= dlx  * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//16
  fzboundaryGPU[offset+i] -= dlx  * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//17
  fzboundaryGPU[offset+i] -= dlx  * dlym * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//18
  fzboundaryGPU[offset+i] -= dlxm * dlyp * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//19
  fzboundaryGPU[offset+i] -= dlxm * dlyp * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//20
  fzboundaryGPU[offset+i] -= dlxm * dlyp * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//21
  fzboundaryGPU[offset+i] -= dlxm * dly  * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//22
  fzboundaryGPU[offset+i] -= dlxm * dly  * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//23
  fzboundaryGPU[offset+i] -= dlxm * dly  * dlzm * deltaV;
  offset += nboundaryGPU+npGPU;//24
  fzboundaryGPU[offset+i] -= dlxm * dlym * dlzp * deltaV;
  offset += nboundaryGPU+npGPU;//25
  fzboundaryGPU[offset+i] -= dlxm * dlym * dlz  * deltaV;
  offset += nboundaryGPU+npGPU;//26
  fzboundaryGPU[offset+i] -= dlxm * dlym * dlzm * deltaV;


}




























__global__ void updateParticlesQuasiNeutrallyBuoyantSemiImplicitTEST4(particlesincell* pc, 
								      int* errorKernel,
								      double* rxboundaryGPU, 
								      double* ryboundaryGPU, 
								      double* rzboundaryGPU,
								      double* vxboundaryGPU, 
								      double* vyboundaryGPU, 
								      double* vzboundaryGPU,
								      double* rxcellGPU,
								      double* rycellGPU,
								      double* rzcellGPU,
								      double* vxPredictionGPU,
								      double* vyPredictionGPU,
								      double* vzPredictionGPU,
								      cufftDoubleComplex* vxZ,
								      cufftDoubleComplex* vyZ,
								      cufftDoubleComplex* vzZ,
								      double* vxboundaryPredictionGPU,
								      double* vyboundaryPredictionGPU,
								      double* vzboundaryPredictionGPU){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * fetch_double(texVxGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVxGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVxGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVxGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVxGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVxGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVxGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVxGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVxGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVxGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVxGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVxGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVxGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVxGPU,icelx) + 
    dlx  * dly  * dlzp * fetch_double(texVxGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVxGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVxGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVxGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVxGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVxGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVxGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVxGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVxGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVxGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVxGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVxGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVxGPU,vecinopxpypz);
  
  double deltaV = dlxm * dlym * dlzm * vxZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vxZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vxZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vxZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vxZ[vecino2].x +
    dlxm * dly  * dlzp * vxZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vxZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vxZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vxZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vxZ[vecinomymz].x +
    dlx  * dlym * dlz  * vxZ[vecino1].x +
    dlx  * dlym * dlzp * vxZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vxZ[vecino0].x + 
    dlx  * dly  * dlz  * vxZ[icelx].x + 
    dlx  * dly  * dlzp * vxZ[vecino5].x + 
    dlx  * dlyp * dlzm * vxZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vxZ[vecino4].x + 
    dlx  * dlyp * dlzp * vxZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vxZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vxZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vxZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vxZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vxZ[vecino3].x +
    dlxp * dly  * dlzp * vxZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vxZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vxZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vxZ[vecinopxpypz].x;

  
  vxboundaryGPU[nboundaryGPU+i] = v + deltaV + vxboundaryPredictionGPU[nboundaryGPU+i];
  
      
  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  
  v = dlxm * dlym * dlzm * fetch_double(texVyGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVyGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVyGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVyGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVyGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVyGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVyGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVyGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVyGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVyGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVyGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVyGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVyGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVyGPU,icely) + 
    dlx  * dly  * dlzp * fetch_double(texVyGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVyGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVyGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVyGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVyGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVyGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVyGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVyGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVyGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVyGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVyGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVyGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVyGPU,vecinopxpypz);
									  
  deltaV = dlxm * dlym * dlzm * vyZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vyZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vyZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vyZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vyZ[vecino2].x +
    dlxm * dly  * dlzp * vyZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vyZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vyZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vyZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vyZ[vecinomymz].x +
    dlx  * dlym * dlz  * vyZ[vecino1].x +
    dlx  * dlym * dlzp * vyZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vyZ[vecino0].x + 
    dlx  * dly  * dlz  * vyZ[icely].x + 
    dlx  * dly  * dlzp * vyZ[vecino5].x + 
    dlx  * dlyp * dlzm * vyZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vyZ[vecino4].x + 
    dlx  * dlyp * dlzp * vyZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vyZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vyZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vyZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vyZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vyZ[vecino3].x +
    dlxp * dly  * dlzp * vyZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vyZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vyZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vyZ[vecinopxpypz].x;

  
  vyboundaryGPU[nboundaryGPU+i] = v + deltaV + vyboundaryPredictionGPU[nboundaryGPU+i];
  
  
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * fetch_double(texVzGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVzGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVzGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVzGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVzGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVzGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVzGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVzGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVzGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVzGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVzGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVzGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVzGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVzGPU,icelz) + 
    dlx  * dly  * dlzp * fetch_double(texVzGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVzGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVzGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVzGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVzGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVzGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVzGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVzGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVzGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVzGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVzGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVzGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVzGPU,vecinopxpypz);
  
  deltaV = dlxm * dlym * dlzm * vzZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vzZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vzZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vzZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vzZ[vecino2].x +
    dlxm * dly  * dlzp * vzZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vzZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vzZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vzZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vzZ[vecinomymz].x +
    dlx  * dlym * dlz  * vzZ[vecino1].x +
    dlx  * dlym * dlzp * vzZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vzZ[vecino0].x + 
    dlx  * dly  * dlz  * vzZ[icelz].x + 
    dlx  * dly  * dlzp * vzZ[vecino5].x + 
    dlx  * dlyp * dlzm * vzZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vzZ[vecino4].x + 
    dlx  * dlyp * dlzp * vzZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vzZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vzZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vzZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vzZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vzZ[vecino3].x +
    dlxp * dly  * dlzp * vzZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vzZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vzZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vzZ[vecinopxpypz].x;

  
  vzboundaryGPU[nboundaryGPU+i] = v + deltaV + vzboundaryPredictionGPU[nboundaryGPU+i];
  

  

}







__global__ void updateParticlesQuasiNeutrallyBuoyantSemiImplicitTEST4_2(particlesincell* pc, 
									int* errorKernel,
									double* rxboundaryGPU, 
									double* ryboundaryGPU, 
									double* rzboundaryGPU,
									double* vxboundaryGPU, 
									double* vyboundaryGPU, 
									double* vzboundaryGPU,
									double* vxboundaryPredictionGPU,
									double* vyboundaryPredictionGPU,
									double* vzboundaryPredictionGPU,
									double* rxcellGPU,
									double* rycellGPU,
									double* rzcellGPU,
									double* vxPredictionGPU,
									double* vyPredictionGPU,
									double* vzPredictionGPU,
									cufftDoubleComplex* vxZ,
									cufftDoubleComplex* vyZ,
									cufftDoubleComplex* vzZ){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * fetch_double(texVxGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVxGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVxGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVxGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVxGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVxGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVxGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVxGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVxGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVxGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVxGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVxGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVxGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVxGPU,icelx) + 
    dlx  * dly  * dlzp * fetch_double(texVxGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVxGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVxGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVxGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVxGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVxGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVxGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVxGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVxGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVxGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVxGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVxGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVxGPU,vecinopxpypz);
  
  double deltaV = dlxm * dlym * dlzm * vxZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vxZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vxZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vxZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vxZ[vecino2].x +
    dlxm * dly  * dlzp * vxZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vxZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vxZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vxZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vxZ[vecinomymz].x +
    dlx  * dlym * dlz  * vxZ[vecino1].x +
    dlx  * dlym * dlzp * vxZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vxZ[vecino0].x + 
    dlx  * dly  * dlz  * vxZ[icelx].x + 
    dlx  * dly  * dlzp * vxZ[vecino5].x + 
    dlx  * dlyp * dlzm * vxZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vxZ[vecino4].x + 
    dlx  * dlyp * dlzp * vxZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vxZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vxZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vxZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vxZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vxZ[vecino3].x +
    dlxp * dly  * dlzp * vxZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vxZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vxZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vxZ[vecinopxpypz].x;


  vxboundaryPredictionGPU[nboundaryGPU+i] = 0.5 * (vxboundaryGPU[nboundaryGPU+i] + v + deltaV);
  
      
  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

  
  v = dlxm * dlym * dlzm * fetch_double(texVyGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVyGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVyGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVyGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVyGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVyGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVyGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVyGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVyGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVyGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVyGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVyGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVyGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVyGPU,icely) + 
    dlx  * dly  * dlzp * fetch_double(texVyGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVyGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVyGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVyGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVyGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVyGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVyGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVyGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVyGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVyGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVyGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVyGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVyGPU,vecinopxpypz);
									  
  deltaV = dlxm * dlym * dlzm * vyZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vyZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vyZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vyZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vyZ[vecino2].x +
    dlxm * dly  * dlzp * vyZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vyZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vyZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vyZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vyZ[vecinomymz].x +
    dlx  * dlym * dlz  * vyZ[vecino1].x +
    dlx  * dlym * dlzp * vyZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vyZ[vecino0].x + 
    dlx  * dly  * dlz  * vyZ[icely].x + 
    dlx  * dly  * dlzp * vyZ[vecino5].x + 
    dlx  * dlyp * dlzm * vyZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vyZ[vecino4].x + 
    dlx  * dlyp * dlzp * vyZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vyZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vyZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vyZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vyZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vyZ[vecino3].x +
    dlxp * dly  * dlzp * vyZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vyZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vyZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vyZ[vecinopxpypz].x;

									  
  vyboundaryPredictionGPU[nboundaryGPU+i] = 0.5 * (vyboundaryGPU[nboundaryGPU+i] + v + deltaV);
  
  
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  //dlx = tex1D(texDelta, fabs(r));
  //dlxp = tex1D(texDelta, fabs(rp));
  //dlxm = tex1D(texDelta, fabs(rm));
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  //dly = tex1D(texDelta, fabs(r));
  //dlyp = tex1D(texDelta, fabs(rp));
  //dlym = tex1D(texDelta, fabs(rm));
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  //dlz = tex1D(texDelta, fabs(r));
  //dlzp = tex1D(texDelta, fabs(rp));
  //dlzm = tex1D(texDelta, fabs(rm));
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * fetch_double(texVzGPU,vecinomxmymz) +
    dlxm * dlym * dlz  * fetch_double(texVzGPU,vecinomxmy) +
    dlxm * dlym * dlzp * fetch_double(texVzGPU,vecinomxmypz) +
    dlxm * dly  * dlzm * fetch_double(texVzGPU,vecinomxmz) + 
    dlxm * dly  * dlz  * fetch_double(texVzGPU,vecino2) +
    dlxm * dly  * dlzp * fetch_double(texVzGPU,vecinomxpz) +
    dlxm * dlyp * dlzm * fetch_double(texVzGPU,vecinomxpymz) +
    dlxm * dlyp * dlz  * fetch_double(texVzGPU,vecinomxpy) +
    dlxm * dlyp * dlzp * fetch_double(texVzGPU,vecinomxpypz) +
    dlx  * dlym * dlzm * fetch_double(texVzGPU,vecinomymz) +
    dlx  * dlym * dlz  * fetch_double(texVzGPU,vecino1) +
    dlx  * dlym * dlzp * fetch_double(texVzGPU,vecinomypz) + 
    dlx  * dly  * dlzm * fetch_double(texVzGPU,vecino0) + 
    dlx  * dly  * dlz  * fetch_double(texVzGPU,icelz) + 
    dlx  * dly  * dlzp * fetch_double(texVzGPU,vecino5) + 
    dlx  * dlyp * dlzm * fetch_double(texVzGPU,vecinopymz) + 
    dlx  * dlyp * dlz  * fetch_double(texVzGPU,vecino4) + 
    dlx  * dlyp * dlzp * fetch_double(texVzGPU,vecinopypz) + 
    dlxp * dlym * dlzm * fetch_double(texVzGPU,vecinopxmymz) + 
    dlxp * dlym * dlz  * fetch_double(texVzGPU,vecinopxmy) + 
    dlxp * dlym * dlzp * fetch_double(texVzGPU,vecinopxmypz) +
    dlxp * dly  * dlzm * fetch_double(texVzGPU,vecinopxmz) + 
    dlxp * dly  * dlz  * fetch_double(texVzGPU,vecino3) +
    dlxp * dly  * dlzp * fetch_double(texVzGPU,vecinopxpz) +
    dlxp * dlyp * dlzm * fetch_double(texVzGPU,vecinopxpymz) +
    dlxp * dlyp * dlz  * fetch_double(texVzGPU,vecinopxpy) +
    dlxp * dlyp * dlzp * fetch_double(texVzGPU,vecinopxpypz);
  
  deltaV = dlxm * dlym * dlzm * vzZ[vecinomxmymz].x +
    dlxm * dlym * dlz  * vzZ[vecinomxmy].x +
    dlxm * dlym * dlzp * vzZ[vecinomxmypz].x +
    dlxm * dly  * dlzm * vzZ[vecinomxmz].x + 
    dlxm * dly  * dlz  * vzZ[vecino2].x +
    dlxm * dly  * dlzp * vzZ[vecinomxpz].x +
    dlxm * dlyp * dlzm * vzZ[vecinomxpymz].x +
    dlxm * dlyp * dlz  * vzZ[vecinomxpy].x +
    dlxm * dlyp * dlzp * vzZ[vecinomxpypz].x +
    dlx  * dlym * dlzm * vzZ[vecinomymz].x +
    dlx  * dlym * dlz  * vzZ[vecino1].x +
    dlx  * dlym * dlzp * vzZ[vecinomypz].x + 
    dlx  * dly  * dlzm * vzZ[vecino0].x + 
    dlx  * dly  * dlz  * vzZ[icelz].x + 
    dlx  * dly  * dlzp * vzZ[vecino5].x + 
    dlx  * dlyp * dlzm * vzZ[vecinopymz].x + 
    dlx  * dlyp * dlz  * vzZ[vecino4].x + 
    dlx  * dlyp * dlzp * vzZ[vecinopypz].x + 
    dlxp * dlym * dlzm * vzZ[vecinopxmymz].x + 
    dlxp * dlym * dlz  * vzZ[vecinopxmy].x + 
    dlxp * dlym * dlzp * vzZ[vecinopxmypz].x +
    dlxp * dly  * dlzm * vzZ[vecinopxmz].x + 
    dlxp * dly  * dlz  * vzZ[vecino3].x +
    dlxp * dly  * dlzp * vzZ[vecinopxpz].x +
    dlxp * dlyp * dlzm * vzZ[vecinopxpymz].x +
    dlxp * dlyp * dlz  * vzZ[vecinopxpy].x +
    dlxp * dlyp * dlzp * vzZ[vecinopxpypz].x;

  
  vzboundaryPredictionGPU[nboundaryGPU+i] = 0.5 * (vzboundaryGPU[nboundaryGPU+i] + v + deltaV);
  
  

}

















__global__ void findNeighborParticlesQuasiNeutrallyBuoyantTEST4_2(particlesincell* pc, 
								  int* errorKernel,
								  const double* rxcellGPU,
								  const double* rycellGPU,
								  const double* rzcellGPU,
								  double* rxboundaryGPU,  //q^{n}
								  double* ryboundaryGPU, 
								  double* rzboundaryGPU,
								  const double* rxboundaryPredictionGPU, //q^{n+1/2}
								  const double* ryboundaryPredictionGPU, 
								  const double* rzboundaryPredictionGPU,
								  const double* vxGPU, //v^{n+1/2}
								  const double* vyGPU, 
								  const double* vzGPU){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i); //rxboundaryPredictionGPU
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * mytGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];

  
  rxboundaryGPU[i] = rxboundaryGPU[i] + dtGPU * v;
  //r = rxboundaryGPU[i] + dtGPU * v;
  //rxboundaryGPU[i] =  r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU;


  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  

  ryboundaryGPU[i] = ryboundaryGPU[i] + dtGPU * v;
  //r = ryboundaryGPU[i] + dtGPU * v;
  //ryboundaryGPU[i] =  r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU;
 
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];


  rzboundaryGPU[i] = rzboundaryGPU[i] + dtGPU * v;
  //r = rzboundaryGPU[i] + dtGPU * v;
  //rzboundaryGPU[i] =  r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU;
  

}







//For the first time step
//Calculate positions at n+1/2
//and save in rxboundaryPrediction
__global__ void findNeighborParticlesQuasiNeutrallyBuoyantTEST4_3(particlesincell* pc, 
								  int* errorKernel,
								  double* rxcellGPU,
								  double* rycellGPU,
								  double* rzcellGPU,
								  double* rxboundaryGPU,  //q^{n}
								  double* ryboundaryGPU, 
								  double* rzboundaryGPU,
								  double* rxboundaryPredictionGPU, //q^{n+1/2}
								  double* ryboundaryPredictionGPU, 
								  double* rzboundaryPredictionGPU,
								  double* vxGPU, //v^{n+1/2}
								  double* vyGPU, 
								  double* vzGPU){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   

  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i); //rxboundaryPredictionGPU
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];

  
  rxboundaryPredictionGPU[i] = rxboundaryGPU[i] + 0.5 * dtGPU * v;


  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  
  
  ryboundaryPredictionGPU[i] = ryboundaryGPU[i] + 0.5 * dtGPU * v;
  
 
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];

  
  rzboundaryPredictionGPU[i] = rzboundaryGPU[i] + 0.5 * dtGPU * v;
  


}










__global__ void calculateVelocityAtHalfTimeStep(double* vxGPU,
						double* vyGPU,
						double* vzGPU,
						double* vxPredictionGPU,
						double* vyPredictionGPU,
						double* vzPredictionGPU,
						cufftDoubleComplex* vxZ,
						cufftDoubleComplex* vyZ,
						cufftDoubleComplex* vzZ){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(ncellsGPU)) return;   

  vxGPU[i] = 0.5 * (vxGPU[i] + vxPredictionGPU[i] + vxZ[i].x / double(ncellsGPU));
  vyGPU[i] = 0.5 * (vyGPU[i] + vyPredictionGPU[i] + vyZ[i].x / double(ncellsGPU));
  vzGPU[i] = 0.5 * (vzGPU[i] + vzPredictionGPU[i] + vzZ[i].x / double(ncellsGPU));

}















//Calculate nu*L*\Delta v^{k=2} and store it
//in vxGPU
__global__ void laplacianDeltaV(const cufftDoubleComplex* vxZ,
				const cufftDoubleComplex* vyZ,
				const cufftDoubleComplex* vzZ,
				double* vxGPU,
				double* vyGPU,
				double* vzGPU){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(ncellsGPU)) return;   

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5; 
  
  double wx, wy, wz;
  double vx, vy, vz;
  double vx0, vx1, vx2, vx3, vx4, vx5;
  double vy0, vy1, vy2, vy3, vy4, vy5;
  double vz0, vz1, vz2, vz3, vz4, vz5;

  vecino0 = tex1Dfetch(texvecino0GPU,i);
  vecino1 = tex1Dfetch(texvecino1GPU,i);
  vecino2 = tex1Dfetch(texvecino2GPU,i);
  vecino3 = tex1Dfetch(texvecino3GPU,i);
  vecino4 = tex1Dfetch(texvecino4GPU,i);
  vecino5 = tex1Dfetch(texvecino5GPU,i);

  vx = vxZ[i].x / double(ncellsGPU);
  vy = vyZ[i].x / double(ncellsGPU);
  vz = vzZ[i].x / double(ncellsGPU);
  vx0 = vxZ[vecino0].x / double(ncellsGPU);
  vx1 = vxZ[vecino1].x / double(ncellsGPU);
  vx2 = vxZ[vecino2].x / double(ncellsGPU);
  vx3 = vxZ[vecino3].x / double(ncellsGPU);
  vx4 = vxZ[vecino4].x / double(ncellsGPU);
  vx5 = vxZ[vecino5].x / double(ncellsGPU);
  vy0 = vyZ[vecino0].x / double(ncellsGPU);
  vy1 = vyZ[vecino1].x / double(ncellsGPU);
  vy2 = vyZ[vecino2].x / double(ncellsGPU);
  vy3 = vyZ[vecino3].x / double(ncellsGPU);
  vy4 = vyZ[vecino4].x / double(ncellsGPU);
  vy5 = vyZ[vecino5].x / double(ncellsGPU);
  vz0 = vzZ[vecino0].x / double(ncellsGPU);
  vz1 = vzZ[vecino1].x / double(ncellsGPU);
  vz2 = vzZ[vecino2].x / double(ncellsGPU);
  vz3 = vzZ[vecino3].x / double(ncellsGPU);
  vz4 = vzZ[vecino4].x / double(ncellsGPU);
  vz5 = vzZ[vecino5].x / double(ncellsGPU);



  //Laplacian part
  wx  = invdxGPU * invdxGPU * (vx3 - 2*vx + vx2);
  wx += invdyGPU * invdyGPU * (vx4 - 2*vx + vx1);
  wx += invdzGPU * invdzGPU * (vx5 - 2*vx + vx0);
  wx  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wx;
  wy  = invdxGPU * invdxGPU * (vy3 - 2*vy + vy2);
  wy += invdyGPU * invdyGPU * (vy4 - 2*vy + vy1);
  wy += invdzGPU * invdzGPU * (vy5 - 2*vy + vy0);
  wy  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wy;
  wz  = invdxGPU * invdxGPU * (vz3 - 2*vz + vz2);
  wz += invdyGPU * invdyGPU * (vz4 - 2*vz + vz1);
  wz += invdzGPU * invdzGPU * (vz5 - 2*vz + vz0);
  wz  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wz;


  vxGPU[i] = wx;
  vyGPU[i] = wy;
  vzGPU[i] = wz;

}























//Calculate 0.5*dt*nu*J*L*\Delta v^{k=2} and store it
//in vxboundaryPredictionGPU
__global__ void interpolateLaplacianDeltaV(double* rxcellGPU,
					   double* rycellGPU,
					   double* rzcellGPU,
					   double* vxGPU,
					   double* vyGPU,
					   double* vzGPU,
					   double* rxboundaryPredictionGPU,
					   double* ryboundaryPredictionGPU,
					   double* rzboundaryPredictionGPU,
					   double* vxboundaryPredictionGPU,
					   double* vyboundaryPredictionGPU,
					   double* vzboundaryPredictionGPU){
  

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;     



  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i); //rxboundaryPredictionGPU
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vxGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vxGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vxGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vxGPU[vecino2] +
    dlxm * dly  * dlzp * vxGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vxGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vxGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vxGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vxGPU[vecinomymz] +
    dlx  * dlym * dlz  * vxGPU[vecino1] +
    dlx  * dlym * dlzp * vxGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vxGPU[vecino0] + 
    dlx  * dly  * dlz  * vxGPU[icelx] + 
    dlx  * dly  * dlzp * vxGPU[vecino5] + 
    dlx  * dlyp * dlzm * vxGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vxGPU[vecino4] + 
    dlx  * dlyp * dlzp * vxGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vxGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vxGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vxGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vxGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vxGPU[vecino3] +
    dlxp * dly  * dlzp * vxGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vxGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vxGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vxGPU[vecinopxpypz];

  
  vxboundaryPredictionGPU[i] = v;

  

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vyGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vyGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vyGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vyGPU[vecino2] +
    dlxm * dly  * dlzp * vyGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vyGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vyGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vyGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vyGPU[vecinomymz] +
    dlx  * dlym * dlz  * vyGPU[vecino1] +
    dlx  * dlym * dlzp * vyGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vyGPU[vecino0] + 
    dlx  * dly  * dlz  * vyGPU[icely] + 
    dlx  * dly  * dlzp * vyGPU[vecino5] + 
    dlx  * dlyp * dlzm * vyGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vyGPU[vecino4] + 
    dlx  * dlyp * dlzp * vyGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vyGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vyGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vyGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vyGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vyGPU[vecino3] +
    dlxp * dly  * dlzp * vyGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vyGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vyGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vyGPU[vecinopxpypz];
  

  vyboundaryPredictionGPU[i] = v;  
 
 
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzGPU[vecinomxmymz] +
    dlxm * dlym * dlz  * vzGPU[vecinomxmy] +
    dlxm * dlym * dlzp * vzGPU[vecinomxmypz] +
    dlxm * dly  * dlzm * vzGPU[vecinomxmz] + 
    dlxm * dly  * dlz  * vzGPU[vecino2] +
    dlxm * dly  * dlzp * vzGPU[vecinomxpz] +
    dlxm * dlyp * dlzm * vzGPU[vecinomxpymz] +
    dlxm * dlyp * dlz  * vzGPU[vecinomxpy] +
    dlxm * dlyp * dlzp * vzGPU[vecinomxpypz] +
    dlx  * dlym * dlzm * vzGPU[vecinomymz] +
    dlx  * dlym * dlz  * vzGPU[vecino1] +
    dlx  * dlym * dlzp * vzGPU[vecinomypz] + 
    dlx  * dly  * dlzm * vzGPU[vecino0] + 
    dlx  * dly  * dlz  * vzGPU[icelz] + 
    dlx  * dly  * dlzp * vzGPU[vecino5] + 
    dlx  * dlyp * dlzm * vzGPU[vecinopymz] + 
    dlx  * dlyp * dlz  * vzGPU[vecino4] + 
    dlx  * dlyp * dlzp * vzGPU[vecinopypz] + 
    dlxp * dlym * dlzm * vzGPU[vecinopxmymz] + 
    dlxp * dlym * dlz  * vzGPU[vecinopxmy] + 
    dlxp * dlym * dlzp * vzGPU[vecinopxmypz] +
    dlxp * dly  * dlzm * vzGPU[vecinopxmz] + 
    dlxp * dly  * dlz  * vzGPU[vecino3] +
    dlxp * dly  * dlzp * vzGPU[vecinopxpz] +
    dlxp * dlyp * dlzm * vzGPU[vecinopxpymz] +
    dlxp * dlyp * dlz  * vzGPU[vecinopxpy] +
    dlxp * dlyp * dlzp * vzGPU[vecinopxpypz];

  
  vzboundaryPredictionGPU[i] = v;




}
























//Calculate nu*L*\Delta v^{k=2} and store it
//in vxGPU
__global__ void laplacianDeltaV_2(cufftDoubleComplex* vxZ,
				  cufftDoubleComplex* vyZ,
				  cufftDoubleComplex* vzZ){
  
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(ncellsGPU)) return;   

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5; 
  
  double wx, wy, wz;
  double vx, vy, vz;
  double vx0, vx1, vx2, vx3, vx4, vx5;
  double vy0, vy1, vy2, vy3, vy4, vy5;
  double vz0, vz1, vz2, vz3, vz4, vz5;

  vecino0 = tex1Dfetch(texvecino0GPU,i);
  vecino1 = tex1Dfetch(texvecino1GPU,i);
  vecino2 = tex1Dfetch(texvecino2GPU,i);
  vecino3 = tex1Dfetch(texvecino3GPU,i);
  vecino4 = tex1Dfetch(texvecino4GPU,i);
  vecino5 = tex1Dfetch(texvecino5GPU,i);

  vx = vxZ[i].x / double(ncellsGPU);
  vy = vyZ[i].x / double(ncellsGPU);
  vz = vzZ[i].x / double(ncellsGPU);
  vx0 = vxZ[vecino0].x / double(ncellsGPU);
  vx1 = vxZ[vecino1].x / double(ncellsGPU);
  vx2 = vxZ[vecino2].x / double(ncellsGPU);
  vx3 = vxZ[vecino3].x / double(ncellsGPU);
  vx4 = vxZ[vecino4].x / double(ncellsGPU);
  vx5 = vxZ[vecino5].x / double(ncellsGPU);
  vy0 = vyZ[vecino0].x / double(ncellsGPU);
  vy1 = vyZ[vecino1].x / double(ncellsGPU);
  vy2 = vyZ[vecino2].x / double(ncellsGPU);
  vy3 = vyZ[vecino3].x / double(ncellsGPU);
  vy4 = vyZ[vecino4].x / double(ncellsGPU);
  vy5 = vyZ[vecino5].x / double(ncellsGPU);
  vz0 = vzZ[vecino0].x / double(ncellsGPU);
  vz1 = vzZ[vecino1].x / double(ncellsGPU);
  vz2 = vzZ[vecino2].x / double(ncellsGPU);
  vz3 = vzZ[vecino3].x / double(ncellsGPU);
  vz4 = vzZ[vecino4].x / double(ncellsGPU);
  vz5 = vzZ[vecino5].x / double(ncellsGPU);



  //Laplacian part
  wx  = invdxGPU * invdxGPU * (vx3 - 2*vx + vx2);
  wx += invdyGPU * invdyGPU * (vx4 - 2*vx + vx1);
  wx += invdzGPU * invdzGPU * (vx5 - 2*vx + vx0);
  wx  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wx;
  wy  = invdxGPU * invdxGPU * (vy3 - 2*vy + vy2);
  wy += invdyGPU * invdyGPU * (vy4 - 2*vy + vy1);
  wy += invdzGPU * invdzGPU * (vy5 - 2*vy + vy0);
  wy  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wy;
  wz  = invdxGPU * invdxGPU * (vz3 - 2*vz + vz2);
  wz += invdyGPU * invdyGPU * (vz4 - 2*vz + vz1);
  wz += invdzGPU * invdzGPU * (vz5 - 2*vz + vz0);
  wz  = 0.5 * dtGPU * (shearviscosityGPU/densfluidGPU) * wz;


  vxZ[i].y = wx;
  vyZ[i].y = wy;
  vzZ[i].y = wz;

}
















//Calculate 0.5*dt*nu*J*L*\Delta v^{k=2} and store it
//in vxboundaryPredictionGPU
__global__ void interpolateLaplacianDeltaV_2(double* rxcellGPU,
					     double* rycellGPU,
					     double* rzcellGPU,
					     cufftDoubleComplex* vxZ,
					     cufftDoubleComplex* vyZ,
					     cufftDoubleComplex* vzZ,
					     double* rxboundaryPredictionGPU,
					     double* ryboundaryPredictionGPU,
					     double* rzboundaryPredictionGPU,
					     double* vxboundaryPredictionGPU,
					     double* vyboundaryPredictionGPU,
					     double* vzboundaryPredictionGPU){
  

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;     



  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i); //rxboundaryPredictionGPU
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);
  
  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  int vecinopxpy, vecinopxmy, vecinopxpz, vecinopxmz;
  int vecinomxpy, vecinomxmy, vecinomxpz, vecinomxmz;
  int vecinopypz, vecinopymz, vecinomypz, vecinomymz;
  int vecinopxpypz, vecinopxpymz, vecinopxmypz, vecinopxmymz;
  int vecinomxpypz, vecinomxpymz, vecinomxmypz, vecinomxmymz;

  double r, rp, rm;
  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;
  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;

  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * myGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*myGPU) % myGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*myGPU) % myGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }
  
  //VELOCITY IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelx);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelx);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelx);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelx);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelx);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelx);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelx);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelx);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelx);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelx);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelx);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelx);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelx);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelx);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelx);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelx);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelx);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelx);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelx);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelx);
  //int vecinopxpxpypz = tex1Dfetch(texvecino3GPU, vecinopxpypz);
  //int vecinopxpxpymz = tex1Dfetch(texvecino3GPU, vecinopxpymz);
  //int vecinopxpxmypz = tex1Dfetch(texvecino3GPU, vecinopxmypz);
  //int vecinopxpxmymz = tex1Dfetch(texvecino3GPU, vecinopxmymz);
  //int vecinopxpx     = tex1Dfetch(texvecino3GPU, vecino3);
  //int vecinopxpxpy   = tex1Dfetch(texvecino3GPU, vecinopxpy);
  //int vecinopxpxmy   = tex1Dfetch(texvecino3GPU, vecinopxmy);
  //int vecinopxpxpz   = tex1Dfetch(texvecino3GPU, vecinopxpz);
  //int vecinopxpxmz   = tex1Dfetch(texvecino3GPU, vecinopxmz);
  
  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);
  
  double v = dlxm * dlym * dlzm * vxZ[vecinomxmymz].y +
    dlxm * dlym * dlz  * vxZ[vecinomxmy].y +
    dlxm * dlym * dlzp * vxZ[vecinomxmypz].y +
    dlxm * dly  * dlzm * vxZ[vecinomxmz].y + 
    dlxm * dly  * dlz  * vxZ[vecino2].y +
    dlxm * dly  * dlzp * vxZ[vecinomxpz].y +
    dlxm * dlyp * dlzm * vxZ[vecinomxpymz].y +
    dlxm * dlyp * dlz  * vxZ[vecinomxpy].y +
    dlxm * dlyp * dlzp * vxZ[vecinomxpypz].y +
    dlx  * dlym * dlzm * vxZ[vecinomymz].y +
    dlx  * dlym * dlz  * vxZ[vecino1].y +
    dlx  * dlym * dlzp * vxZ[vecinomypz].y + 
    dlx  * dly  * dlzm * vxZ[vecino0].y + 
    dlx  * dly  * dlz  * vxZ[icelx].y + 
    dlx  * dly  * dlzp * vxZ[vecino5].y + 
    dlx  * dlyp * dlzm * vxZ[vecinopymz].y + 
    dlx  * dlyp * dlz  * vxZ[vecino4].y + 
    dlx  * dlyp * dlzp * vxZ[vecinopypz].y + 
    dlxp * dlym * dlzm * vxZ[vecinopxmymz].y + 
    dlxp * dlym * dlz  * vxZ[vecinopxmy].y + 
    dlxp * dlym * dlzp * vxZ[vecinopxmypz].y +
    dlxp * dly  * dlzm * vxZ[vecinopxmz].y + 
    dlxp * dly  * dlz  * vxZ[vecino3].y +
    dlxp * dly  * dlzp * vxZ[vecinopxpz].y +
    dlxp * dlyp * dlzm * vxZ[vecinopxpymz].y +
    dlxp * dlyp * dlz  * vxZ[vecinopxpy].y +
    dlxp * dlyp * dlzp * vxZ[vecinopxpypz].y;

  
  vxboundaryPredictionGPU[i] = v;

  

  //VELOCITY IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icely);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icely);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icely);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icely);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icely);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icely);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icely);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icely);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icely);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icely);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icely);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icely);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icely);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icely);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icely);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icely);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icely);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icely);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icely);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icely);  
  //DEFINE MORE NEIGHBORS
  int vecinopymxpymz = tex1Dfetch(texvecino4GPU, vecinomxpymz);
  int vecinopymxpy   = tex1Dfetch(texvecino4GPU, vecinomxpy);
  int vecinopymxpypz = tex1Dfetch(texvecino4GPU, vecinomxpypz);
  int vecinopypymz   = tex1Dfetch(texvecino4GPU, vecinopymz);
  int vecinopypy     = tex1Dfetch(texvecino4GPU, vecino4);
  int vecinopypypz   = tex1Dfetch(texvecino4GPU, vecinopypz);
  int vecinopypxpymz = tex1Dfetch(texvecino4GPU, vecinopxpymz);
  int vecinopypxpy   = tex1Dfetch(texvecino4GPU, vecinopxpy);
  int vecinopypxpypz = tex1Dfetch(texvecino4GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);


  v = dlxm * dlym * dlzm * vyZ[vecinomxmymz].y +
    dlxm * dlym * dlz  * vyZ[vecinomxmy].y +
    dlxm * dlym * dlzp * vyZ[vecinomxmypz].y +
    dlxm * dly  * dlzm * vyZ[vecinomxmz].y + 
    dlxm * dly  * dlz  * vyZ[vecino2].y +
    dlxm * dly  * dlzp * vyZ[vecinomxpz].y +
    dlxm * dlyp * dlzm * vyZ[vecinomxpymz].y +
    dlxm * dlyp * dlz  * vyZ[vecinomxpy].y +
    dlxm * dlyp * dlzp * vyZ[vecinomxpypz].y +
    dlx  * dlym * dlzm * vyZ[vecinomymz].y +
    dlx  * dlym * dlz  * vyZ[vecino1].y +
    dlx  * dlym * dlzp * vyZ[vecinomypz].y + 
    dlx  * dly  * dlzm * vyZ[vecino0].y + 
    dlx  * dly  * dlz  * vyZ[icely].y + 
    dlx  * dly  * dlzp * vyZ[vecino5].y + 
    dlx  * dlyp * dlzm * vyZ[vecinopymz].y + 
    dlx  * dlyp * dlz  * vyZ[vecino4].y + 
    dlx  * dlyp * dlzp * vyZ[vecinopypz].y + 
    dlxp * dlym * dlzm * vyZ[vecinopxmymz].y + 
    dlxp * dlym * dlz  * vyZ[vecinopxmy].y + 
    dlxp * dlym * dlzp * vyZ[vecinopxmypz].y +
    dlxp * dly  * dlzm * vyZ[vecinopxmz].y + 
    dlxp * dly  * dlz  * vyZ[vecino3].y +
    dlxp * dly  * dlzp * vyZ[vecinopxpz].y +
    dlxp * dlyp * dlzm * vyZ[vecinopxpymz].y +
    dlxp * dlyp * dlz  * vyZ[vecinopxpy].y +
    dlxp * dlyp * dlzp * vyZ[vecinopxpypz].y;
  

  vyboundaryPredictionGPU[i] = v;  
 
 
  //VELOCITY IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
  vecinopxpy = tex1Dfetch(texvecinopxpyGPU, icelz);
  vecinopxmy = tex1Dfetch(texvecinopxmyGPU, icelz);
  vecinopxpz = tex1Dfetch(texvecinopxpzGPU, icelz);
  vecinopxmz = tex1Dfetch(texvecinopxmzGPU, icelz);
  vecinomxpy = tex1Dfetch(texvecinomxpyGPU, icelz);
  vecinomxmy = tex1Dfetch(texvecinomxmyGPU, icelz);
  vecinomxpz = tex1Dfetch(texvecinomxpzGPU, icelz);
  vecinomxmz = tex1Dfetch(texvecinomxmzGPU, icelz);
  vecinopypz = tex1Dfetch(texvecinopypzGPU, icelz);
  vecinopymz = tex1Dfetch(texvecinopymzGPU, icelz);
  vecinomypz = tex1Dfetch(texvecinomypzGPU, icelz);
  vecinomymz = tex1Dfetch(texvecinomymzGPU, icelz);
  vecinopxpypz = tex1Dfetch(texvecinopxpypzGPU, icelz);
  vecinopxpymz = tex1Dfetch(texvecinopxpymzGPU, icelz);
  vecinopxmypz = tex1Dfetch(texvecinopxmypzGPU, icelz);
  vecinopxmymz = tex1Dfetch(texvecinopxmymzGPU, icelz);
  vecinomxpypz = tex1Dfetch(texvecinomxpypzGPU, icelz);
  vecinomxpymz = tex1Dfetch(texvecinomxpymzGPU, icelz);
  vecinomxmypz = tex1Dfetch(texvecinomxmypzGPU, icelz);
  vecinomxmymz = tex1Dfetch(texvecinomxmymzGPU, icelz);  
  //DEFINE MORE NEIGHBORS
  int vecinopzmxmypz = tex1Dfetch(texvecino5GPU, vecinomxmypz);
  int vecinopzmxpz   = tex1Dfetch(texvecino5GPU, vecinomxpz);
  int vecinopzmxpypz = tex1Dfetch(texvecino5GPU, vecinomxpypz);
  int vecinopzmypz   = tex1Dfetch(texvecino5GPU, vecinomypz);
  int vecinopzpz     = tex1Dfetch(texvecino5GPU, vecino5);
  int vecinopzpypz   = tex1Dfetch(texvecino5GPU, vecinopypz);
  int vecinopzpxmypz = tex1Dfetch(texvecino5GPU, vecinopxmypz);
  int vecinopzpxpz   = tex1Dfetch(texvecino5GPU, vecinopxpz);
  int vecinopzpxpypz = tex1Dfetch(texvecino5GPU, vecinopxpypz);

  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = delta(1.5*r);
  dlxp = delta(1.5*rp);
  dlxm = delta(1.5*rm);

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = delta(1.5*r);
  dlyp = delta(1.5*rp);
  dlym = delta(1.5*rm);

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = delta(1.5*r);
  dlzp = delta(1.5*rp);
  dlzm = delta(1.5*rm);

 
  v = dlxm * dlym * dlzm * vzZ[vecinomxmymz].y +
    dlxm * dlym * dlz  * vzZ[vecinomxmy].y +
    dlxm * dlym * dlzp * vzZ[vecinomxmypz].y +
    dlxm * dly  * dlzm * vzZ[vecinomxmz].y + 
    dlxm * dly  * dlz  * vzZ[vecino2].y +
    dlxm * dly  * dlzp * vzZ[vecinomxpz].y +
    dlxm * dlyp * dlzm * vzZ[vecinomxpymz].y +
    dlxm * dlyp * dlz  * vzZ[vecinomxpy].y +
    dlxm * dlyp * dlzp * vzZ[vecinomxpypz].y +
    dlx  * dlym * dlzm * vzZ[vecinomymz].y +
    dlx  * dlym * dlz  * vzZ[vecino1].y +
    dlx  * dlym * dlzp * vzZ[vecinomypz].y + 
    dlx  * dly  * dlzm * vzZ[vecino0].y + 
    dlx  * dly  * dlz  * vzZ[icelz].y + 
    dlx  * dly  * dlzp * vzZ[vecino5].y + 
    dlx  * dlyp * dlzm * vzZ[vecinopymz].y + 
    dlx  * dlyp * dlz  * vzZ[vecino4].y + 
    dlx  * dlyp * dlzp * vzZ[vecinopypz].y + 
    dlxp * dlym * dlzm * vzZ[vecinopxmymz].y + 
    dlxp * dlym * dlz  * vzZ[vecinopxmy].y + 
    dlxp * dlym * dlzp * vzZ[vecinopxmypz].y +
    dlxp * dly  * dlzm * vzZ[vecinopxmz].y + 
    dlxp * dly  * dlz  * vzZ[vecino3].y +
    dlxp * dly  * dlzp * vzZ[vecinopxpz].y +
    dlxp * dlyp * dlzm * vzZ[vecinopxpymz].y +
    dlxp * dlyp * dlz  * vzZ[vecinopxpy].y +
    dlxp * dlyp * dlzp * vzZ[vecinopxpypz].y;

  
  vzboundaryPredictionGPU[i] = v;




}



















// Spread the particles thermal drift
// drift = (m_f_tilde/(m_e+m_f_tilde) * kT * ( \delta S(q) / \delta q )
__global__ void kernelSpreadParticlesThermalDrift(const double *rxcellGPU,
						  const double *rycellGPU,
						  const double *rzcellGPU,
						  double *fxboundaryGPU,
						  double *fyboundaryGPU,
						  double *fzboundaryGPU,
						  particlesincell *pc,
						  int *errorKernel){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>=(npGPU)) return;   
  
  double rx = fetch_double(texrxboundaryGPU,nboundaryGPU+i);
  double ry = fetch_double(texryboundaryGPU,nboundaryGPU+i);
  double rz = fetch_double(texrzboundaryGPU,nboundaryGPU+i);

  int vecino0, vecino1, vecino2, vecino3, vecino4, vecino5;
  
  double r, rp, rm;

  double dlx, dlxp, dlxm;
  double dly, dlyp, dlym;
  double dlz, dlzp, dlzm;
  int icelx, icely, icelz;

  {
    int mxmy = mxGPU * mytGPU;
    r = rx;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jx   = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;
    r = rx - 0.5*dxGPU;
    r = r - (int(r*invlxGPU + 0.5*((r>0)-(r<0)))) * lxGPU;
    int jxdx = int(r * invdxGPU + 0.5*mxGPU) % mxGPU;

    r = ry;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jy   = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;
    r = ry - 0.5*dyGPU;
    r = r - (int(r*invlyGPU + 0.5*((r>0)-(r<0)))) * lyGPU;
    int jydy = int(r * invdyGPU + 0.5*mytGPU) % mytGPU;

    r = rz;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jz   = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;
    r = rz - 0.5*dzGPU;
    r = r - (int(r*invlzGPU + 0.5*((r>0)-(r<0)))) * lzGPU;
    int jzdz = int(r * invdzGPU + 0.5*mzGPU) % mzGPU;

    icelx  = jxdx;
    icelx += jy * mxGPU;
    icelx += jz * mxmy;

    icely  = jx;
    icely += jydy * mxGPU;
    icely += jz * mxmy;

    icelz  = jx;
    icelz += jy * mxGPU;
    icelz += jzdz * mxmy;
  }

  double auxdx = invdxGPU/1.5;
  double auxdy = invdyGPU/1.5;
  double auxdz = invdzGPU/1.5;

  //FORCE IN THE X DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelx);
  vecino1 = tex1Dfetch(texvecino1GPU, icelx);
  vecino2 = tex1Dfetch(texvecino2GPU, icelx);
  vecino3 = tex1Dfetch(texvecino3GPU, icelx);
  vecino4 = tex1Dfetch(texvecino4GPU, icelx);
  vecino5 = tex1Dfetch(texvecino5GPU, icelx);

  double dlS, dlpS, dlmS;

  r =  (rx - rxcellGPU[icelx] - dxGPU*0.5);
  rp = (rx - rxcellGPU[vecino3] - dxGPU*0.5);
  rm = (rx - rxcellGPU[vecino2] - dxGPU*0.5);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = tex1D(texDelta, fabs(r));
  dlxp = tex1D(texDelta, fabs(rp));
  dlxm = tex1D(texDelta, fabs(rm));
  dlS = functionDeltaDerived(1.5*r);
  dlpS = functionDeltaDerived(1.5*rp);
  dlmS = functionDeltaDerived(1.5*rm);


  r =  (ry - rycellGPU[icelx]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = tex1D(texDelta, fabs(r));
  dlyp = tex1D(texDelta, fabs(rp));
  dlym = tex1D(texDelta, fabs(rm));

  r =  (rz - rzcellGPU[icelx]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = tex1D(texDelta, fabs(r));
  dlzp = tex1D(texDelta, fabs(rp));
  dlzm = tex1D(texDelta, fabs(rm));

  double spxpypz, spxpy, spxpymz, spxpz, spx, spxmz, spxmypz, spxmy, spxmymz;
  double spypz, spy, spymz, spz, s, smz, smypz, smy, smymz;
  double smxpypz, smxpy, smxpymz, smxpz, smx, smxmz, smxmypz, smxmy, smxmymz;

  double massfluid = 1.5 * volumeGPU*volumeParticleGPU*densfluidGPU;
  
  spxpypz = dlpS * dlyp * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxpy   = dlpS * dlyp * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxpymz = dlpS * dlyp * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxpz   = dlpS * dly  * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;
  
  spx     = dlpS * dly  * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;//

  spxmz   = dlpS * dly  * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxmypz = dlpS * dlym * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxmy   = dlpS * dlym * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spxmymz = dlpS * dlym * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spypz   = dlS  * dlyp * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spy     = dlS  * dlyp * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spymz   = dlS  * dlyp * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  spz     = dlS  * dly  * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;
  
  s       = dlS  * dly  * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;//
  
  smz     = dlS  * dly  * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smypz   = dlS  * dlym * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smy     = dlS  * dlym * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smymz   = dlS  * dlym * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxpypz = dlmS * dlyp * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxpy   = dlmS * dlyp * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxpymz = dlmS * dlyp * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxpz   = dlmS * dly  * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;
  
  smx     = dlmS * dly  * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;//

  smxmz   = dlmS * dly  * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxmypz = dlmS * dlym * dlzp * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxmy   = dlmS * dlym * dlz  * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;

  smxmymz = dlmS * dlym * dlzm * temperatureGPU * massfluid *  invdxGPU / 
    (massParticleGPU+massfluid) ;



  int offset = nboundaryGPU;
  fxboundaryGPU[offset+i]        += spxpypz;
  offset += nboundaryGPU+npGPU;//1
  fxboundaryGPU[offset+i] +=  spxpy;
  offset += nboundaryGPU+npGPU;//2
  fxboundaryGPU[offset+i] +=  spxpymz;
  offset += nboundaryGPU+npGPU;//3
  fxboundaryGPU[offset+i] +=  spxpz;
  offset += nboundaryGPU+npGPU;//4
  fxboundaryGPU[offset+i] +=  spx;
  offset += nboundaryGPU+npGPU;//5
  fxboundaryGPU[offset+i] +=  spxmz;
  offset += nboundaryGPU+npGPU;//6
  fxboundaryGPU[offset+i] +=  spxmypz;
  offset += nboundaryGPU+npGPU;//7
  fxboundaryGPU[offset+i] +=  spxmy;
  offset += nboundaryGPU+npGPU;//8
  fxboundaryGPU[offset+i] +=  spxmymz;
  offset += nboundaryGPU+npGPU;//9
  fxboundaryGPU[offset+i] +=  spypz;
  offset += nboundaryGPU+npGPU;//10
  fxboundaryGPU[offset+i] +=  spy;
  offset += nboundaryGPU+npGPU;//11
  fxboundaryGPU[offset+i] +=  spymz;
  offset += nboundaryGPU+npGPU;//12
  fxboundaryGPU[offset+i] +=  spz;
  offset += nboundaryGPU+npGPU;//13
  fxboundaryGPU[offset+i] +=  s;
  offset += nboundaryGPU+npGPU;//14
  fxboundaryGPU[offset+i] +=  smz;
  offset += nboundaryGPU+npGPU;//15
  fxboundaryGPU[offset+i] +=  smypz;
  offset += nboundaryGPU+npGPU;//16
  fxboundaryGPU[offset+i] +=  smy;
  offset += nboundaryGPU+npGPU;//17
  fxboundaryGPU[offset+i] +=  smymz;
  offset += nboundaryGPU+npGPU;//18
  fxboundaryGPU[offset+i] +=  smxpypz;
  offset += nboundaryGPU+npGPU;//19
  fxboundaryGPU[offset+i] +=  smxpy;
  offset += nboundaryGPU+npGPU;//20
  fxboundaryGPU[offset+i] +=  smxpymz;
  offset += nboundaryGPU+npGPU;//21
  fxboundaryGPU[offset+i] +=  smxpz;
  offset += nboundaryGPU+npGPU;//22
  fxboundaryGPU[offset+i] +=  smx;
  offset += nboundaryGPU+npGPU;//23
  fxboundaryGPU[offset+i] +=  smxmz;
  offset += nboundaryGPU+npGPU;//24
  fxboundaryGPU[offset+i] +=  smxmypz;
  offset += nboundaryGPU+npGPU;//25
  fxboundaryGPU[offset+i] +=  smxmy;
  offset += nboundaryGPU+npGPU;//26
  fxboundaryGPU[offset+i] +=  smxmymz;
  







  //FORCE IN THE Y DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icely);
  vecino1 = tex1Dfetch(texvecino1GPU, icely);
  vecino2 = tex1Dfetch(texvecino2GPU, icely);
  vecino3 = tex1Dfetch(texvecino3GPU, icely);
  vecino4 = tex1Dfetch(texvecino4GPU, icely);
  vecino5 = tex1Dfetch(texvecino5GPU, icely);
  
  r =  (rx - rxcellGPU[icely]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = tex1D(texDelta, fabs(r));
  dlxp = tex1D(texDelta, fabs(rp));
  dlxm = tex1D(texDelta, fabs(rm));

  r =  (ry - rycellGPU[icely] - dyGPU*0.5);
  rp = (ry - rycellGPU[vecino4] - dyGPU*0.5);
  rm = (ry - rycellGPU[vecino1] - dyGPU*0.5); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = tex1D(texDelta, fabs(r));
  dlyp = tex1D(texDelta, fabs(rp));
  dlym = tex1D(texDelta, fabs(rm));
  dlS = functionDeltaDerived(1.5*r);
  dlpS = functionDeltaDerived(1.5*rp);
  dlmS = functionDeltaDerived(1.5*rm);

  r =  (rz - rzcellGPU[icely]);
  rp = (rz - rzcellGPU[vecino5]);
  rm = (rz - rzcellGPU[vecino0]);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = tex1D(texDelta, fabs(r));
  dlzp = tex1D(texDelta, fabs(rp));
  dlzm = tex1D(texDelta, fabs(rm));




  spxpypz = dlxp * dlpS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxpy   = dlxp * dlpS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxpymz = dlxp * dlpS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxpz   = dlxp * dlS  * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;
  
  spx     = dlxp * dlS  * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;//

  spxmz   = dlxp * dlS  * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxmypz = dlxp * dlmS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxmy   = dlxp * dlmS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spxmymz = dlxp * dlmS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spypz   = dlx * dlpS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spy     = dlx  * dlpS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spymz   = dlx  * dlpS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  spz     = dlx  * dlS  * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;
  
  s       = dlx  * dlS  * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;//
  
  smz     = dlx  * dlS  * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smypz   = dlx  * dlmS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smy     = dlx  * dlmS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smymz   = dlx  * dlmS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxpypz = dlxm * dlpS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxpy   = dlxm * dlpS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxpymz = dlxm * dlpS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxpz   = dlxm * dlS  * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;
  
  smx     = dlxm * dlS  * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;//

  smxmz   = dlxm * dlS  * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxmypz = dlxm * dlmS * dlzp * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxmy   = dlxm * dlmS * dlz  * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;

  smxmymz = dlxm * dlmS * dlzm * temperatureGPU * massfluid *  invdyGPU / 
    (massParticleGPU+massfluid) ;




  offset = nboundaryGPU;
  fyboundaryGPU[offset+i]        +=  spxpypz;
  offset += nboundaryGPU+npGPU;//1
  fyboundaryGPU[offset+i] +=  spxpy;
  offset += nboundaryGPU+npGPU;//2
  fyboundaryGPU[offset+i] +=  spxpymz;
  offset += nboundaryGPU+npGPU;//3
  fyboundaryGPU[offset+i] +=  spxpz;
  offset += nboundaryGPU+npGPU;//4
  fyboundaryGPU[offset+i] +=  spx;
  offset += nboundaryGPU+npGPU;//5
  fyboundaryGPU[offset+i] +=  spxmz;
  offset += nboundaryGPU+npGPU;//6
  fyboundaryGPU[offset+i] +=  spxmypz;
  offset += nboundaryGPU+npGPU;//7
  fyboundaryGPU[offset+i] +=  spxmy;
  offset += nboundaryGPU+npGPU;//8
  fyboundaryGPU[offset+i] += + spxmymz;
  offset += nboundaryGPU+npGPU;//9
  fyboundaryGPU[offset+i] +=  spypz;
  offset += nboundaryGPU+npGPU;//10
  fyboundaryGPU[offset+i] +=  spy;
  offset += nboundaryGPU+npGPU;//11
  fyboundaryGPU[offset+i] +=  spymz;
  offset += nboundaryGPU+npGPU;//12
  fyboundaryGPU[offset+i] +=  spz;
  offset += nboundaryGPU+npGPU;//13
  fyboundaryGPU[offset+i] +=  s;
  offset += nboundaryGPU+npGPU;//14
  fyboundaryGPU[offset+i] +=  smz;
  offset += nboundaryGPU+npGPU;//15
  fyboundaryGPU[offset+i] +=  smypz;
  offset += nboundaryGPU+npGPU;//16
  fyboundaryGPU[offset+i] +=  smy;
  offset += nboundaryGPU+npGPU;//17
  fyboundaryGPU[offset+i] +=  smymz;
  offset += nboundaryGPU+npGPU;//18
  fyboundaryGPU[offset+i] +=  smxpypz;
  offset += nboundaryGPU+npGPU;//19
  fyboundaryGPU[offset+i] +=  smxpy;
  offset += nboundaryGPU+npGPU;//20
  fyboundaryGPU[offset+i] +=  smxpymz;
  offset += nboundaryGPU+npGPU;//21
  fyboundaryGPU[offset+i] +=  smxpz;
  offset += nboundaryGPU+npGPU;//22
  fyboundaryGPU[offset+i] +=  smx;
  offset += nboundaryGPU+npGPU;//23
  fyboundaryGPU[offset+i] +=  smxmz;
  offset += nboundaryGPU+npGPU;//24
  fyboundaryGPU[offset+i] +=  smxmypz;
  offset += nboundaryGPU+npGPU;//25
  fyboundaryGPU[offset+i] +=  smxmy;
  offset += nboundaryGPU+npGPU;//26
  fyboundaryGPU[offset+i] +=  smxmymz;


  
  
  //FORCE IN THE Z DIRECTION
  vecino0 = tex1Dfetch(texvecino0GPU, icelz);
  vecino1 = tex1Dfetch(texvecino1GPU, icelz);
  vecino2 = tex1Dfetch(texvecino2GPU, icelz);
  vecino3 = tex1Dfetch(texvecino3GPU, icelz);
  vecino4 = tex1Dfetch(texvecino4GPU, icelz);
  vecino5 = tex1Dfetch(texvecino5GPU, icelz);
 
  r =  (rx - rxcellGPU[icelz]);
  rp = (rx - rxcellGPU[vecino3]);
  rm = (rx - rxcellGPU[vecino2]);
  r =  auxdx * (r - int(r*invlxGPU + 0.5*((r>0)-(r<0)))*lxGPU);
  rm = auxdx * (rm - int(rm*invlxGPU + 0.5*((rm>0)-(rm<0)))*lxGPU);
  rp = auxdx * (rp - int(rp*invlxGPU + 0.5*((rp>0)-(rp<0)))*lxGPU);
  dlx = tex1D(texDelta, fabs(r));
  dlxp = tex1D(texDelta, fabs(rp));
  dlxm = tex1D(texDelta, fabs(rm));

  r =  (ry - rycellGPU[icelz]);
  rp = (ry - rycellGPU[vecino4]);
  rm = (ry - rycellGPU[vecino1]); 
  r =  auxdy * (r - int(r*invlyGPU + 0.5*((r>0)-(r<0)))*lyGPU);
  rp = auxdy * (rp - int(rp*invlyGPU + 0.5*((rp>0)-(rp<0)))*lyGPU);
  rm = auxdy * (rm - int(rm*invlyGPU + 0.5*((rm>0)-(rm<0)))*lyGPU);
  dly = tex1D(texDelta, fabs(r));
  dlyp = tex1D(texDelta, fabs(rp));
  dlym = tex1D(texDelta, fabs(rm));

  r =  (rz - rzcellGPU[icelz] - dzGPU*0.5);
  rp = (rz - rzcellGPU[vecino5] - dzGPU*0.5);
  rm = (rz - rzcellGPU[vecino0] - dzGPU*0.5);
  r = auxdz * (r - int(r*invlzGPU + 0.5*((r>0)-(r<0)))*lzGPU);
  rp = auxdz * (rp - int(rp*invlzGPU + 0.5*((rp>0)-(rp<0)))*lzGPU);
  rm = auxdz * (rm - int(rm*invlzGPU + 0.5*((rm>0)-(rm<0)))*lzGPU);
  dlz = tex1D(texDelta, fabs(r));
  dlzp = tex1D(texDelta, fabs(rp));
  dlzm = tex1D(texDelta, fabs(rm));
  dlS = functionDeltaDerived(1.5*r);
  dlpS = functionDeltaDerived(1.5*rp);
  dlmS = functionDeltaDerived(1.5*rm);



  spxpypz = dlxp * dlyp * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;
  
  spxpy   = dlxp * dlyp * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spxpymz = dlxp * dlyp * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spxpz   = dlxp * dly  * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;
  
  spx     = dlxp * dly  * dlS  * temperatureGPU * massfluid *  invdzGPU / 
     (massParticleGPU+massfluid) ;//

  spxmz   = dlxp * dly  * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spxmypz = dlxp * dlym * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spxmy   = dlxp * dlym * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spxmymz = dlxp * dlym * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spypz   = dlx  * dlyp * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spy     = dlx  * dlyp * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spymz   = dlx  * dlyp * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  spz     = dlx  * dly  * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;
  
  s       = dlx  * dly  * dlS  * temperatureGPU * massfluid *  invdzGPU / 
     (massParticleGPU+massfluid) ;//
  
  smz     = dlx  * dly  * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smypz   = dlx  * dlym * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smy     = dlx  * dlym * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smymz   = dlx  * dlym * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxpypz = dlxm * dlyp * dlpS * temperatureGPU * massfluid *  invdzGPU / 
     (massParticleGPU+massfluid) ;

  smxpy   = dlxm * dlyp * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxpymz = dlxm * dlyp * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxpz   = dlxm * dly  * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;
  
  smx     = dlxm * dly  * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;//

  smxmz   = dlxm * dly  * dlmS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxmypz = dlxm * dlym * dlpS * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxmy   = dlxm * dlym * dlS  * temperatureGPU * massfluid *  invdzGPU / 
    (massParticleGPU+massfluid) ;

  smxmymz = dlxm * dlym * dlmS * temperatureGPU * massfluid *  invdzGPU / 
     (massParticleGPU+massfluid) ;
  


  offset = nboundaryGPU;
  fzboundaryGPU[offset+i] += spxpypz;
  offset += nboundaryGPU+npGPU;//1
  fzboundaryGPU[offset+i] += spxpy;
  offset += nboundaryGPU+npGPU;//2
  fzboundaryGPU[offset+i] += spxpymz;
  offset += nboundaryGPU+npGPU;//3
  fzboundaryGPU[offset+i] += spxpz;
  offset += nboundaryGPU+npGPU;//4
  fzboundaryGPU[offset+i] += spx;
  offset += nboundaryGPU+npGPU;//5
  fzboundaryGPU[offset+i] += spxmz;
  offset += nboundaryGPU+npGPU;//6
  fzboundaryGPU[offset+i] += spxmypz;
  offset += nboundaryGPU+npGPU;//7
  fzboundaryGPU[offset+i] += spxmy;
  offset += nboundaryGPU+npGPU;//8
  fzboundaryGPU[offset+i] += spxmymz;
  offset += nboundaryGPU+npGPU;//9
  fzboundaryGPU[offset+i] += spypz;
  offset += nboundaryGPU+npGPU;//10
  fzboundaryGPU[offset+i] += spy;
  offset += nboundaryGPU+npGPU;//11
  fzboundaryGPU[offset+i] += spymz;
  offset += nboundaryGPU+npGPU;//12
  fzboundaryGPU[offset+i] += spz;
  offset += nboundaryGPU+npGPU;//13
  fzboundaryGPU[offset+i] += s;
  offset += nboundaryGPU+npGPU;//14
  fzboundaryGPU[offset+i] += smz;
  offset += nboundaryGPU+npGPU;//15
  fzboundaryGPU[offset+i] += smypz;
  offset += nboundaryGPU+npGPU;//16
  fzboundaryGPU[offset+i] += smy;
  offset += nboundaryGPU+npGPU;//17
  fzboundaryGPU[offset+i] += smymz;
  offset += nboundaryGPU+npGPU;//18
  fzboundaryGPU[offset+i] += smxpypz;
  offset += nboundaryGPU+npGPU;//19
  fzboundaryGPU[offset+i] += smxpy;
  offset += nboundaryGPU+npGPU;//20
  fzboundaryGPU[offset+i] += smxpymz;
  offset += nboundaryGPU+npGPU;//21
  fzboundaryGPU[offset+i] += smxpz;
  offset += nboundaryGPU+npGPU;//22
  fzboundaryGPU[offset+i] += smx;
  offset += nboundaryGPU+npGPU;//23
  fzboundaryGPU[offset+i] += smxmz;
  offset += nboundaryGPU+npGPU;//24
  fzboundaryGPU[offset+i] += smxmypz;
  offset += nboundaryGPU+npGPU;//25
  fzboundaryGPU[offset+i] += smxmy;
  offset += nboundaryGPU+npGPU;//26
  fzboundaryGPU[offset+i] += smxmymz;




  
}
