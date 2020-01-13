// Filename: runSchemeParticlesWall.cu
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


bool runSchemeParticlesWall(){
  int threadsPerBlock = 512;
  if((ncells/threadsPerBlock) < 512) threadsPerBlock = 256;
  if((ncells/threadsPerBlock) < 256) threadsPerBlock = 128;
  if((ncells/threadsPerBlock) < 64) threadsPerBlock = 64;
  if((ncells/threadsPerBlock) < 64) threadsPerBlock = 32;
  int numBlocks = (ncells-1)/threadsPerBlock + 1;

  int threadsPerBlockBoundary = 512;
  if((nboundary/threadsPerBlockBoundary) < 512) threadsPerBlockBoundary = 256;
  if((nboundary/threadsPerBlockBoundary) < 256) threadsPerBlockBoundary = 128;
  if((nboundary/threadsPerBlockBoundary) < 64) threadsPerBlockBoundary = 64;
  if((nboundary/threadsPerBlockBoundary) < 64) threadsPerBlockBoundary = 32;
  int numBlocksBoundary = (nboundary-1)/threadsPerBlockBoundary + 1;

  int threadsPerBlockPartAndBoundary = 512;
  if(((np+nboundary)/threadsPerBlockPartAndBoundary) < 512) threadsPerBlockPartAndBoundary = 256;
  if(((np+nboundary)/threadsPerBlockPartAndBoundary) < 256) threadsPerBlockPartAndBoundary = 128;
  if(((np+nboundary)/threadsPerBlockPartAndBoundary) < 64) threadsPerBlockPartAndBoundary = 64;
  if(((np+nboundary)/threadsPerBlockPartAndBoundary) < 64) threadsPerBlockPartAndBoundary = 32;
  int numBlocksPartAndBoundary = (np+nboundary-1)/threadsPerBlockPartAndBoundary + 1;

  int threadsPerBlockParticles = 512;
  if((np/threadsPerBlockParticles) < 512) threadsPerBlockParticles = 256;
  if((np/threadsPerBlockParticles) < 256) threadsPerBlockParticles = 128;
  if((np/threadsPerBlockParticles) < 64) threadsPerBlockParticles = 64;
  if((np/threadsPerBlockParticles) < 64) threadsPerBlockParticles = 32;
  int numBlocksParticles = (np-1)/threadsPerBlockParticles + 1;

  int threadsPerBlockNeighbors, numBlocksNeighbors;
  if(ncellst>numNeighbors){
    threadsPerBlockNeighbors = 512;
    if((ncellst/threadsPerBlockNeighbors) < 512) threadsPerBlockNeighbors = 256;
    if((ncellst/threadsPerBlockNeighbors) < 256) threadsPerBlockNeighbors = 128;
    if((ncellst/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 64;
    if((ncellst/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 32;
    numBlocksNeighbors = (ncellst-1)/threadsPerBlockNeighbors + 1;
  }
  else{
    threadsPerBlockNeighbors = 512;
    if((numNeighbors/threadsPerBlockNeighbors) < 512) threadsPerBlockNeighbors = 256;
    if((numNeighbors/threadsPerBlockNeighbors) < 256) threadsPerBlockNeighbors = 128;
    if((numNeighbors/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 64;
    if((numNeighbors/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 32;
    numBlocksNeighbors = (numNeighbors-1)/threadsPerBlockNeighbors + 1;
  }
  int nGhost = ncellst - ncells;
  int threadsPerBlockGhost = 512;
  if((nGhost/threadsPerBlockGhost) < 512) threadsPerBlockGhost = 256;
  if((nGhost/threadsPerBlockGhost) < 256) threadsPerBlockGhost = 128;
  if((nGhost/threadsPerBlockGhost) < 64) threadsPerBlockGhost = 64;
  if((nGhost/threadsPerBlockGhost) < 64) threadsPerBlockGhost = 32;
  int numBlocksGhost = (nGhost-1)/threadsPerBlockGhost + 1;

  int threadsPerBlockCellst = 512;
  if((ncellst/threadsPerBlockCellst) < 512) threadsPerBlockCellst = 256;
  if((ncellst/threadsPerBlockCellst) < 256) threadsPerBlockCellst = 128;
  if((ncellst/threadsPerBlockCellst) < 64) threadsPerBlockCellst = 64;
  if((ncellst/threadsPerBlockCellst) < 64) threadsPerBlockCellst = 32;
  int numBlocksCellst = (ncellst-1)/threadsPerBlockCellst + 1;


  //initialize random numbers
  size_t numberRandom = 12 * ncellst;
  if(!initializeRandomNumbersGPU(numberRandom,seed)) return 0;


  //Initialize textures cells
  if(!texturesCells()) return 0;
  
  //Inilialize ghost index
  if(!initGhostIndexParticlesWallGPU()) return 0;


  initializeVecinos<<<numBlocksCellst,threadsPerBlockCellst>>>(vecino1GPU,
							       vecino2GPU,
							       vecino3GPU,
							       vecino4GPU,
							       vecinopxpyGPU,
							       vecinopxmyGPU,
							       vecinopxpzGPU,
							       vecinopxmzGPU,
							       vecinomxpyGPU,
							       vecinomxmyGPU,
							       vecinomxpzGPU,
							       vecinomxmzGPU,
							       vecinopypzGPU,
							       vecinopymzGPU,
							       vecinomypzGPU,
							       vecinomymzGPU,
							       vecinopxpypzGPU,
							       vecinopxpymzGPU,
							       vecinopxmypzGPU,
							       vecinopxmymzGPU,
							       vecinomxpypzGPU,
							       vecinomxpymzGPU,
							       vecinomxmypzGPU,
							       vecinomxmymzGPU);
  
  initializeVecinos2<<<numBlocksCellst,threadsPerBlockCellst>>>(vecino0GPU,
								vecino1GPU,
								vecino2GPU,
								vecino3GPU,
								vecino4GPU,
								vecino5GPU);

  /********* Added by Adolfo **********/
  // The first step is outputted
  if(!gpuToHostParticles()) return 0;
  if(!saveFunctionsSchemeParticlesWall(1,0)) return 0;

  /******** Added by Adolfo **********/
  double freq    = freqWall;
  // VxWall0 has taken the value given at the input file
  double vxWall  = vxWall0;
  double xWall = 0.0; // initial position of the wall
  double *Fatt;
  Fatt = new double [nboundary + np];
  *Fatt = 0;
  double *counter;
  counter = new double [np];
  *counter = 0;
  double polymer_force_wall = 0;
  double mass_wall = massWall;
  // double k_wall = mass_wall * freq * freq;
  cout << "mass_wall = " << mass_wall << endl;
  // cout << "k_wall = " << k_wall << endl;
  /**********************************/

  int substep = 0;
  step = -numstepsRelaxation;

  while(step<numsteps){
    /***** Added by Adolfo ****/
    double time = dt * step;
    /**************************/

    /**** Added by Adolfo ****/
    /** The shear rate is calculated later at kernelDpParticlesWall_shearRate **/
    cutilSafeCall(cudaMemcpy(cdydvx,dydvxGPU,ncellst*sizeof(double),cudaMemcpyDeviceToHost));

    // Calculation of the total viscous force 
    double visc_force = 0;
    if (step != 0)
      for(int j=0;j<ncells;j++){
	int i, fx, fy, fz;
	fx = j % mx;
	fy = (j % (mx*my)) / mx;
	fz = j / (mx*my);
	//go to index tanking into account ghost cells
	fy++;
	i = fx + fy*mxt + fz*mxt*myt;
	if (fy == 1)
	  {
	    visc_force = visc_force + cdydvx[i];
	  }
      }
    else // step = 0
       visc_force = 0;

    // visc_force before here is the sumation of the shear rate along all the cells.
    // That is why it is multiplied by the area of the cell, because 
    visc_force = visc_force * shearviscosity * lx/mx * lz/mz; 
    
    cutilSafeCall(cudaMemcpy(Fatt,FattGPU,np*sizeof(double),cudaMemcpyDeviceToHost));

    
    /*** The total force on the wall is calculated ***/
    polymer_force_wall = 0;
    if (step != 0)
      for(int aa=0;aa<np;aa++)
	polymer_force_wall += Fatt[aa]; 

    /*
    if (step > 1) 
      for(int aa=0;aa<np;aa++)
	cout << "AAA " << step << " " << aa << " " << Fatt[aa] << endl;
    */
	
	
    //double Force_total = visc_force + polymer_force_wall - k_wall * xWall;
    double Force_total = visc_force + polymer_force_wall;

    // cout << visc_force << " " << polymer_force_wall << endl;

    /* The next was for the ring-down */
    /* vxWall = vxWall + Force_total / mass_wall * dt; 
       xWall = xWall + vxWall * dt; */

    /** Enforcing the wall velocity **/
    vxWall = vxWall0 * cos(freq * time); 
    xWall = vxWall0 / freq * sin(freq * time);


    if(!(step%samplefreq)){
	cout << "wall " << time << " " << xWall << " " << vxWall << endl;
	cout << "force " << time << " " << visc_force << " " << polymer_force_wall << " " << Force_total << endl;
      }

    /********* End of Added by Adolfo ****************/

    //Generate random numbers
    generateRandomNumbers(numberRandom);

    //Provide data to ghost cells
    /*kernelFeedGhostCellsParticlesWall<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand);*/

    /********* Commented by Adolfo ******/
    /*
    kernelFeedGhostCellsParticlesWall2<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);
    */
    /*********** Added by Adolfo *********/
    kernelFeedGhostCellsParticlesWall2_time<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand,
       vxWall,
       time);

    //Boundaries and particles 
    //Update particles position to q^{n+1/2}
    //Spread and interpolate particles force
    /******* Changed by Adolfo ********/
    /*
    boundaryParticlesFunctionParticlesWall(0,
					   numBlocksBoundary,
					   threadsPerBlockBoundary,
					   numBlocksNeighbors,
					   threadsPerBlockNeighbors,
					   numBlocksPartAndBoundary,
					   threadsPerBlockPartAndBoundary,
					   numBlocksParticles,
					   threadsPerBlockParticles,
					   numBlocks,
					   threadsPerBlock,
					   numBlocksCellst,
					   threadsPerBlockCellst);
    */
    boundaryParticlesFunctionParticlesWall(0,
					   numBlocksBoundary,
					   threadsPerBlockBoundary,
					   numBlocksNeighbors,
					   threadsPerBlockNeighbors,
					   numBlocksPartAndBoundary,
					   threadsPerBlockPartAndBoundary,
					   numBlocksParticles,
					   threadsPerBlockParticles,
					   numBlocks,
					   threadsPerBlock,
					   numBlocksCellst,
					   threadsPerBlockCellst,
					   xWall);

    //Provide data to ghost cells
    /*kernelFeedGhostCellsParticlesWall<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand);*/

    /******** Commented by Adolfo *******/
    /*
    kernelFeedGhostCellsParticlesWall2<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);
    */

    /*********** Added by Adolfo *********/
    kernelFeedGhostCellsParticlesWall2_time<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand,
       vxWall,
       time);


    //First substep RK3
    kernelDpParticlesWall<<<numBlocksCellst,threadsPerBlockCellst>>>(densityGPU,
								     densityGPU,
								     vxGPU,
								     vyGPU,
								     vzGPU,
								     dmGPU,
								     dpxGPU,
								     dpyGPU,
								     dpzGPU,
								     dRand,
								     fxboundaryGPU,
								     fyboundaryGPU,
								     fzboundaryGPU,
								     ghostIndexGPU, 
								     realIndexGPU,
								     substep,
								     0,1,-sqrt(3));

    kernelDpParticlesWall_2<<<numBlocksCellst,threadsPerBlockCellst>>>(densityPredictionGPU,
								       vxPredictionGPU,
								       vyPredictionGPU,
								       vzPredictionGPU,
								       dmGPU,
								       dpxGPU,
								       dpyGPU,
								       dpzGPU,
								       ghostIndexGPU, 
								       realIndexGPU);
    
    cutilSafeCall( cudaBindTexture(0,texVxGPU,vxPredictionGPU,ncellst*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVyGPU,vyPredictionGPU,ncellst*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVzGPU,vzPredictionGPU,ncellst*sizeof(double)));

    //Provide data to ghost cells
    /*kernelFeedGhostCellsParticlesWall<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand);*/

    /************ Commented by Adolfo **********/
    /*
    kernelFeedGhostCellsParticlesWall2<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);
    */

    /*********** Added by Adolfo *********/
    kernelFeedGhostCellsParticlesWall2_time<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand,
       vxWall,
       time);
 
    //Second substep RK3
    kernelDpParticlesWall<<<numBlocksCellst,threadsPerBlockCellst>>>(densityPredictionGPU,
								     densityGPU,
								     vxGPU,
								     vyGPU,
								     vzGPU,
								     dmGPU,
								     dpxGPU,
								     dpyGPU,
								     dpzGPU,
								     dRand,
								     fxboundaryGPU,
								     fyboundaryGPU,
								     fzboundaryGPU,
								     ghostIndexGPU, 
								     realIndexGPU,
								     substep,
								     0.75,0.25,sqrt(3));
    
    kernelDpParticlesWall_2<<<numBlocksCellst,threadsPerBlockCellst>>>(densityPredictionGPU,
								       vxPredictionGPU,
								       vyPredictionGPU,
								       vzPredictionGPU,
								       dmGPU,
								       dpxGPU,
								       dpyGPU,
								       dpzGPU,
								       ghostIndexGPU, 
								       realIndexGPU);
    
    //Provide data to ghost cells
    /*kernelFeedGhostCellsParticlesWall<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand);*/

    /*********** Commented by Adolfo *********/
    /*
    kernelFeedGhostCellsParticlesWall2<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);
    */

    /*********** Added by Adolfo *********/
    kernelFeedGhostCellsParticlesWall2_time<<<numBlocksGhost,threadsPerBlockGhost>>>
      (ghostToPIGPU,
       ghostToGhostGPU,
       densityGPU,
       densityPredictionGPU,
       vxGPU,
       vyGPU,
       vzGPU,
       vxPredictionGPU,
       vyPredictionGPU,
       vzPredictionGPU,
       dRand,
       vxWall,
       time);

    //Third substep RK3
    kernelDpParticlesWall<<<numBlocksCellst,threadsPerBlockCellst>>>(densityPredictionGPU,
								     densityGPU,
								     vxGPU,
								     vyGPU,
								     vzGPU,
								     dmGPU,
								     dpxGPU,
								     dpyGPU,
								     dpzGPU,
								     dRand,
								     fxboundaryGPU,
								     fyboundaryGPU,
								     fzboundaryGPU,
								     ghostIndexGPU, 
								     realIndexGPU,
								     substep,
								     1./3.,2./3.,0);

    //Copy v^n to vxPredictionGPU
    //We need it for the particle update
    copyField<<<numBlocksCellst,threadsPerBlockCellst>>>(vxGPU,
							 vyGPU,
							 vzGPU,
							 vxPredictionGPU,
							 vyPredictionGPU,
							 vzPredictionGPU);
    
    kernelDpParticlesWall_2<<<numBlocksCellst,threadsPerBlockCellst>>>(densityGPU,
								       vxGPU,
								       vyGPU,
								       vzGPU,
								       dmGPU,
								       dpxGPU,
								       dpyGPU,
								       dpzGPU,
								       ghostIndexGPU, 
								       realIndexGPU);

    /**** Added by Adolfo ***/
    kernelDpParticlesWall_shearRate<<<numBlocksCellst,threadsPerBlockCellst>>>(densityGPU,
									       vxGPU,
									       vyGPU,
									       vzGPU,
									       dmGPU,
									       dpxGPU,
									       dpyGPU,
									       dpzGPU,
									       dydvxGPU,
									       vxWall,
									       ghostIndexGPU, 
									       realIndexGPU);

    cutilSafeCall( cudaBindTexture(0,texVxGPU,vxGPU,ncellst*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVyGPU,vyGPU,ncellst*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVzGPU,vzGPU,ncellst*sizeof(double)));


    //Boundaries and particles part start
    /******* Changed by Adolfo ********/
    /*
    boundaryParticlesFunctionParticlesWall(1,
					   numBlocksBoundary,
					   threadsPerBlockBoundary,
					   numBlocksNeighbors,
					   threadsPerBlockNeighbors,
					   numBlocksPartAndBoundary,
					   threadsPerBlockPartAndBoundary,
					   numBlocksParticles,
					   threadsPerBlockParticles,
					   numBlocks,
					   threadsPerBlock,
					   numBlocksCellst,
					   threadsPerBlockCellst);
    */
    boundaryParticlesFunctionParticlesWall(1,
					   numBlocksBoundary,
					   threadsPerBlockBoundary,
					   numBlocksNeighbors,
					   threadsPerBlockNeighbors,
					   numBlocksPartAndBoundary,
					   threadsPerBlockPartAndBoundary,
					   numBlocksParticles,
					   threadsPerBlockParticles,
					   numBlocks,
					   threadsPerBlock,
					   numBlocksCellst,
					   threadsPerBlockCellst,
					   xWall);


    step++;
        
    if(!(step%samplefreq)&&(step>0)){
      cout << "Particles Wall  " << step << endl;

      //Provide data to ghost cells
      /*kernelFeedGhostCellsParticlesWall<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);*/

      /********* Commented by Adolfo ***********/
      /*      
      kernelFeedGhostCellsParticlesWall2<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand);
      */

      /*********** Added by Adolfo *********/
      kernelFeedGhostCellsParticlesWall2_time<<<numBlocksGhost,threadsPerBlockGhost>>>
	(ghostToPIGPU,
	 ghostToGhostGPU,
	 densityGPU,
	 densityPredictionGPU,
	 vxGPU,
	 vyGPU,
	 vzGPU,
	 vxPredictionGPU,
	 vyPredictionGPU,
	 vzPredictionGPU,
	 dRand,
	 vxWall,
	 time);


      if(!gpuToHostParticles()) return 0;
      if(!saveFunctionsSchemeParticlesWall(1,step)) return 0;

    }
    
  }


  freeRandomNumbersGPU();

  return 1;
}
