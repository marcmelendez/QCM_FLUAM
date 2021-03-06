// Filename: runSchemeCompressibleParticles.cu
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


bool runSchemeMomentumCoupling(){
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
  if(ncells>numNeighbors){
    threadsPerBlockNeighbors = 512;
    if((ncells/threadsPerBlockNeighbors) < 512) threadsPerBlockNeighbors = 256;
    if((ncells/threadsPerBlockNeighbors) < 256) threadsPerBlockNeighbors = 128;
    if((ncells/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 64;
    if((ncells/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 32;
    numBlocksNeighbors = (ncells-1)/threadsPerBlockNeighbors + 1;
  }
  else{
    threadsPerBlockNeighbors = 512;
    if((numNeighbors/threadsPerBlockNeighbors) < 512) threadsPerBlockNeighbors = 256;
    if((numNeighbors/threadsPerBlockNeighbors) < 256) threadsPerBlockNeighbors = 128;
    if((numNeighbors/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 64;
    if((numNeighbors/threadsPerBlockNeighbors) < 64) threadsPerBlockNeighbors = 32;
    numBlocksNeighbors = (numNeighbors-1)/threadsPerBlockNeighbors + 1;
  }



  step = -numstepsRelaxation;

  //initialize random numbers
  size_t numberRandom = 12 * ncells;
  if(!initializeRandomNumbersGPU(numberRandom,seed)) return 0;

  //Initialize textures cells
  if(!texturesCells()) return 0;

  initializeVecinos<<<numBlocks,threadsPerBlock>>>(vecino1GPU,
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

  initializeVecinos2<<<numBlocks,threadsPerBlock>>>(vecino0GPU,
						    vecino1GPU,
						    vecino2GPU,
						    vecino3GPU,
						    vecino4GPU,
						    vecino5GPU);


  while(step<numsteps){
    //Generate random numbers
    generateRandomNumbers(numberRandom);

    //Boundaries and particles 
    //Update particles position to q^{n+1/2}
    //Spread and interpolate particles force
    boundaryParticlesFunctionMomentumCoupling(0,
					      //boundaryParticlesFunctionFreeEnergyCompressibleParticles(0,
					      numBlocksBoundary,
					      threadsPerBlockBoundary,
					      numBlocksNeighbors,
					      threadsPerBlockNeighbors,
					      numBlocksPartAndBoundary,
					      threadsPerBlockPartAndBoundary,
					      numBlocksParticles,
					      threadsPerBlockParticles,
					      numBlocks,
					      threadsPerBlock);
   

    //First substep RK3
    kernelDpFreeEnergyCompressibleParticles<<<numBlocks,threadsPerBlock>>>(densityGPU,
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
									   omegaGPU,
									   step,
									   0,1,-sqrt(3));

    kernelDpThermostat_2<<<numBlocks,threadsPerBlock>>>(densityPredictionGPU,
							vxPredictionGPU,
							vyPredictionGPU,
							vzPredictionGPU,
							dmGPU,
							dpxGPU,
							dpyGPU,
							dpzGPU);

    cutilSafeCall( cudaBindTexture(0,texVxGPU,vxPredictionGPU,ncells*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVyGPU,vyPredictionGPU,ncells*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVzGPU,vzPredictionGPU,ncells*sizeof(double)));

    //Second substep RK3

    //Boundaries and particles 
    //Interpolate fluid density
    boundaryParticlesFunctionMomentumCoupling(1,
    //boundaryParticlesFunctionFreeEnergyCompressibleParticles(1,
					      numBlocksBoundary,
					      threadsPerBlockBoundary,
					      numBlocksNeighbors,
					      threadsPerBlockNeighbors,
					      numBlocksPartAndBoundary,
					      threadsPerBlockPartAndBoundary,
					      numBlocksParticles,
					      threadsPerBlockParticles,
					      numBlocks,
					      threadsPerBlock);

    kernelDpFreeEnergyCompressibleParticles<<<numBlocks,threadsPerBlock>>>(densityPredictionGPU,
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
									   omegaGPU,
									   step,
									   0.75,0.25,sqrt(3));
    
    kernelDpThermostat_2<<<numBlocks,threadsPerBlock>>>(densityPredictionGPU,
							vxPredictionGPU,
							vyPredictionGPU,
							vzPredictionGPU,
							dmGPU,
							dpxGPU,
							dpyGPU,
							dpzGPU);

    //Third substep RK3

    //Boundaries and particles 
    //Interpolate fluid density
    boundaryParticlesFunctionMomentumCoupling(1,
    //boundaryParticlesFunctionFreeEnergyCompressibleParticles(2,
					      numBlocksBoundary,
					      threadsPerBlockBoundary,
					      numBlocksNeighbors,
					      threadsPerBlockNeighbors,
					      numBlocksPartAndBoundary,
					      threadsPerBlockPartAndBoundary,
					      numBlocksParticles,
					      threadsPerBlockParticles,
					      numBlocks,
					      threadsPerBlock);

    kernelDpFreeEnergyCompressibleParticles<<<numBlocks,threadsPerBlock>>>(densityPredictionGPU,
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
									   omegaGPU,
									   step,
									   1./3.,2./3.,0);
    //Copy v^n to vxPredictionGPU
    //We need it for the particle update
    copyField<<<numBlocks,threadsPerBlock>>>(vxGPU,
					     vyGPU,
					     vzGPU,
					     vxPredictionGPU,
					     vyPredictionGPU,
					     vzPredictionGPU);

    kernelDpThermostat_2<<<numBlocks,threadsPerBlock>>>(densityGPU,
							vxGPU,
							vyGPU,
							vzGPU,
							dmGPU,
							dpxGPU,
							dpyGPU,
							dpzGPU);

    cutilSafeCall( cudaBindTexture(0,texVxGPU,vxGPU,ncells*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVyGPU,vyGPU,ncells*sizeof(double)));
    cutilSafeCall( cudaBindTexture(0,texVzGPU,vzGPU,ncells*sizeof(double)));


    //Boundaries and particles part start
    boundaryParticlesFunctionMomentumCoupling(2,
    //boundaryParticlesFunctionFreeEnergyCompressibleParticles(2,
					      numBlocksBoundary,
					      threadsPerBlockBoundary,
					      numBlocksNeighbors,
					      threadsPerBlockNeighbors,
					      numBlocksPartAndBoundary,
					      threadsPerBlockPartAndBoundary,
					      numBlocksParticles,
					      threadsPerBlockParticles,
					      numBlocks,
					      threadsPerBlock);
    
    
    step++;
    //cout << step << endl;
    
    if(!(step%samplefreq)&&(step>0)){
      cout << "Compressible  " << step << endl;

      //We want to save interplate density
      //and save it in rxboundaryPrediction
      if(1){//We save density
	//if(0){//We don't save density
	interpolateDensity<<<numBlocksParticles,threadsPerBlockParticles>>>(rxcellGPU,
									    rycellGPU,
									    rzcellGPU,
									    densityGPU,
									    rxboundaryGPU,
									    ryboundaryGPU,
									    rzboundaryGPU,
									    rxboundaryPredictionGPU);//J*rho
	cutilSafeCall(cudaMemcpy(vxParticleI,
				 &rxboundaryPredictionGPU[nboundary],np*sizeof(double),cudaMemcpyDeviceToHost));
      }


      if(!gpuToHostParticles()) return 0;
      if(!saveFunctionsSchemeBoundary(1,step)) return 0;
      
      /*if(!(step%40000)){
	cout  << endl << endl << "OLD TEMPERATURE " << temperature << endl;
	temperature -= 1e-4;
	double fact1 = sqrt((4.*temperature*shearviscosity)/(dt*cVolume));
	double fact2 = sqrt((2.*temperature*bulkviscosity)/(3.*dt*cVolume));
	double fact4 = sqrt((2.*temperature*shearviscosity)/(dt*cVolume));
	cutilSafeCall(cudaMemcpyToSymbol(temperatureGPU,&temperature,sizeof(double)));
	cutilSafeCall(cudaMemcpyToSymbol(fact1GPU,&fact1,sizeof(double)));
	cutilSafeCall(cudaMemcpyToSymbol(fact2GPU,&fact2,sizeof(double)));
	cutilSafeCall(cudaMemcpyToSymbol(fact4GPU,&fact4,sizeof(double)));
	cout  << "NEW TEMPERATURE " << temperature << endl << endl;
	}*/
      cudaThreadSynchronize();     
    }
  
  }



  freeRandomNumbersGPU();




  return 1;
}
