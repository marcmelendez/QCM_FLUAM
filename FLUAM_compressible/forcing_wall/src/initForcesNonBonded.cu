// Filename: initForcesNonBonded.cu
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


float functionForceNonBonded1(double r){
  float sigma, epsilon;
  sigma = 2 * lx / double(mx);
  epsilon = temperature ;
  //return -epsilon * ( r - sigma);
  /********* Changed by Adolfo *********/
  // return 48. * epsilon * (pow(sigma/r,12) - 0.5*pow(sigma/r,6))/r;
  return 0;


}
float functionForceNonBonded1(double r, double a, double b){
  //  float sigma, epsilon;
  //sigma = 2 * lx / double(mx);
  //epsilon = temperature ;
  //return -epsilon * ( r - sigma);
  return (a*pow(1.0/r,12) - b*pow(1.0/r,6))/r;
}

bool initForcesNonBonded(){
  texforceNonBonded1.normalized = true;
  texforceNonBonded1.addressMode[0] = cudaAddressModeClamp;//Wrap and Clamp
  texforceNonBonded1.filterMode = cudaFilterModeLinear;//cudaFilterModeLinear and cudaFilterModePoint

  float *h_data;

  //!*R initialization of the particle type interaction
  double cutoffnear = 1.0; 
  int ntypes = 1;
  double *Aij_param;
  double *Bij_param;
  if(true){ //HELIX: Dont use types, helices interact only when their types are different, but in the same way
    Aij_param = new double;
    Bij_param = new double;
    //    float sigma, epsilon;
    //sigma = 2 * lx / float(mx); //READ FROM FILE
    //epsilon = temperature;
    *Aij_param = 0.0;//48.0f * pow(sigma,12)*epsilon;
    *Bij_param = 0.0;// 48.0f * pow(sigma,6)*0.5*epsilon;
  }
  else{    
     ifstream in("LJ.in");
     in>>ntypes;
     Aij_param = new double[ntypes*ntypes];
     Bij_param = new double[ntypes*ntypes];
     for(int i=0; i<ntypes; i++)for(int j=0; j<ntypes; j++)
				  in>>Aij_param[i+ntypes*j];
     
     
     for(int i=0; i<ntypes; i++)for(int j=0; j<ntypes; j++)
				  in>>Bij_param[i+ntypes*j];
     in>>cutoffnear;
   }
      
  //!*R Upload all the information to the GPU
  cudaMalloc((void **)&Aij_paramGPU, ntypes*ntypes*sizeof(double));
  cudaMalloc((void **)&Bij_paramGPU, ntypes*ntypes*sizeof(double));

  cudaMemcpy(Aij_paramGPU, Aij_param, ntypes*ntypes*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(Bij_paramGPU, Bij_param, ntypes*ntypes*sizeof(double), cudaMemcpyHostToDevice);
  
  //HELIX: be carefull with this
  for(int i=1; i<np; i++){
    if(particle_types[i]!=particle_types[i-1]) ntypes++;
  }
  cudaMemcpyToSymbol(ntypesGPU, &ntypes, sizeof(int));
  cudaMemcpyToSymbol(cutoffnearGPU, &cutoffnear, sizeof(double));
  



  int size = 4096;
  h_data = new float[size];
  float r, dr;
  float cutoff2 = cutoff * cutoff;
  dr = cutoff2/float(size);
  r = 0.5 * dr;
  for(int i=0;i<size;i++){
    h_data[i] = functionForceNonBonded1(sqrt(r))/sqrt(r);
    r += dr;
  }
  h_data[size-1] = 0.;
  h_data[0] = 0.;
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  cutilSafeCall( cudaMallocArray( &forceNonBonded1, &channelDesc, size, 1 )); 
  cutilSafeCall( cudaMemcpyToArray( forceNonBonded1, 0, 0, h_data, size*sizeof(float), cudaMemcpyHostToDevice));
  cutilSafeCall( cudaBindTextureToArray( texforceNonBonded1, forceNonBonded1, channelDesc));

  /*
  r = 0.5 * dr;
    for(int i=0;i<size;i++){
    cout << r << " " << h_data[i] << endl;
    r += dr;
    }
    exit(0);
*/
  cout << "INIT FORCE NON-BONDED 1 COMPLETED" << endl;
  delete[] h_data;
  //!*R clean up
  delete[] Aij_param;
  delete[] Bij_param;
return 1;
}
