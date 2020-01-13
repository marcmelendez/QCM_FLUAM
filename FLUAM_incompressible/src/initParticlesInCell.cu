// Filename: initParticlesInCell.cu
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


__global__ void initParticlesInCell(int* partincellX, 
				    int* partincellY, 
				    int* partincellZ,
				    int* countparticlesincellX, 
				    int* countparticlesincellY, 
				    int* countparticlesincellZ,
				    int* countPartInCellNonBonded, 
				    int* partInCellNonBonded,
				    particlesincell *pc){

  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>0) return;   
  pc->countparticlesincellX = countparticlesincellX;
  pc->countparticlesincellY = countparticlesincellY;
  pc->countparticlesincellZ = countparticlesincellZ;
  pc->partincellX = partincellX;
  pc->partincellY = partincellY;
  pc->partincellZ = partincellZ;
  pc->countPartInCellNonBonded = countPartInCellNonBonded;
  pc->partInCellNonBonded = partInCellNonBonded;


  return;

}



__global__ void initParticlesInCellOmega(int* partincell,
                                         int* countparticlesincell,
                                         particlesincell *pc){


  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>0) return;

  pc->countparticlesincell = countparticlesincell;
  pc->partincell = partincell;

  return;

}

//!*R
/*!
*@brief Links particle_type pt GPU struct with all particle type information stored in the GPU
*@param *pt pointer to GPU's particle_type struct
*@param *types GPU array storing particle types (from 0 to np)
*@param *Aij First parameters of the force, DIM ntypes*ntypes
. It is a matrix stored in a 1D array
*/
__global__ void initParticleTypes(particle_type *pt, int *types, double *Aij, double *Bij){
  int i = blockDim.x * blockIdx.x + threadIdx.x;
  if(i>0) return;

  pt->types = types;
  pt->Aij_param = Aij;
  pt->Bij_param = Bij;

}

