// Filename: particles.h
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

#ifdef GLOBALS_PARTICLES
#define EXTERN_PARTICLES
#else
#define EXTERN_PARTICLES extern
#endif

EXTERN_PARTICLES int np;

typedef struct {
    double r[3], v[3], f[3];
    //mass;
} particle;

//EXTERN_PARTICLES particle* p;
EXTERN_PARTICLES int mxrcut, myrcut, mzrcut;
EXTERN_PARTICLES int dtn;
EXTERN_PARTICLES double rcut, rcut2;
EXTERN_PARTICLES double sigma, zeta, rho;
EXTERN_PARTICLES double mass;
EXTERN_PARTICLES double sigmap;
EXTERN_PARTICLES int* head;
EXTERN_PARTICLES int* list;
//EXTERN_PARTICLES int** neighbor1;
//EXTERN_PARTICLES int** neighbor2;
EXTERN_PARTICLES double *rxParticle, *ryParticle, *rzParticle;
//!*R particle_type information for the CPU
EXTERN_PARTICLES int *particle_types;
EXTERN_PARTICLES double *vxParticle, *vyParticle, *vzParticle;
EXTERN_PARTICLES double *vxParticleI, *vyParticleI, *vzParticleI;
EXTERN_PARTICLES double volumeParticle;
EXTERN_PARTICLES bool loadparticles;
EXTERN_PARTICLES bool loadcolors;
EXTERN_PARTICLES double cutoff;
EXTERN_PARTICLES int numNeighbors;
EXTERN_PARTICLES bool setVolumeParticle;
EXTERN_PARTICLES int greenStart, greenEnd;
//EXTERN_PARTICLES int nparticles;

//NEW bonded forces
EXTERN_PARTICLES bool bondedForces;
EXTERN_PARTICLES bool bondedForcesVersion;
EXTERN_PARTICLES string bondedForcesFile;
EXTERN_PARTICLES int nbondsParticleParticle;
EXTERN_PARTICLES int nbondsParticleFixedPoint;
EXTERN_PARTICLES int *bondsParticleParticle, *bondsParticleParticleOffset;
EXTERN_PARTICLES int *bondsParticleFixedPoint, *bondsParticleFixedPointOffset;
EXTERN_PARTICLES int *bondsIndexParticleParticle;
EXTERN_PARTICLES int *bondsIndexParticleFixedPoint;
EXTERN_PARTICLES double* kSpringParticleParticle;
EXTERN_PARTICLES double* r0ParticleParticle;
EXTERN_PARTICLES double* kSpringParticleFixedPoint;
EXTERN_PARTICLES double* r0ParticleFixedPoint;
EXTERN_PARTICLES double* rxFixedPoint;
EXTERN_PARTICLES double* ryFixedPoint;
EXTERN_PARTICLES double* rzFixedPoint;

EXTERN_PARTICLES bool threeBondedForces;
EXTERN_PARTICLES string threeBondedForcesFile;
EXTERN_PARTICLES int NbondsThreeParticle;
EXTERN_PARTICLES int *threebondList;
EXTERN_PARTICLES double *threekSprings;
EXTERN_PARTICLES double *threer0Springs;
EXTERN_PARTICLES int *threeNbonds;
EXTERN_PARTICLES int *threeisinbonds;
EXTERN_PARTICLES int *threeCumulativeIndex;




EXTERN_PARTICLES string LJParameterFile;
EXTERN_PARTICLES bool LJParameterFileProvided;
