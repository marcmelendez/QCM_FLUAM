// Filename: saveFluidFinalConfiguration.cu
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


#include "header.h"
#include "headerOtherFluidVariables.h"
#include "fluid.h"
#include "cells.h"
#include "particles.h"
//#include "visit_writer.h"
//#include "visit_writer.c"
#include "hydroAnalysis.h"

/***** Added by Adolfo *****/
#include <sstream>
/**************************/

bool saveFluidFinalConfiguration(){
  ofstream file;
  
  string savefile;
  savefile = outputname + ".fluidFinalConfiguration";
  file.open(savefile.c_str());
  file.precision(15);
  file << "Number of cells " << ncells << endl;
  file << "mx my mz" << endl;
  file << mx << " " << my << " " << mz << endl;
  file << "lx ly lz" << endl;
  file << lx << " " << ly << " " << lz << endl;
  file << "Density " << densfluid << endl;
  file << "Pressure parameters " << pressurea0 << " " << pressurea1 << " " << pressurea2 << endl;
  file << "Thermostat " << thermostat << endl;
  file << "Temperature " << temperature << endl;
  file << "Cells properties " << endl;
  file << "Density " << endl;
  file << "Velocity" << endl;


  int i;
  for(int j=0;j<ncells;j++){
    if(particlesWall){
      //go to index tanking into account ghost cells
      int fx = j % mx;
      int fy = (j % (my*mx)) / mx;
      int fz = j / (mx*my);
      fy++;
      i = fx + fy*mxt + fz*mxt*myt;
    }
    else{
      i = j;
    }
    file << cDensity[i] << endl;
    if(incompressibleBinaryMixture || incompressibleBinaryMixtureMidPoint)
      file << c[i] << endl;
    file << cvx[i] << " " << cvy[i] << " " << cvz[i] << endl;
  }
    
  file.close();


  return 1;

}


/****** Subroutine added by Adolfo *******/
 bool saveFluid_column(int option, long long step){
  ofstream file;

 // Sending a number as a stream into output
  // string
  ostringstream str1;
  str1 << step;
  
  // the str() coverts number into string
  string step_strg = str1.str();
  
  string savefile;
  savefile = outputname + ".fluid_profile" + step_strg;
  file.open(savefile.c_str());
  file.precision(15);

  int i;
  int fx, fy, fz;
  for(int j=0;j<ncells;j++){
    fx = j % mx;
    fy = (j % (my*mx)) / mx;
    fz = j / (mx*my);
    i = fx + fy*mxt + fz*mxt*myt;
    //    file << fx << " " << fy << " " << fz << endl;
    if ( (fx == 0) && (fy == 0) )
      file << fz << " " << cvx[i]  << endl;
  }
    
  file.close();

  return 1;

}


/******** Added by Adolfo **********/
bool saveFluid_gradv_average(long long step){
  ofstream file;

 // Sending a number as a stream into output
  // string
  ostringstream str1;
  str1 << step;
  
  // the str() coverts number into string
  string step_strg = str1.str();
  
  string savefile;
  savefile = outputname + ".gradv" + step_strg;
  file.open(savefile.c_str());
  file.precision(15);

  double gradv[11];

  for(int j=0;j<11;j++)
    gradv[j] = 0.0;

 double dh = -0.1;
  
  long long counter = 0;
  for(int j=0;j<ncells;j++){
    //go to index taking into account ghost cells
    int fz = j / (mx*my);
    double invdx = double(mx)/lx;
    if (fz == mz/2+1) //close to the wall (but not exactly on it)
      {
	counter = counter + 1;
	int fx = j % mx;
	int fy = (j % (my*mx)) / mx;
	// int i = fx + fy*mxt + fz*mxt*myt;
	int vec1 = fx + fy*mxt + (fz+1)*mxt*myt; 
	int vec2 = fx + fy*mxt + (fz+2)*mxt*myt; 
	int vec3 = fx + fy*mxt + (fz+3)*mxt*myt; 
  
	/*** We look for a polynomia v(h) = A3 h^3 + A2 h^2 + A1 h + A0.
         The interpolation points are at h = 0, 1, 2, 3;
	 The velocity gradient at h = 0 is A1 * invdx     *****/
	 double A3 = 1.0/2.0 * (1.0/3.0 * (cvx[vec3] - cvx[j]) + cvx[vec1] - cvx[vec2]);
	 double A2 = 1.0/2.0 * (fz*(cvx[j] - cvx[vec3] - 3.0*cvx[vec1] + 3.0*cvx[vec2]) + 2.0*cvx[j] -5.0*cvx[vec1] + 4.0*cvx[vec2] -cvx[vec3]);
	 double A1 = 1.0/2.0 * (fz*fz*(cvx[vec3] - cvx[j] + 3.0*cvx[vec1] - 3.0*cvx[vec2]) + 2.0*fz*(-2.0*cvx[j] + 5.0*cvx[vec1] - 4.0*cvx[vec2] + cvx[vec3]) + (- 11.0/3.0*cvx[j] + 6.0*cvx[vec1] - 3.0*cvx[vec2] + 2.0/3.0*cvx[vec3] ) );
//	 double A0 = 1.0/2.0 * (1.0/3.0 * (fz*fz*fz + 6.0*fz*fz + 11.0*fz + 6.0)*cvx[j] +fz*(cvx[vec2]*(fz*fz+4.0*fz +3.0) - 1.0/3.0*cvx[vec3]*(fz*fz + 3.0*fz + 2.0) - cvx[vec1]*(fz*fz+5.0*fz+6.0)));
	
	for(int k=0;k<11;k++){
 	  double zcoord = fz + dh * (double) k;
	  gradv[k] = gradv[k] + (3.0*A3*zcoord*zcoord + 2.0*A2*zcoord + A1) * invdx;
	}
      }
    
      }

      
    for(int j=0;j<11;j++)
    {
      gradv[j] = gradv[j] / counter;
      double h = 0.5*lz + 1 + dh * (double) j;
      file << step << " " << h << " " << gradv[j] << endl; 
    }

  file.close();
    return 1;

}

/******** Added by Adolfo **********/
bool saveFluid_gradv(long long step){
  ofstream file;

 // Sending a number as a stream into output
  // string
  ostringstream str1;
  str1 << step;
  
  // the str() coverts number into string
  string step_strg = str1.str();
  
  string savefile;
  savefile = outputname + ".gradv" + step_strg;
  file.open(savefile.c_str());
  file.precision(15);

  long long counter = 0;
  for(int j=0;j<ncells;j++){
    //go to index taking into account ghost cells
    int fz = j / (mx*my);
    double invdx = double(mx)/lx;
    if (fz == mz/2+1) //close to the wall (but not exactly on it)
      {
	counter = counter + 1;
	int fx = j % mx;
	int fy = (j % (my*mx)) / mx;
	// int i = fx + fy*mxt + fz*mxt*myt;
	int vec1 = fx + fy*mxt + (fz+1)*mxt*myt; 
	int vec2 = fx + fy*mxt + (fz+2)*mxt*myt; 
	int vec3 = fx + fy*mxt + (fz+3)*mxt*myt; 
  
	/*** We look for a polynomia v(h) = A3 h^3 + A2 h^2 + A1 h + A0.
         The interpolation points are at h = 1, 2, 3, 4;
	 The velocity gradient at h = 0 is A1 * invdx     *****/
	 double A3 = 1.0/2.0 * (1.0/3.0 * (cvx[vec3] - cvx[j]) + cvx[vec1] - cvx[vec2]);
	 double A2 = 1.0/2.0 * (fz*(cvx[j] - cvx[vec3] - 3.0*cvx[vec1] + 3.0*cvx[vec2]) + 2.0*cvx[j] -5.0*cvx[vec1] + 4.0*cvx[vec2] -cvx[vec3]);
	 double A1 = 1.0/2.0 * (fz*fz*(cvx[vec3] - cvx[j] + 3.0*cvx[vec1] - 3.0*cvx[vec2]) + 2.0*fz*(-2.0*cvx[j] + 5.0*cvx[vec1] - 4.0*cvx[vec2] + cvx[vec3]) + (- 11.0/3.0*cvx[j] + 6.0*cvx[vec1] - 3.0*cvx[vec2] + 2.0/3.0*cvx[vec3] ) );
//	 double A0 = 1.0/2.0 * (1.0/3.0 * (fz*fz*fz + 6.0*fz*fz + 11.0*fz + 6.0)*cvx[j] +fz*(cvx[vec2]*(fz*fz+4.0*fz +3.0) - 1.0/3.0*cvx[vec3]*(fz*fz + 3.0*fz + 2.0) - cvx[vec1]*(fz*fz+5.0*fz+6.0)));
	
	 double zcoord = fz - 1.0 + 0.115;
	 double gradv  = (3.0*A3*zcoord*zcoord + 2.0*A2*zcoord + A1) * invdx;
	 file << fx << " " << fy << " " << gradv << endl;
      }
    
  }

  file.close();
  return 1;

}

/****** Subroutine added by Adolfo *******/
 bool saveFluid(int option, long long step){
  ofstream file;

 // Sending a number as a stream into output
    // string
    ostringstream str1;
    str1 << step;
 
    // the str() coverts number into string
    string step_strg = str1.str();

  string savefile;
  savefile = outputname + ".fluid" + step_strg;
  file.open(savefile.c_str());
  file.precision(15);

  int i;
  int fx, fy, fz;
  for(int j=0;j<ncells;j++){
    if(particlesWall){
      //go to index tanking into account ghost cells
      fx = j % mx;
      fy = (j % (my*mx)) / mx;
      fz = j / (mx*my);
      fy++;
      i = fx + fy*mxt + fz*mxt*myt;
    }
    else{
      i = j;
    }
      file << fx << " " << fy << " " << fz << " " << " " << cvx[i] << " " << cvy[i] << " " << cvz[i] << " " << cdydvx[i] << " " << cDensity[i] << endl;
//    file << fx << " " << fy << " " << fz << " " << " " << cvx[i] << " " << cvy[i] << " " << cvz[i]  << endl;
  }
    
  file.close();

  return 1;

}
