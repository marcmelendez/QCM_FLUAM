//Use
//
//radialDistribution.exe fileInput lx ly lz > fileOutput
//
//


#include <iostream>
#include <string>
#include <stdlib.h> //for pseudorandom numbers
#include <time.h>
#include <fstream>
#include <math.h>
using namespace std;


int main(int argc, char* argv[]){
  string word;
  ifstream fileinput(argv[1]);
  double lx = atof(argv[2]);
  double ly = atof(argv[3]);
  int sample = 1;
  if(argc > 4){
    sample  = atoi(argv[4]);
  }
  int numberOfBins = 1000;
  if(argc > 5){
    numberOfBins  = atoi(argv[5]);
  }
  
  getline(fileinput,word);
  cout << word << endl;
  word = word.substr(18,10);
  int np;
  np = atoi(word.c_str());

  double r,dr,t,xij,yij,zij;
  dr = lx / (2.0 * numberOfBins);
  double x[np], y[np], z[np];
  int bin;
  int n=0;
  int hist[numberOfBins];
  for(int i=0;i<numberOfBins;i++)
    hist[i]=0;

  int skip=0;
  int n_count = 0;
  while(!fileinput.eof()){
    fileinput >> t;
    n++;
    for(int i=0;i<np;i++){
      fileinput >> x[i] >> y[i] >> z[i];
    }
    if(n>skip and ((n % sample) == 0)){
      n_count++;
      for(int i=0;i<np-1;i++){
	for(int j=i+1;j<np;j++){
	  xij = x[i] - x[j];
	  xij = xij - (int(xij/lx + 0.5*((xij>0)-(xij<0)))) * lx;
	  yij = y[i] - y[j];
	  yij = yij - (int(yij/ly + 0.5*((yij>0)-(yij<0)))) * ly;
	  
	  r = sqrt(xij*xij + yij*yij); //2D
	  bin = int(r/dr);
	  if(bin<numberOfBins)
	    hist[bin] = hist[bin] + 2;
	}
      }
    }
  }


  
  double rLow, rUp, nIdeal, dens, pi, constant;
  pi = 4 * atan(1);
  dens = np / (lx*ly);
  constant = pi * dens;
  for(int i=0;i<numberOfBins;i++){
    rLow = i * dr;
    rUp = rLow + dr;
    nIdeal = constant * (pow(rUp,2) - pow(rLow,2));
    // cout << (i+0.5)*dr << " " << hist[i]/((n-skip)*np*nIdeal / (1.0 * sample)) << " " << hist[i] << endl;
    cout << (i+0.5)*dr << " " << hist[i]/(n_count*np*nIdeal) << " " << hist[i] << endl;
  }
}


