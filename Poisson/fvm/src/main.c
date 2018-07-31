#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <time.h>
#include <string.h>
#include "mesh.h"
#include "fvmgud.h"
int main(int argc, char *argv[]){
  int n; 
  if(argc>1) n=atoi(argv[1]); else n=10; 

  double *lamalla=genera_malla1D(n,0,1);
  /*
  for(int i=0;i<n;++i){
    printf("Nodo %d, coordenada %lf\n",i,lamalla[i]);
  }*/
  interfase_Poisson(n,lamalla,0.5,fa1,1,0.0,1,1.0);
  //printf("Sus diferencias finitas han terminado.\nTenga un buen dia!\n");
return(0);}
