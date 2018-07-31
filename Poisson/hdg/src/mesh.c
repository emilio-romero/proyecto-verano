#include "mesh.h" 

double *genera_malla1D(int nnodos, double x0, double xf){
  double h=(xf-x0)/(double)(nnodos-1);  
  double *mimalla=(double*)malloc(nnodos*sizeof(double));
  for(int i=0;i<nnodos;++i){
    mimalla[i]=x0+h*(double)i;
  } 
return(mimalla);}
