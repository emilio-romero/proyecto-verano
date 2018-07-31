#include "fvmgud.h"

int interfase_Poisson(int nnodos, double *malla, double pinterfase,
    double(*a)(double,double),
    int tipofrontera1, double frontera1, 
    int tipofrontera2, double frontera2){
  gsl_vector *e=gsl_vector_calloc(nnodos-1); //Arriba
  gsl_vector *d=gsl_vector_calloc(nnodos); 
  gsl_vector *f=gsl_vector_calloc(nnodos-1);//Abajo
  gsl_vector *b=gsl_vector_calloc(nnodos); 
  gsl_vector *res=gsl_vector_calloc(nnodos); 
  double h=(malla[nnodos-1]-malla[0])/((double)nnodos);
  double aim1,ai,aip1,wjm1,wj; 
  double hi,him1,hip1; 
  double t1,t2,t3; 
  int manejofronteras=0; 
  for(int i=1;i<nnodos-1;++i){
    aim1=a(malla[i-1],pinterfase);
    ai=a(malla[i],pinterfase);
    aip1=a(malla[i+1],pinterfase);
    wjm1=2.0*aim1*ai/(aim1+ai);
    wj=2.0*ai*aip1/(ai+aip1);
    hi=0.5*(malla[i+1]-malla[i-1]);
    if(i>1){
      him1=0.5*(malla[i]-malla[i-2]);
    }else{
      him1=(malla[1]-malla[0]);
    }
    if(i<nnodos-2){
      hip1=0.5*(malla[i+2]-malla[i]);
    }else{
      hip1=malla[nnodos-1]-malla[nnodos-2];
    }
  t1=-2.0*wjm1/(him1+hi);
  t2=(2.0*wj/(hi+hip1)+2.0*wjm1/(him1+hi)); 
  t3=-2.0*wj/(hi+hip1);
    gsl_vector_set(f,i-1,t1);
    gsl_vector_set(d,i,t2);
    gsl_vector_set(e,i,t3);
  }
  manejofronteras=condiciones_Frontera(tipofrontera1,frontera1,tipofrontera2,frontera2,
      e,d,f,b,nnodos,h);
  if(manejofronteras==-1){
    printf("No hay suficientes condiciones para resolver\n");
  } 
 gsl_linalg_solve_tridiag(d,e,f,b,res); 
  double phit,at=1.0/(pinterfase-0.001*pinterfase+0.001);
  for(int i=0;i<nnodos;++i){
    if(malla[i]<pinterfase){
      phit=at*malla[i];
    }else{
      phit=0.001*at*malla[i]+1-0.001*at;
    }
    printf("%lf %lf %lf\n",malla[i],gsl_vector_get(res,i),phit);
  }
  gsl_vector_free(d);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
return(0);}

int condiciones_Frontera(int tipofrontera1, double frontera1, int tipofrontera2, double frontera2,
    gsl_vector *e, gsl_vector *d, gsl_vector *f,gsl_vector *b, int nnodos,double h){
 if(tipofrontera1==1){
    /*Fronteras tipo Dirichlet*/ 
   gsl_vector_set(e,0,0.0);
   gsl_vector_set(d,0,1.0);
   gsl_vector_set(b,0,frontera1);
  }else if(tipofrontera1==2){
    /*Fronteras tipo Neumann*/ 
    gsl_vector_set(e,0,gsl_vector_get(e,1));
    gsl_vector_set(d,0,gsl_vector_get(d,1)/2.0);
    gsl_vector_set(b,0,gsl_vector_get(b,0)-frontera1/h);
  }else{
    //printf("No hay condiciones suficientes!\n");
    return(-1);
  }

  if(tipofrontera2==1){
    /*Fronteras tipo Dirichlet*/ 
    gsl_vector_set(d,nnodos-1,1.0);
    gsl_vector_set(f,nnodos-2,0.0);
   gsl_vector_set(b,nnodos-1,frontera2);
  }else if(tipofrontera2==2){
    /*Fronteras tipo Neumann*/ 
    gsl_vector_set(f,nnodos-2,gsl_vector_get(f,nnodos-3)); 
    gsl_vector_set(d,nnodos-1,gsl_vector_get(d,nnodos-2)/2.0); 
    gsl_vector_set(b,nnodos-1,gsl_vector_get(b,nnodos-1)-frontera2/h);
  }else{
    //printf("No hay condiciones suficientes!\n");
    return(-1);
  }

return(0);}


double fa1(double x, double pinterfase){
//  if(x<pinterfase){
    return(1.0); 
 /* }
  else{
    return(1.0);
  }*/
}


