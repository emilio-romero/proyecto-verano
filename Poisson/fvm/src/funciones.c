#include "funciones.h"

double rosenbrock(gsl_vector *x){
  double aux=0.0,xip1, xi;
  int n=(int)x->size; 
  for(int i=0;i<(n-1);++i){
    xip1=gsl_vector_get(x,i+1);
    xi=gsl_vector_get(x,i);
    aux=aux+100.0*gsl_pow_2(xip1-xi*xi) + gsl_pow_2(1.0-xi); 
  }
return(aux);}

int grosenbrock(gsl_vector *x, gsl_vector *g){
  int n=(int)x->size;
  double xip1,xi,xim1;
  double t1,t2,t3; 
  xi=gsl_vector_get(x,0);
  xip1=gsl_vector_get(x,1);
  gsl_vector_set(g,0,-400.0*xi*(xip1-xi*xi)-2.0*(1-xi));
  for(int i=1;i<n-1;++i){
    xim1=gsl_vector_get(x,i-1);
    xi=gsl_vector_get(x,i);
    xip1=gsl_vector_get(x,i+1);
    t1=200.0*(xi-xim1*xim1);
    t2=-400.0*xi*(xip1-xi*xi);
    t3=-2.0*(1.0-xi);
    gsl_vector_set(g,i,t1+t2+t3);
  }
  xim1=gsl_vector_get(x,n-2);
  xi=gsl_vector_get(x,n-1);
  gsl_vector_set(g,n-1,200.0*(xi-xim1*xim1));
return(1);}

double sphere(gsl_vector *x){
  double aux=0,xi; 
  int n=(int)x->size; 
  for(int i=0;i<n;++i){
   xi=gsl_vector_get(x,i);
   aux+=xi*xi;
  }
return(aux);}

int gsphere(gsl_vector *x, gsl_vector *g){
int n=(int)x->size; 
double xi; 
for(int i=0;i<n;i++){
  xi=gsl_vector_get(x,i);
  gsl_vector_set(g,i,2.0*xi);
}
return(1);}

double sumsquares(gsl_vector *x){
  double aux=0,xi; 
  int n=(int)x->size; 
  for(int i=0;i<n;++i){
   xi=gsl_vector_get(x,i);
   aux+=xi*xi*(double)(i+1);
  }
return(aux);}

int gsumsquares(gsl_vector *x, gsl_vector *g){
int n=(int)x->size; 
double xi; 
for(int i=0;i<n;i++){
  xi=gsl_vector_get(x,i);
  gsl_vector_set(g,i,2.0*xi*(double)(i+1));
}
return(1);}




double rastrigin(gsl_vector *x){
  double aux=0,xi; 
  int n=(int)x->size; 
  for(int i=0;i<n;++i){
    xi=gsl_vector_get(x,i);
    aux+=xi*xi-10.0*gsl_sf_cos(2.0*xi*M_PI);
  }
return(aux+10*n);}

int grastrigin(gsl_vector *x, gsl_vector *g){
int n=(int)x->size; 
double xi; 
for(int i=0;i<n;i++){
  xi=gsl_vector_get(x,i);
  gsl_vector_set(g,i,2.0*xi+20.0*M_PI*gsl_sf_sin(2.0*xi*M_PI));
}
return(1);}

double booth(gsl_vector *x){
  double aux=0,xi,y; 
  xi=gsl_vector_get(x,0);
  y=gsl_vector_get(x,1);
  aux=gsl_pow_2(xi+2.0*y-7.0)+gsl_pow_2(2.0*xi+y-5);
return(aux);}

int gbooth(gsl_vector *x, gsl_vector *g){
  double xi,yi; 
  xi=gsl_vector_get(x,0);
  yi=gsl_vector_get(x,1);
  gsl_vector_set(g,0,2.0*(xi+2.0*yi-7)+4.0*(2.0*xi+yi-5)); 
  gsl_vector_set(g,1,4.0*(xi+2.0*yi-7)+2.0*(2.0*xi+yi-5)); 
return(1);}

double himmelblau(gsl_vector *x){
  double aux=0,xi,y; 
  xi=gsl_vector_get(x,0);
  y=gsl_vector_get(x,1);
  aux=gsl_pow_2(xi*xi+y-11.0)+gsl_pow_2(xi+y*y-7);
return(aux);}

int ghimmelblau(gsl_vector *x, gsl_vector *g){
  double xi,yi; 
  xi=gsl_vector_get(x,0);
  yi=gsl_vector_get(x,1);
  gsl_vector_set(g,0,4.0*xi*(xi*xi+yi-11.0)+2.0*(xi+yi*yi-7.0)); 
  gsl_vector_set(g,1,2.0*(xi*xi+yi-11.0)+4.0*yi*(xi+yi*yi-7.0)); 
return(1);}

double wood(gsl_vector* x){
double aux=0,x1,x2,x3,x4; 
x1=gsl_vector_get(x,0);
x2=gsl_vector_get(x,1);
x3=gsl_vector_get(x,2);
x4=gsl_vector_get(x,3);
aux+=100.0*gsl_pow_2(x1*x1-x2)+gsl_pow_2(x1-1.0)+gsl_pow_2(x3-1.0); 
aux+=90.0*gsl_pow_2(x3*x3-x4)+10.1*(gsl_pow_2(x2-1.0)+gsl_pow_2(x4-1.0))+19.8*(x2-1.0)*(x4-1.0);

return(aux);}

int gwood(gsl_vector* x,gsl_vector *g){
double aux=0,x1,x2,x3,x4; 
x1=gsl_vector_get(x,0);
x2=gsl_vector_get(x,1);
x3=gsl_vector_get(x,2);
x4=gsl_vector_get(x,3);

gsl_vector_set(g,0,400.0*x1*(x1*x1-x2)+2.0*(x1-1.0));
gsl_vector_set(g,1,-200.0*(x1*x1-x2)+20.2*(x1-1.0)+19.8*(x4-1.0));
gsl_vector_set(g,2,2.0*(x3-1.0)+360.0*x3*(x3*x3-x4));
gsl_vector_set(g,3,-180.0*(x3*x3-x4)+20.2*(x4-1.0)+19.8*(x2-1.0));

return(aux);}

