#include "gradienteconjugado.h"

int GC_NTPA(gsl_vector *x0, double epsi, double t0, double(*f)(gsl_vector*),
    int(*g)(gsl_vector*, gsl_vector*),int maxiter,gsl_vector *sol){
  int n=(int)x0->size;
  gsl_vector *xk=gsl_vector_calloc(n);
  gsl_vector *dk=gsl_vector_calloc(n);
  gsl_vector *gk=gsl_vector_calloc(n);
  gsl_vector *xkm1=gsl_vector_calloc(n);
  gsl_vector *dkm1=gsl_vector_calloc(n);
  gsl_vector *gkm1=gsl_vector_calloc(n);
  gsl_vector *bd=gsl_vector_calloc(n);
  gsl_vector *dely=gsl_vector_calloc(n);
  gsl_vector *y=gsl_vector_calloc(n);
  gsl_vector *s=gsl_vector_calloc(n);
  //fin de creacion de vectores
  int cont=0;
  double alphak,normg,beta,delta; 
  double ak, tk,ss,yy,sy; 
  double dty,gty,gts;
  gsl_vector_memcpy(xk,x0);
  g(xk,gk);
  gsl_vector_memcpy(dk,gk);
  gsl_vector_scale(dk,-1.0);
  tk=t0;
  while(cont<maxiter){
    alphak=Wolfe_conditions(xk,dk,gk,f,g); //Cambiar por condiciones de Wolfe
    /*Calcular x_{k+1}*/
    gsl_vector_memcpy(xkm1,dk);
    gsl_vector_scale(xkm1,alphak);
    gsl_vector_add(xkm1,xk);
    /*Calcular g_{k+1}*/
    g(xkm1,gkm1);
    /*Criterio de paro*/
    normg=gsl_blas_dnrm2(gkm1);
    if(normg<epsi){
      gsl_vector_memcpy(xk,xkm1);
      break;
    }
    /*Calculo de y_k y s_k*/
    gsl_vector_memcpy(y,gkm1);
    gsl_vector_sub(y,gk);
    gsl_vector_memcpy(s,xkm1);
    gsl_vector_sub(s,xk);
    /*Calcular a_k*/
    ss=gsl_blas_dnrm2(s);
    yy=gsl_blas_dnrm2(y);
    gsl_blas_ddot(s,y,&sy);
    ak=ss*ss*yy*yy/(sy*sy);
    /*Calcular d_{k+1}*/
    gsl_vector_memcpy(dkm1,gkm1);
    gsl_vector_scale(dkm1,-1.0);
    //Calcular beta
    /**/
    gsl_blas_ddot(dk,y,&dty);
    gsl_blas_ddot(gkm1,y,&gty);
    gsl_blas_ddot(gkm1,s,&gts);
    beta=(tk*gty-gts)/dty;
    if(beta<0) beta=0.0; 
    /**/
    gsl_vector_memcpy(bd,dk);
    gsl_vector_scale(bd,beta);
    //Calcular delta
    /**/
    delta=tk*gts/sy; 
    /**/
    gsl_vector_memcpy(dely,y);
    gsl_vector_scale(dely,delta);
    gsl_vector_add(dkm1,bd);
    gsl_vector_add(dkm1,dely);
    /*Actualizacion de x_k, g_k y d_k*/
    gsl_vector_memcpy(dk,dkm1);
    gsl_vector_memcpy(xk,xkm1);
    gsl_vector_memcpy(gk,gkm1);
    tk=1.0/(1+ak);
   /* 
    if(sy/(yy*yy)<tk){
      tk=sy/(yy*yy);
    }
   /* */
    cont++;
  
  }
    printf("k:%d f:%g |g|=%g \n",cont,f(xk),normg);
  gsl_vector_memcpy(sol,xk);
  //Liberacion de memoria
  gsl_vector_free(xk);
  gsl_vector_free(dk);
  gsl_vector_free(gk);
  gsl_vector_free(xkm1);
  gsl_vector_free(dkm1);
  gsl_vector_free(gkm1);
  gsl_vector_free(bd);
  gsl_vector_free(dely);
  gsl_vector_free(y);
  gsl_vector_free(s);
return(1);}

double Wolfe_conditions(gsl_vector *xk, gsl_vector *dk, gsl_vector *gk,
    double(*f)(gsl_vector*),int(*g)(gsl_vector*,gsl_vector*)){
  double a=0.0,c1=0.0001,c2=0.1,t=1.0,b=1000000; 
  double fk,ft,gtd,gttd;
  int n=xk->size; 
  gsl_vector *xt=gsl_vector_alloc(n);
  gsl_vector *gt=gsl_vector_alloc(n);
  gsl_vector_memcpy(xt,dk);
  gsl_blas_ddot(gk,dk,&gtd);
  fk=f(xk);
  int cont=0;
  do{
    gsl_vector_scale(xt,t);
    gsl_vector_add(xt,xk);
    ft=f(xt);
    g(xt,gt);
    gsl_blas_ddot(gt,dk,&gttd);
    if(ft>fk+c1*t*gtd){
      b=t; 
      t=0.5*(a+b);
    }
    else if(gttd<c2*gtd){
      a=t; 
      if(b>=1000000) t=2*a; 
      else t=0.5*(a+b);
    }
    else{
      break;
    }
    /*reset*/
    gsl_vector_memcpy(xt,dk);
    cont++;
  }while(cont<10);
  gsl_vector_free(xt);
  gsl_vector_free(gt);
return(t);}


