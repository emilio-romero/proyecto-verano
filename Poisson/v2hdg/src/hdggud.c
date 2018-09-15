#include "hdggud.h"

int interfase_Poisson(int nnodos, double *malla, double pinterfase,
    double(*a)(double,double),
    int tipofrontera1, double frontera1, 
    int tipofrontera2, double frontera2){
  int nelem=nnodos-1;
  gsl_matrix *F=gsl_matrix_calloc(2*(nelem-2),2*(nelem-2));
  gsl_vector *b=gsl_vector_calloc(2*(nelem-2)); 
  gsl_vector *resmedios=gsl_vector_calloc(2*(nelem-2)); 
  double ai,aip1,bi=0.5,bip1=0.5; 
  double dxi,dxip1;
  int nodoactl, nodoactlp1; 
  int manejofronteras=0; 
 
  for(int i=1;i<nelem-2;++i){
    nodoactl=2*i; 
    nodoactlp1=2*i+1; 
    dxi=0.5*(malla[i+1]-malla[i-1]);
    dxip1=0.5*(malla[i+2]-malla[i]);
    ai=1.0/(2.0+dxi/(a(malla[i],pinterfase)));
    aip1=1.0/(2.0+dxip1/(a(malla[i+1],pinterfase)));
   /*============ Matriz tridiagonal por bloques ==============*/ 
    /*Submatriz 1 (izquierda)*/
    gsl_matrix_set(F,nodoactl,nodoactl-2,ai); 
    gsl_matrix_set(F,nodoactl,nodoactl-1,bi); 
    gsl_matrix_set(F,nodoactlp1,nodoactl-2,ai); 
    gsl_matrix_set(F,nodoactlp1,nodoactl-1,bi); 
    /*Submatriz 2 (central)*/ 
    gsl_matrix_set(F,nodoactl,nodoactl,-2.0+ai+aip1);
    gsl_matrix_set(F,nodoactl,nodoactlp1,bi-bip1);
    gsl_matrix_set(F,nodoactlp1,nodoactl,ai-aip1);
    gsl_matrix_set(F,nodoactlp1,nodoactlp1,-2.0+bi+bip1);
    /*Submatriz 3 (derecha)*/
    gsl_matrix_set(F,nodoactl,nodoactlp1+1,aip1);
    gsl_matrix_set(F,nodoactl,nodoactlp1+2,-bip1);
    gsl_matrix_set(F,nodoactlp1,nodoactlp1+1,-aip1);
    gsl_matrix_set(F,nodoactlp1,nodoactlp1+2,bip1);
  }
  /*===== Punto x_3/2    =====*/
    dxi=(malla[1]-malla[0]);
    dxip1=0.5*(malla[2]-malla[0]);
    ai=1.0/(2.0+dxi/(a(malla[0],pinterfase)));
    aip1=1.0/(2.0+dxip1/(a(malla[1],pinterfase)));
    /*Submatriz 2 (central)*/ 
    gsl_matrix_set(F,0,0,ai+aip1-2.0);
    gsl_matrix_set(F,0,1,1.0-ai-bi);
    gsl_matrix_set(F,1,0,ai-aip1);//ai-aip1);
    gsl_matrix_set(F,1,1,bip1-ai-1.0);//-2.0+bi+bip1);
    /*Submatriz 3 (derecha)*/
    gsl_matrix_set(F,0,2,aip1);//aip1);
    gsl_matrix_set(F,0,3,-bip1);//-bip1);
    gsl_matrix_set(F,1,2,-aip1);//-aip1);
    gsl_matrix_set(F,1,3,bip1);//bip1);
    /*Submatriz 3 (derecha)*/
    /*gsl_matrix_set(F,0,4,-aip1);
    gsl_matrix_set(F,0,5,-bip1);
    gsl_matrix_set(F,1,4,-aip1);
    gsl_matrix_set(F,1,5,-bip1);
*/
  /*===== Punto x_n-1/2  =====*/
    dxi=0.5*(malla[nnodos-1]-malla[nnodos-3]);
    dxip1=(malla[nnodos-1]-malla[nnodos-2]);
    ai=1.0/(2.0+dxi/(a(malla[nnodos-2],pinterfase)));
    aip1=1.0/(2.0+dxip1/(a(malla[nnodos-1],pinterfase)));
    /*Submatriz 1 (mas izquierda)*/
    /*gsl_matrix_set(F,2*nnodos-4,2*nnodos-8,-ai); 
    gsl_matrix_set(F,2*nnodos-4,2*nnodos-7,bi); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-8,ai); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-7,-bi); */
    /*Submatriz 1 (izquierda)*/
    gsl_matrix_set(F,2*nnodos-4,2*nnodos-6,ai); 
    //gsl_matrix_set(F,2*nnodos-4,2*nnodos-5,bi-bip1); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-6,ai-aip1); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-5,3.0*bip1-bi); 
    /*Submatriz 2 (central)*/ 
    gsl_matrix_set(F,2*nnodos-4,2*nnodos-4,ai-1.0);
    //gsl_matrix_set(F,2*nnodos-4,2*nnodos-3,-bi);
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-4,ai);
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-3,3.0*bi-2.0);

  /*===== Agregar condiciones de frontera =====*/
  manejofronteras=condiciones_Frontera(tipofrontera1,frontera1,tipofrontera2,frontera2,
      F,b,nnodos);
  if(manejofronteras==-1){
    printf("No hay suficientes condiciones para resolver\n");
  }
/*
  //Imprimir matriz 
  for(int i=0;i<2*(nnodos-1);++i){
    for(int j=0;j<2*(nnodos-1);++j){
      printf("%4.2lf ",gsl_matrix_get(F,i,j));
    }
    printf("\n");
  }
  //
*/

  /*===Resolver el sistema de ecuaciones ===*/ 
  int s; 
  gsl_permutation *p=gsl_permutation_alloc(2*(nnodos-1));
  gsl_linalg_LU_decomp(F,p,&s);
  gsl_linalg_LU_solve(F,p,b,resmedios);
  /*========================================*/ 
  /*
  double phit,at=1.0/(pinterfase-0.001*pinterfase+0.001);
  */
  double resi; 
  printf("%lf %lf\n",malla[0],frontera1);
  for(int i=0;i<nnodos-2;++i){
    /*
    if(malla[i]<pinterfase){
      phit=at*malla[i];
    }else{
      phit=0.001*at*malla[i]+1-0.001*at;
    }
*/
    resi=0.5*(gsl_vector_get(resmedios,2*i+1)+gsl_vector_get(resmedios,2*i+3));
    printf("%lf %lf\n",malla[i+1],resi);
  }
  printf("%lf %lf\n",malla[nnodos-1],frontera2);
 /* 
  printf("Medios\n");

  for(int i=0;i<2*(nnodos-1);++i){
    if(i%2!=0){
      printf("%lf\n",gsl_vector_get(resmedios,i));
    }
  }
  printf("Fin medios\n");
  */
  gsl_matrix_free(F);
  gsl_vector_free(b);
  gsl_vector_free(resmedios);
return(0);}

int condiciones_Frontera(int tipofrontera1, double frontera1, int tipofrontera2, double frontera2,
    gsl_matrix *F, gsl_vector *b, int nnodos){
 if(tipofrontera1==1){
    /*===Fronteras tipo Dirichlet===*/
    /*Submatriz 2 (central)*/
    /*for(int i=0;i<2*(nnodos-1);++i){
      gsl_matrix_set(F,i,1,0);
    } */
    gsl_matrix_set(F,1,0,0.0);
    gsl_matrix_set(F,1,1,1.0);
    /*Submatriz 3 (derecha)*/
    gsl_matrix_set(F,1,2,0.0);
    gsl_matrix_set(F,1,3,0.0);
    /*Submatriz 3 (mas derecha)*/
    gsl_matrix_set(F,1,4,0.0);
    gsl_matrix_set(F,1,5,0.0);
   
    gsl_vector_set(b,1,frontera1);
  }
 else if(tipofrontera1==2){
   /*===Fronteras tipo Neumann===*/
   /*Submatriz 2 (central)*/
   gsl_matrix_set(F,0,0,1);
   gsl_matrix_set(F,0,1,0);
   /*Submatriz 3 (derecha)*/
   gsl_matrix_set(F,0,2,0);
   gsl_matrix_set(F,0,3,0);
   gsl_vector_set(b,0,frontera1);
  }else{
    return(-1);
  }

  if(tipofrontera2==1){
    /*===Fronteras tipo Dirichlet===*/ 
    /*Submatriz 1 (mas izquierda)*/
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-8,0.0); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-7,0.0); 
    /*Submatriz 1 (izquierda)*/
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-6,0.0); 
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-5,0.0); 
    /*Submatriz 2 (central)*/ 
    /*for(int i=0;i<2*(nnodos-1);++i){
      gsl_matrix_set(F,i,2*nnodos-3,0);
    } */
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-4,0.0);
    gsl_matrix_set(F,2*nnodos-3,2*nnodos-3,1.0);
   
    gsl_vector_set(b,2*nnodos-3,frontera2);
  }else if(tipofrontera2==2){
    /*===Fronteras tipo Neumann===*/
    /*Submatriz 1 (izquierda)*/
    gsl_matrix_set(F,nnodos-4,nnodos-6,0); 
    gsl_matrix_set(F,nnodos-4,nnodos-5,0); 
    /*Submatriz 2 (central)*/ 
    gsl_matrix_set(F,2*nnodos-4,2*nnodos-4,1);
    gsl_matrix_set(F,2*nnodos-4,2*nnodos-3,0);
    
    gsl_vector_set(b,2*nnodos-4,frontera2);
  }else{
    return(-1);
  }
return(0);}



