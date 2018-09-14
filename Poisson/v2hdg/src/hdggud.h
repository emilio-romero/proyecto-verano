#ifndef HDGGUD_H 
#define HDGGUD_H
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
int interfase_Poisson(int nnodos, double *malla, double pinterfase,
    double(*a)(double,double),
    int tipofrontera1, double frontera1,
    int tipofrontera2, double frontera2);

int condiciones_Frontera(int tipofrontera1, double frontera1, int tipofrontera2, double frontera2,
    gsl_matrix *F,gsl_vector *b, int nnodos);
//double fa1(double x,double pinterfase);
#endif 
