#ifndef GRADIENTE_CONJUGADO_H
#define GRADIENTE_CONJUGADO_H
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <math.h>
int GC_NTPA(gsl_vector *x0, double epsi, double t0, double(*f)(gsl_vector*),
    int(*g)(gsl_vector*, gsl_vector*),int maxiter, gsl_vector *sol);

double Wolfe_conditions(gsl_vector *xk, gsl_vector *dk, gsl_vector *gk,
    double(*f)(gsl_vector*),int(*g)(gsl_vector*,gsl_vector*));

#endif 
