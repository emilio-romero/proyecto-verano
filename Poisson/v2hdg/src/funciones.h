#ifndef FUNCIONES_H
#define FUNCIONES_H
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_trig.h>
double rosenbrock(gsl_vector *x);
int grosenbrock(gsl_vector *x, gsl_vector *g);
double sphere(gsl_vector *x); 
int gsphere(gsl_vector *x, gsl_vector *g);
double sumsquares(gsl_vector *x); 
int gsumsquares(gsl_vector *x, gsl_vector *g);

double rastrigin(gsl_vector *x);
int grastrigin(gsl_vector *x, gsl_vector *g);
double booth(gsl_vector *x);
int gbooth(gsl_vector *x,gsl_vector *g);
double himmelblau(gsl_vector *x);
int ghimmelblau(gsl_vector *x,gsl_vector *g);
double wood(gsl_vector *x);
int gwood(gsl_vector *x,gsl_vector *g);
#endif 
