#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mcbsp.h"
#include "mcbsp-affinity.h"

#define SZDBL (sizeof(double))
#define SZINT (sizeof(int))
#define TRUE (1)
#define FALSE (0)
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

double  *vecallocd(int n);
int     *vecalloci(int n);
double  **matallocd(int m, int n);
void    vecfreed(double *pd);
void    vecfreei(int *pi);
void    matfreed(double **ppd);
