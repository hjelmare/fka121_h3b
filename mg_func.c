#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

void MakeGaussRandom(double* x, double* y)
{
  double u = *x;
  double v = *y;
  
  *x = sqrt( -2*log(u) ) * cos(2*PI*v);
  *y = sqrt( -2*log(u) ) * sin(2*PI*v);

  return;
}
