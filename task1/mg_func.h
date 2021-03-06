#ifndef _mg_func_h
#define _mg_func_h

extern double GaussSeidel(int size, double cellLength, \
    double grid[size][size], double rho[size][size]);
extern double IncreaseGridDensity(int n, double in[n][n], \
    double out[2*n-1][2*n-1]);
extern double DecreaseGridDensity(int n, double in[n][n], \
    double out[n/2+1][n/2+1]);
extern void ComputeResidual(int n, double c, double grid[n][n], \
    double rho[n][n], double res[n][n]);

#endif
