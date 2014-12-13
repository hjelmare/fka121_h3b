#ifndef _fft_func_h
#define _fft_func_h

extern void Multigrid(int size, double cellLength, double** grid, double** source);
extern double GaussSeidel(int size, double cellLength, double** grid, double** rho);
extern void IncreaseGridDensity(int n, double** in, double** out);
extern void DecreaseGridDensity(int n, double** in, double** out);
extern void ComputeResidual(int n, double c, double** grid, double** rho, double** res);

#endif
