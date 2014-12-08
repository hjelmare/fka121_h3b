#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

void GaussSeidel(int nPoints, double cellLength, double grid[nPoints][nPoints], double rho[nPoints][nPoints], double tolerance) {
  double maxDiff = 2*tolerance; // value is not important, just to get the loop started
  double oldValue, diff;
  int x,y;

  while( maxDiff > tolerance ) {
    maxDiff = 0;
    for ( x = 1 ; x < nPoints - 1 ; x++ ) {
      for ( y = 1 ; y < nPoints - 1 ; y++ ) {
        oldValue = grid[x][y];
        grid[x][y] = 1.0/4.0 * (grid[x+1][y] + grid[x-1][y] + grid[x][y+1] + grid[x][y-1] - pow(cellLength,2)*rho[x][y]);
        diff = fabs(oldValue - grid[x][y]);
        maxDiff = diff > maxDiff ? diff : maxDiff;
      }
    }
  }




  return;
}
