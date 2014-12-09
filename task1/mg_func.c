#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define PI 3.141592653589

double GaussSeidel(int nPoints, double cellLength, double grid[nPoints][nPoints], double rho[nPoints][nPoints]) {
  double maxDiff = 0;
  double oldValue, diff;
  int x,y;

  for ( x = 1 ; x < nPoints - 1 ; x++ ) {
    for ( y = 1 ; y < nPoints - 1 ; y++ ) {
      oldValue = grid[x][y];
      grid[x][y] = 1.0/4.0 * (grid[x+1][y] + grid[x-1][y] + grid[x][y+1] + grid[x][y-1] - pow(cellLength,2)*rho[x][y]);
      diff = fabs(oldValue - grid[x][y]);
      maxDiff = diff > maxDiff ? diff : maxDiff;
    }
  }

  return(maxDiff);
}

void IncreaseGridDensity(int nPoints, double inGrid[nPoints][nPoints], double outGrid[2*nPoints-1][2*nPoints-1]) {
  int x, y;

  // the exactly matching points
  for ( x = 0 ; x < nPoints ; x++ ) {
    for ( y = 0 ; y < nPoints ; y++ ) {
      outGrid[2*x][2*y] = inGrid[x][y];
    }
  }
  // points with four nearest neighbours
  for ( x = 1 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 1 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/4.0 * (inGrid[x/2][y/2] + inGrid[x/2+1][y/2] + inGrid[x/2][y/2+1] + inGrid[x/2+1][y/2+1]);
    }
  }
  // the half points in y
  for ( x = 0 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 1 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/2.0 * (inGrid[x/2][(y-1)/2] + inGrid[x/2][(y+1)/2]);
    }
  } // and in x
  for ( x = 1 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 0 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/2.0 * (inGrid[(x-1)/2][y/2] + inGrid[(x+1)/2][y/2]);
    }
  }

  return;
}

void DecreaseGridDensity(int nPoints, double inGrid[nPoints][nPoints], double outGrid[nPoints/2+1][nPoints/2+1]) {
  int x,y;

  for ( x = 1 ; x < nPoints/2+1 ; x++ ) {
    for ( y = 1 ; y < nPoints/2+1 ; y++ ) {
      outGrid[x][y] = 1.0/4.0 * inGrid[2*x][2*y] \
                    + 1.0/8.0 * (inGrid[2*x-1][2*y] + inGrid[2*x+1][2*y] + inGrid[2*x][2*y-1] + inGrid[2*x][2*y+1]) \
                    + 1.0/16.0 * (inGrid[2*x-1][2*y-1] + inGrid[2*x+1][2*y-1] + inGrid[2*x+1][2*y+1] + inGrid[2*x-1][2*y+1]);
    }
  }
  
  return;
}

void ComputeResidual(int nPoints, double cellLength, double grid[nPoints][nPoints], \
                          double rho[nPoints][nPoints], double residual[nPoints][nPoints]) {
  double oldValue, diff;
  int x,y;
  double tolerance = 0.00001;
  double maxDiff = 2*tolerance;

  for ( x = 1 ; x < nPoints - 1 ; x++ ) {
    for ( y = 1 ; y < nPoints - 1 ; y++ ) {
      residual[x][y] = rho[x][y] \
                     - (grid[x-1][y] - 2*grid[x][y] + grid[x+1][y]) / pow(cellLength,2) \
                     - (grid[x][y-1] - 2*grid[x][y] + grid[x][y+1]) / pow(cellLength,2);
    }
  }

  return;
}

