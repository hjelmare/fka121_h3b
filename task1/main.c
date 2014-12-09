#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nPoints = 11;
  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 * (nPoints - 1)/totalLength;

  double grid[nPoints][nPoints];
  double denseGrid[2*nPoints-1][2*nPoints-1];
  double rho[nPoints][nPoints];
  int x, y;
  double diff, oldValue;
  double tolerance = 0.00001;
  double maxDiff = 2*tolerance;  // just some high number

  for ( x = 0 ; x < nPoints ; x++ ) {
    for ( y = 0 ; y < nPoints ; y++ ) {
      grid[x][y] = 0;
      rho[x][y] = 0;
    }
  }

  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);

  // setting up the problem still belongs here, as does the structuring of the multigrid
  // so what needs to happen here is defining the two grid sizes, running GS on one of them
  // then, we need an interpolating function in mg_func.c to change grids, then we run GS again
  // from here, and that needs to be in a loop

  while(maxDiff > tolerance) {
    maxDiff = GaussSeidel( nPoints, cellLength, grid, rho);
  }

  end = clock();    // we're not really interested in the time to save the data

  IncreaseGridDensity(nPoints, grid, denseGrid);
  DecreaseGridDensity(2*nPoints-1, denseGrid, grid);

  FILE *fGrid = fopen("grid.data","w");

  for ( x = 0 ; x < nPoints ; x++ ) {
    for (  y = 0 ; y < nPoints ; y++ ) {
      fprintf(fGrid,"%e\t",grid[x][y]);
    }
    fprintf(fGrid,"\n");
  }

  printf("Done! (%e s)\n",((double)(end-start) / CLOCKS_PER_SEC));

  return 0;
}
