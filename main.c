#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main() {
  clock_t end, start = clock();

  int nPoints = 101;
  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 * (nPoints - 1)/totalLength;

  double grid[nPoints][nPoints];
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

  while ( maxDiff > tolerance ) {
    maxDiff = 0;
    for ( x = 1 ; x < nPoints-1 ; x++ ) {
      for ( y = 1 ; y < nPoints-1; y++ ) {
        oldValue = grid[x][y];
        grid[x][y] = 1.0/4.0 * (grid[x+1][y] + grid[x-1][y] + grid[x][y+1] + grid[x][y-1] - pow(cellLength,2)*rho[x][y]);
        diff = fabs(oldValue - grid[x][y]);
        maxDiff = diff > maxDiff ? diff : maxDiff;
      }
    }
  }

  end = clock();    // we're not really interested in the time to save the data

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
