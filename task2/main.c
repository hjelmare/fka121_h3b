#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nPoints = 21;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 * (nPoints - 1)/totalLength;

  double **grid;
  double **rho;
  grid = (double**) malloc(nPoints * sizeof(double*));
  rho = (double**) malloc(nPoints * sizeof(double*));
  
  int x, y, i;

  // Initializing
  for ( x = 0 ; x < nPoints ; x++ ) {
    grid[x] = (double*) calloc(nPoints, sizeof(double));
    rho[x] = (double*) calloc(nPoints, sizeof(double));
  }
  
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
  // End of init part

  Multigrid(nPoints, cellLength, grid, rho);

  FILE *fGrid = fopen("grid.data","w");


  int nPlotPoints = nPoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",grid[x][y]);

    }
    fprintf(fGrid,"\n");

  }

  end = clock();

  printf("Done! (%e s)\n",((double)(end-start) / CLOCKS_PER_SEC));

  return 0;
}
