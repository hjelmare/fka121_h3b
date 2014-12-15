#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nMaxPoints = 1281;
  int nMinPoints = 11;
  int nPoints = nMinPoints;
  int nFinerPoints = 2*(nPoints - 1) + 1;
  
  FILE *fLog = fopen("log.data","w");
  FILE *fGrid = fopen("grid.data","w");
  double tolerance = 0.00001;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 / cellLength;

  double maxDiff;
  double diff;

  double **grid;
  double **newGrid;
  double **oldGrid;
  double **saveGrid;
  double **rho;
  
  int x, y;
  
  fprintf(fLog, "Points\tGS its\n");

  // initiera grid och rho på grövsta gridsize
  Allocate2dSq(nPoints, &grid);
  Allocate2dSq(nMaxPoints, &oldGrid);
  Allocate2dSq(nMaxPoints, &saveGrid);
  Allocate2dSq(nPoints, &rho);
  
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);

  maxDiff = 2*tolerance;
  while(maxDiff > tolerance) {
    maxDiff = 0;
    nPoints = nMinPoints;
    nFinerPoints = 2*(nPoints-1)+1;

    while( nPoints <= nMaxPoints){
      // Run multigrid
      Multigrid(nPoints, totalLength, grid, rho, fLog);

      if ( nPoints == nMaxPoints ) {
        for ( x = 0 ; x < nPoints ; x++ ) {
          for ( y = 0 ; y < nPoints ; y++ ) {
            saveGrid[x][y] = grid[x][y];
          }
        }
      }

      // Interpolate the result to higher gridsize
      Allocate2dSq(nFinerPoints, &newGrid);
      IncreaseGridDensity(nPoints, grid, newGrid);
      Free2dSq(nPoints, grid);
      grid = newGrid;
      newGrid = NULL;

      // Free memory
      Free2dSq(nPoints,rho);

      // Update nPoints and derived numbers
      nPoints = nFinerPoints;
      nFinerPoints = 2*(nPoints-1)+1;
      cellLength = totalLength / (nPoints - 1);
      chargeOffset = chargeSeparation / 2 / cellLength;

      // Create a new rho at current grid size
      Allocate2dSq(nPoints, &rho);
      rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
      rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);

    }
//    fprintf(fLog, "%d\t0\n",nPoints);
    nPoints = (nPoints-1)/2 + 1;

    for(x = 0 ; x < nPoints ; x++ ) {
      for( y = 0 ; y < nPoints ; y++) {
        diff = fabs(oldGrid[x][y] - grid[x][y]);
        oldGrid[x][y] = grid[x][y];
        maxDiff = maxDiff > diff ? maxDiff : diff;
      }
    }
  }

  int nPlotPoints = nMaxPoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",saveGrid[x][y]);

    }
    fprintf(fGrid,"\n");

  }

  end = clock();

  printf("Done! (%e s)\n",((double)(end-start) / CLOCKS_PER_SEC));

  return 0;
}
