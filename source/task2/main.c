#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mg_func.h"

int main() {
  // Setting to be changed for each run
  int nPoints = 81;
  FILE *fGrid = fopen("grid81w.data","w");
  int gamma = 2;
  FILE *fLog = fopen("log81w.data","w");

  fprintf(fLog, "Points\t GS its\n");

  double tolerance = 0.00001;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 / cellLength;
  double maxDiff = 2*tolerance;
  double diff;

  double **grid;
  double **oldGrid;
  double **rho;
  Allocate2dSq(nPoints, &grid);
  Allocate2dSq(nPoints, &oldGrid);
  Allocate2dSq(nPoints, &rho);
  
  int x, y, i;
  // Set up charge distribution
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  
  // Iterate multigrid to convergence
  while ( maxDiff > tolerance ) {
    Multigrid(gamma, nPoints, totalLength, grid, rho, fLog);

    // Compare before and after
    maxDiff = 0;
    for ( x = 0 ; x < nPoints; x++ ) {
      for ( y = 0 ; y < nPoints; y++ ) {
        diff = fabs( grid[x][y] - oldGrid[x][y] );
        maxDiff = maxDiff > diff ? maxDiff : diff;
        
        oldGrid[x][y] = grid[x][y];
      }
    }
  }
  
  fclose(fLog);

  // Save solution
  int nPlotPoints = nPoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",grid[x][y]);

    }
    fprintf(fGrid,"\n");
  }

  Free2dSq(nPoints, grid);
  Free2dSq(nPoints, oldGrid);
  Free2dSq(nPoints, rho);
  
  printf("Done!\n");

  return 0;
}
