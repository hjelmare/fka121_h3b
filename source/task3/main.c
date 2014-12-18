#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mg_func.h"

int main() {
  // Settings to change
  int nMaxPoints = 81;
  int nMinPoints = 11;
  FILE *fLog = fopen("log81.data","w");
  FILE *fGrid = fopen("grid81.data","w");
  double tolerance = 0.00001;

  int nPoints = nMinPoints;
  int nFinerPoints = 2*(nPoints - 1) + 1;
  
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

  // Initialize on the coarsest grid
  Allocate2dSq(nPoints, &grid);
  Allocate2dSq(nMaxPoints, &oldGrid);
  Allocate2dSq(nMaxPoints, &saveGrid);
  Allocate2dSq(nPoints, &rho);
  
  // Set up source distribution
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);

  // Iterate until our best approximation changes less than tolerance
  maxDiff = 2*tolerance;  // Just to get the loop started
  while(maxDiff > tolerance) {
    maxDiff = 0;
    nPoints = nMinPoints;
    nFinerPoints = 2*(nPoints-1)+1;

    // Start from coarsest grid, keep going to desired size
    while( nPoints <= nMaxPoints){
      // Run multigrid
      Multigrid(nPoints, totalLength, grid, rho, fLog);
      
      // If finest grid, save data now, before interpolating
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

      // Free stuff we don't need anymore
      Free2dSq(nPoints,rho);

      // Update nPoints and derived numbers
      nPoints = nFinerPoints;
      nFinerPoints = 2*(nPoints-1)+1;
      cellLength = totalLength / (nPoints - 1);
      chargeOffset = chargeSeparation / 2 / cellLength;

      // Create a new rho at current grid size
      Allocate2dSq(nPoints, &rho);
      rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = \
      -1.0 / pow(cellLength,2);
      rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = \
      1.0 / pow(cellLength,2);
    }
    
    // Since we're at the finest grid, revert last change of nPoints
    nPoints = (nPoints-1)/2 + 1;

    // Compare previous solution to new one
    for(x = 0 ; x < nPoints ; x++ ) {
      for( y = 0 ; y < nPoints ; y++) {
        diff = fabs(oldGrid[x][y] - grid[x][y]);
        oldGrid[x][y] = grid[x][y];
        maxDiff = maxDiff > diff ? maxDiff : diff;
      }
    }
  }

  // Save solution
  int nPlotPoints = nMaxPoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",saveGrid[x][y]);

    }
    fprintf(fGrid,"\n");
  }

  printf("Done\n");

  return 0;
}
