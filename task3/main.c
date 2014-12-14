#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nMaxPoints = 81;
  int nPoints = 11;
  int nFinerPoints = 2*(nPoints - 1) + 1;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nPoints - 1);
  int chargeOffset = chargeSeparation / 2 / cellLength;

  double **grid;
  double **newGrid;
  double **rho;
  
  int x, y;

  FILE *fLog = fopen("log.data","w");

  // initiera grid och rho på grövsta gridsize
  Allocate2dSq(nPoints, &grid);
  Allocate2dSq(nPoints, &rho);
  
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);

  while( nPoints < nMaxPoints){
    // här ska följande hända:
    // kör multigrid
    Multigrid(nPoints, totalLength, grid, rho, fLog);
    // interpolera upp grid
    Allocate2dSq(nFinerPoints, &newGrid);
    IncreaseGridDensity(nPoints, grid, newGrid);
    Free2dSq(nPoints, grid);
    grid = newGrid;
    newGrid = NULL;
    // freea rho
    Free2dSq(nPoints,rho);
    // uppdatera nPoints och alla parametrar som hänger på den (cellLength och chargeOffset t.ex... fler?)
    nPoints = nFinerPoints;
    nFinerPoints = 2*(nPoints-1)+1;
    cellLength = totalLength / (nPoints - 1);
    chargeOffset = chargeSeparation / 2 / cellLength;
    // skapa ny rho
    Allocate2dSq(nPoints, &rho);
    rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
    rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
    // var det allt?
  }

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
