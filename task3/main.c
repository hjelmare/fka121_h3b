#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nMaxPoints = 81;
  int nCoarsestPoints = 11;
  int nCoarcePoints, nFinePoints;
  int nPoints;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (coarsestGrid - 1);
  int chargeOffset = chargeSeparation / 2 * (coarsestGrid - 1)/totalLength;

  double **grid;
  double **rho;
  
  int x, y, i;

  FILE *logFile = fopen("log.data","w");

  i=0;
  nPoints = nCoarsestPoints * pow(2,i);
  // initiera grid och rho på grövsta gridsize
  for ( x = 0 ; x < nPoints ; x++ ) {
    grid[x] = (double*) calloc(nPoints, sizeof(double));
    rho[x] = (double*) calloc(nPoints, sizeof(double));
  }
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);

  Multigrid(coarsestGrid, totalLength, grid, rho, fLog);

// NYTT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  while( nPoints < nMaxPoints){
    i++;
    nCourseGridPoints = nCoarsestPoints*pow(2,i-1);
    nFinePoints = nCoarsestPoints*pow(2,i);
    // här ska följande hända:
    // kör multigrid
    // interpolera upp grid
    // uppdatera nPoints och alla parametrar som hänger på den (cellLength och chargeOffset t.ex... fler?)
    // freea rho
    // skapa ny rho
    // var det allt?
    
    //skapa minnes plats, används längre ner
    fineGrid = (double**) malloc(nFinePoints * sizeof(double*));
    for ( x = 0 ; x < finerGridPoints ; x++ ) {
      finerGrid[x] = (double*) calloc(nFinePoints, sizeof(double));
    }

    // Gör om till finare grid
    IncreaseGridDensity(nCoarsePoints, grid, fineGrid);
    
    //Kör multigrid
    Multigrid(nPoints, totalLength, fineGrid, rho, logFile);
    

    cellLength = totalLength / (nPoints - 1);
    chargeOffset = (chargeSeparation / 2) / cellLength;
  
    Free2DSq(nPoints, grid);  
    Free2DSq(nPoints, rho);  
    
    grid = (double**) malloc(nPoints * sizeof(double*));
    rho = (double**) malloc(nPoints * sizeof(double*));

    // Initializing
    for ( x = 0 ; x < nPoints ; x++ ) {
      grid[x] = (double*) calloc(nPoints, sizeof(double));
      rho[x] = (double*) calloc(nPoints, sizeof(double));
    }
  
    rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
    rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
    // end Calculate Rho
  }

// SLUT NYTT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
