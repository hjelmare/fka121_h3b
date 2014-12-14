#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nMaxPoints = 81;
  int nCoarsestPoints = 11;
  int nCoarsePoints, nFinePoints;
  int nPoints;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (nCoarsestPoints - 1);
  int chargeOffset = chargeSeparation / 2 * (nCoarsestPoints - 1)/totalLength;

  double** grid;
  double** rho;
  
  int x, y, i;

  FILE *logFile = fopen("log.data","w");

  i=0;
  nPoints = nCoarsestPoints * pow(2,i);
  // initiera grid och rho på grövsta gridsize
printf("%d\n", nPoints);
  Allocate2Dsq(nPoints, grid);
  Allocate2Dsq(nPoints, rho);
  
  printf("%e", nPoints);
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);

  //Kör den första multigriden
//  Multigrid(nCoarsestPoints, totalLength, grid, rho, logFile);

// NYTT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
  while( nPoints < nMaxPoints){
    i++;
    nCoarsePoints = nCoarsestPoints*pow(2,i-1);
    nFinePoints = nCoarsestPoints*pow(2,i);
    // här ska följande hända:
    // kör multigrid
    // interpolera upp grid
    // uppdatera nPoints och alla parametrar som hänger på den (cellLength och chargeOffset t.ex... fler?)
    // freea rho
    // skapa ny rho
    // var det allt?
    
    //skapa minnes plats, används längre ner
    double** fineGrid;
    Allocate2Dsq(nFinePoints, fineGrid);

    // Gör om till finare grid
    IncreaseGridDensity(nCoarsePoints, grid, fineGrid);

    //Uppdatera rho & grid
    Free2DSq(nCoarsePoints, rho); //free memory space 
    Free2DSq(nCoarsePoints, grid);  
    double** grid;
    double** rho;
    Allocate2Dsq(nFinePoints, grid);
    Allocate2Dsq(nFinePoints, rho);
    rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
    rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
    
    for(x=0; x<nFinePoints; x++){
      for(y=0; y<nFinePoints; y++){
        grid[x][y] = fineGrid[x][y];
      }
    }
    
    //Kör multigrid
    Multigrid(nFinePoints, totalLength, grid, rho, logFile);
    

    cellLength = totalLength / (nFinePoints - 1);
    chargeOffset = (chargeSeparation / 2) / cellLength;
  
printf("nPoints:\t nCource = %e\t nFine =%e\n", nCoarsePoints, nFinePoints);  
  
  }

// SLUT NYTT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  FILE *fGrid = fopen("grid.data","w");
  int nPlotPoints = nMaxPoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",grid[x][y]);

    }
    fprintf(fGrid,"\n");

  }

  end = clock();
*/
  printf("Done! (%e s)\n",((double)(end-start) / CLOCKS_PER_SEC));

  return 0;
}
