#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nPoints = 81;
  int coarsestGrid = 11;

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double cellLength = totalLength / (coarsestGrid - 1);
  int chargeOffset = chargeSeparation / 2 * (coarsestGrid - 1)/totalLength;

  double **grid;
  double **rho;
  grid = (double**) malloc(coarsestGrid * sizeof(double*));
  rho = (double**) malloc(coarsestGrid * sizeof(double*));
  
  int x, y, i;

  // Initializing
  for ( x = 0 ; x < nPoints ; x++ ) {
    grid[x] = (double*) calloc(nPoints, sizeof(double));
    rho[x] = (double*) calloc(nPoints, sizeof(double));
  }
  
  rho[nPoints / 2 + chargeOffset][nPoints / 2 ] = 1.0 / pow(cellLength,2);
  rho[nPoints / 2 - chargeOffset][nPoints / 2 ] = -1.0 / pow(cellLength,2);
  // End of init part

//  Multigrid(nPoints, totalLength, grid, rho);

// NYTT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  i = 0;
  while((nPoints-1)/(2*i) > 11){
    i+=i;
  }
  nbrOfMultigridCallings = i;

  i=0;
  while( (coarsestGrid-1)*pow(2,i) < nPoints){

  
  
  // Calculates Rho at a finer grid
  cellLength = totalLength / (coarsestGrid*pow(2,i) - 1);
  chargeOffset = chargeSeparation / 2 * (coarsestGrid*pow(2,i) - 1)/totalLength;
  
  Free2DSq(coarsestGrid*pow(2,i), grid);  
  Free2DSq(coarsestGrid*pow(2,i), rho);  
  double **grid;
  double **rho;
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
