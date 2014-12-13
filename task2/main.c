#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "mg_func.h"

int main() {
  clock_t end, start = clock();

  int nCoarsePoints = 11;
  int nFinePoints = 21; // must be N = 2*(n-1)+1

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double fineCellLength = totalLength / (nFinePoints - 1);
  double coarseCellLength = totalLength / (nCoarsePoints - 1);
  int fineChargeOffset = chargeSeparation / 2 * (nFinePoints - 1)/totalLength;

  double fineGrid[nFinePoints][nFinePoints];
  double oldFineGrid[nFinePoints][nFinePoints];
  double fineRho[nFinePoints][nFinePoints];
  double fineResidual[nFinePoints][nFinePoints];
  double coarseResidual[nCoarsePoints][nCoarsePoints];
  double coarseRho[nCoarsePoints][nCoarsePoints];
  double coarseError[nCoarsePoints][nCoarsePoints];
  double fineError[nFinePoints][nFinePoints];
  
  int x, y, i;
  double diff;
  int nPresmooth = 2;
  int nPostsmooth = 2;
  double innerTolerance = 0.00001;
  double innerMaxDiff = 2*innerTolerance;  // just some high number
  double outerTolerance = 0.00001;
  double outerMaxDiff = 2*outerTolerance;

  // Initializing
  for ( x = 0 ; x < nFinePoints ; x++ ) {
    for ( y = 0 ; y < nFinePoints ; y++ ) {
      fineGrid[x][y] = 0;
      fineRho[x][y] = 0;
      fineResidual[x][y] = 0;
      fineError[x][y] = 0;
    }
  }
  for ( x = 0 ; x < nCoarsePoints ; x++ ) {
    for ( y = 0 ; y < nCoarsePoints ; y++ ) {
      coarseResidual[x][y] = 0;
      coarseRho[x][y] = 0;
      coarseError[x][y] = 0;
    }
  }

  fineRho[nFinePoints / 2 + fineChargeOffset][nFinePoints / 2 ] = 1.0 / pow(fineCellLength,2);
  fineRho[nFinePoints / 2 - fineChargeOffset][nFinePoints / 2 ] = -1.0 / pow(fineCellLength,2);
  // End of init part

  while( outerMaxDiff > outerTolerance) {
    for( x = 0 ; x < nFinePoints ; x++ ) {
      for( y = 0 ; y < nFinePoints ; y++ ) {
        oldFineGrid[x][y] = fineGrid[x][y];
      }
    }

    for( i = 0 ; i < nPresmooth ; i++ ) {
      GaussSeidel( nFinePoints, fineCellLength, fineGrid, fineRho);
    }

    ComputeResidual(nFinePoints, fineCellLength, fineGrid, fineRho, fineResidual);

    DecreaseGridDensity(nFinePoints, fineResidual, coarseResidual);

    innerMaxDiff = 2*innerTolerance;
    while (innerMaxDiff > innerTolerance ) {
      innerMaxDiff = GaussSeidel( nCoarsePoints, coarseCellLength, coarseError, coarseResidual);
    }

    IncreaseGridDensity(nCoarsePoints, coarseError, fineError);

    for ( x = 0 ; x < nFinePoints ; x++ ) {
      for ( y =0 ; y < nFinePoints ; y++ ) {
        fineGrid[x][y] = fineGrid[x][y] + fineError[x][y];
      }
    }

    for ( i = 0 ; i < nPostsmooth ; i++ ) {
      GaussSeidel( nFinePoints, fineCellLength, fineGrid, fineRho);
    }

    outerMaxDiff = 0;
    for( x = 0 ; x < nFinePoints ; x++ ) {
      for( y = 0 ; y < nFinePoints ; y++ ) {
        diff = fabs(oldFineGrid[x][y] - fineGrid[x][y]);
        outerMaxDiff = diff > outerMaxDiff ? diff : outerMaxDiff;
      }
    }
}

  FILE *fGrid = fopen("grid.data","w");

  int nPlotPoints = nFinePoints;
//  nPlotPoints = nCoarsePoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",fineGrid[x][y]);
    }
    fprintf(fGrid,"\n");
  }

  end = clock();    // we're not really interested in the time to save the data

  printf("Done! (%e s)\n",((double)(end-start) / CLOCKS_PER_SEC));

  return 0;
}
