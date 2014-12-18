#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mg_func.h"

int main() {

  int nCoarsePoints = 11;
  int nFinePoints = 2*(nCoarsePoints - 1) +1; // must be N = 2*(n-1)+1

  double chargeSeparation = 0.2;
  double totalLength = 1;
  double fineCellLength = totalLength / (nFinePoints - 1);
  double coarseCellLength = totalLength / (nCoarsePoints - 1);
  int fineChargeOffset = chargeSeparation / 2 * \
  (nFinePoints - 1)/totalLength;

  double fineGrid[nFinePoints][nFinePoints];
  double oldFineGrid[nFinePoints][nFinePoints];
  double fineRho[nFinePoints][nFinePoints];
  double fineResidual[nFinePoints][nFinePoints];
  double coarseResidual[nCoarsePoints][nCoarsePoints];
  double coarseRho[nCoarsePoints][nCoarsePoints];
  double coarseError[nCoarsePoints][nCoarsePoints];
  double fineError[nFinePoints][nFinePoints];
  
  int x, y, i, z,j;
  double diff;
  int nPresmooth = 2;   // "a few"
  int nPostsmooth = 2;
  double innerTolerance = 0.00001;
  double innerMaxDiff = 2*innerTolerance;  // high enough to start looping
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

  // Setting up the source distribution
  fineRho[nFinePoints / 2 + fineChargeOffset][nFinePoints / 2 ] = \
  -1.0 / pow(fineCellLength,2);
  fineRho[nFinePoints / 2 - fineChargeOffset][nFinePoints / 2 ] = \
  1.0 / pow(fineCellLength,2);
  // End of init part

  while( outerMaxDiff > outerTolerance) {
    // Keep a copy of the old solution to check convergence
    for( x = 0 ; x < nFinePoints ; x++ ) {
      for( y = 0 ; y < nFinePoints ; y++ ) {
        oldFineGrid[x][y] = fineGrid[x][y];
      }
    }
    
    for( i = 0 ; i < nPresmooth ; i++ ) {
      GaussSeidel( nFinePoints, fineCellLength, fineGrid, fineRho);
    }

    ComputeResidual(nFinePoints, fineCellLength, fineGrid, \
    fineRho, fineResidual);

    DecreaseGridDensity(nFinePoints, fineResidual, coarseResidual);

    // Solve residual eq on coarse grid "exactly" i.e. to within tolerance
    innerMaxDiff = 2*innerTolerance;
    while (innerMaxDiff > innerTolerance ) {
      innerMaxDiff = GaussSeidel( nCoarsePoints, coarseCellLength, \
      coarseError, coarseResidual);
    }

    IncreaseGridDensity(nCoarsePoints, coarseError, fineError);

    // Apply interpolated correction
    for ( x = 0 ; x < nFinePoints ; x++ ) {
      for ( y =0 ; y < nFinePoints ; y++ ) {
        fineGrid[x][y] = fineGrid[x][y] + fineError[x][y];
      }
    }

    for ( i = 0 ; i < nPostsmooth ; i++ ) {
      GaussSeidel( nFinePoints, fineCellLength, fineGrid, fineRho);
    }

    // Check size of largest change made this iteration
    outerMaxDiff = 0;
    for( x = 0 ; x < nFinePoints ; x++ ) {
      for( y = 0 ; y < nFinePoints ; y++ ) {
        diff = fabs(oldFineGrid[x][y] - fineGrid[x][y]);
        outerMaxDiff = diff > outerMaxDiff ? diff : outerMaxDiff;
      }
    }
  }

  // Save solution
  FILE *fGrid = fopen("grid.data","w");

  int nPlotPoints = nFinePoints;
  for ( x = 0 ; x < nPlotPoints ; x++ ) {
    for (  y = 0 ; y < nPlotPoints ; y++ ) {
      fprintf(fGrid,"%e\t",fineGrid[x][y]);
    }
    fprintf(fGrid,"\n");
  }

  printf("Done!\n");

  return 0;
}
