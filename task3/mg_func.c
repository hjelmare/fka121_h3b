#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mg_func.h"
#define PI 3.141592653589
#define COARSESTPOINTS 11

void Allocate2dSq(int size, double ***array)
{
  int i;
  *array = (double**) malloc(size * sizeof(double));
  for ( i = 0 ; i < size ; i++ ) {
    (*array)[i] = (double*) calloc(size, sizeof(double));
  }
  
  return;
}

void Free2dSq(int size, double **array)
{
  int i;
  for ( i = 0 ; i < size ; i++) {
    free(array[i]);
    array[i] = NULL;
  }
  free(array);
  array = NULL;

  return;
}


void Multigrid(int nPoints, double totalLength, double** grid, double** source, FILE *fLog)
{
  int i,j,k,x,y;
  int gamma = 2;
  int nPresmooth = 2, nPostsmooth = 2;
  int nCoarsePoints = nPoints/2+1;
  double cellLength = totalLength / (double) (nPoints-1);
  double** residual;
  double** coarseResidual;
  double** v;
  double** fineV;
  double tolerance = 0.00001;
  double maxDiff = 2*tolerance;

  residual = (double**) malloc(nPoints * sizeof(double*));
  fineV = (double**) malloc(nPoints * sizeof(double*));
  coarseResidual = (double**) malloc(nCoarsePoints * sizeof(double*));
  v = (double**) malloc(nCoarsePoints * sizeof(double*));

  for (i = 0 ; i < nPoints ; i++ ) {
    residual[i] = (double*) calloc(nPoints, sizeof(double));
    fineV[i] = (double*) calloc(nPoints, sizeof(double));
  }
  for ( i = 0 ; i < nCoarsePoints ; i++ ) {
    coarseResidual[i] = (double*) calloc(nCoarsePoints, sizeof(double));
    v[i] = (double*) calloc(nCoarsePoints, sizeof(double));
  }

  if ( nPoints <= COARSESTPOINTS ) {
    while (maxDiff > tolerance ) {
      maxDiff = GaussSeidel(nPoints, totalLength, grid, source);
    }
    fprintf(fLog, "%d\n",nPoints);
  } else {
    for( i = 0 ; i < nPresmooth ; i++ ) {
      GaussSeidel( nPoints, totalLength, grid, source);
    }

    ComputeResidual(nPoints, totalLength, grid, source, residual);

    DecreaseGridDensity(nPoints, residual, coarseResidual);

    for ( i = 0 ; i < nCoarsePoints ; i++ ) {
      for ( j = 0 ; j < nCoarsePoints ; j++ ) {
        v[i][j] = 0;
      }
    }
    
    fprintf(fLog,"%d\n", nPoints);
    for ( k = 0 ; k < gamma ; k++) {
      Multigrid(nCoarsePoints, totalLength, v, coarseResidual, fLog);
      fprintf(fLog,"%d\n", nPoints);
    }
    IncreaseGridDensity(nCoarsePoints, v, fineV);

    for ( x = 0 ; x < nPoints ; x++ ) {
      for ( y =0 ; y < nPoints ; y++ ) {
        grid[x][y] = grid[x][y] + fineV[x][y];
      }
    }

    for ( i = 0 ; i < nPostsmooth ; i++ ) {
      GaussSeidel( nPoints, totalLength, grid, source);
    }
  }

  FILE *fV = fopen("res.data","w");
  for ( i = 0 ; i < nPoints ; i++ ) {
    for ( j = 0 ; j < nPoints ; j++ ) {
      fprintf(fV, "%e\t",fineV[i][j]);
    }
    fprintf(fV,"\n");
  }
  fclose(fV);

  return; // grid is "returned" as an argument...
}


double GaussSeidel(int nPoints, double totalLength, double** grid, double** rho) {
  double maxDiff = 0;
  double cellLength = totalLength/ (double) (nPoints-1);
  double oldValue, diff;
  int x,y;

  for ( x = 1 ; x < nPoints - 1 ; x++ ) {
    for ( y = 1 ; y < nPoints - 1 ; y++ ) {
      oldValue = grid[x][y];
      grid[x][y] = 1.0/4.0 * (grid[x+1][y] + grid[x-1][y] + grid[x][y+1] + grid[x][y-1] - pow(cellLength,2)*rho[x][y]);
      diff = fabs(oldValue - grid[x][y]);
      maxDiff = diff > maxDiff ? diff : maxDiff;
    }
  }

  return(maxDiff);
}

void IncreaseGridDensity(int nPoints, double** inGrid, double** outGrid) {
  int x, y;
  // the exactly matching points
  for ( x = 0 ; x < nPoints ; x++ ) {
    for ( y = 0 ; y < nPoints ; y++ ) {
      outGrid[2*x][2*y] = inGrid[x][y];
    }
  }
  // points with four nearest neighbours
  for ( x = 1 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 1 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/4.0 * (inGrid[x/2][y/2] + inGrid[x/2+1][y/2] + inGrid[x/2][y/2+1] + inGrid[x/2+1][y/2+1]);
    }
  }
  // the half points in y
  for ( x = 0 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 1 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/2.0 * (inGrid[x/2][(y-1)/2] + inGrid[x/2][(y+1)/2]);
    }
  } // and in x
  for ( x = 1 ; x < 2*nPoints-1 ; x += 2 ) {
    for ( y = 0 ; y < 2*nPoints-1 ; y += 2 ) {
      outGrid[x][y] = 1.0/2.0 * (inGrid[(x-1)/2][y/2] + inGrid[(x+1)/2][y/2]);
    }
  }

  return;
}

void DecreaseGridDensity(int nPoints, double** inGrid, double** outGrid) {
  int x,y;

  for ( x = 1 ; x < nPoints/2 ; x++ ) {
    for ( y = 1 ; y < nPoints/2 ; y++ ) {
      outGrid[x][y] = 1.0/4.0 * inGrid[2*x][2*y] \
                    + 1.0/8.0 * (inGrid[2*x-1][2*y] + inGrid[2*x+1][2*y] + inGrid[2*x][2*y-1] + inGrid[2*x][2*y+1]) \
                    + 1.0/16.0 * (inGrid[2*x-1][2*y-1] + inGrid[2*x+1][2*y-1] + inGrid[2*x+1][2*y+1] + inGrid[2*x-1][2*y+1]);
    }
  }
  
  return;
}

void ComputeResidual(int nPoints, double totalLength, double** grid, \
                          double** rho, double** residual) {
  double oldValue, diff;
  int x,y;
  double cellLength = totalLength / (double) (nPoints-1);
  double tolerance = 0.00001;
  double maxDiff = 2*tolerance;

  for ( x = 1 ; x < nPoints - 1 ; x++ ) {
    for ( y = 1 ; y < nPoints - 1 ; y++ ) {
      residual[x][y] = rho[x][y] \
                     - (grid[x-1][y] - 2*grid[x][y] + grid[x+1][y]) / pow(cellLength,2) \
                     - (grid[x][y-1] - 2*grid[x][y] + grid[x][y+1]) / pow(cellLength,2);
    }
  }

  return;
}

