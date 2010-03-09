#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


/* The function in this file is supposed to get a matrix 
// and return the maximum of each row. The out put is a vector.*/


// The index (i,j) in matrix M having n rows. i and j start from 0.
#define INDEX(M,n,i,j) (M)[(i) + (j)*(n)]


SEXP maximum_of_rows(SEXP matrix_from_R){

	double *matrix;
	int columnNum, rowNum, i, j;
	SEXP matrixDim, returnVal;
	double max;

	// To set the type from R to C.
	PROTECT(matrix_from_R = coerceVector(matrix_from_R, REALSXP));	

	// Computing the values of the indices.
	matrixDim = getAttrib(matrix_from_R, R_DimSymbol);
    rowNum = INTEGER(matrixDim)[0];
    columnNum = INTEGER(matrixDim)[1];
	
	// Getting values of indeces in C, matrix is a double*.
	matrix = REAL(matrix_from_R);					
	// This is the matrix in C language. A vector of length columnNum * rowNum.


	// Allocate memory for returning value to R which is kept in a matrix..	
   	PROTECT(returnVal = allocVector(REALSXP, rowNum));	
	for (i=0; i< rowNum; i++) {							// For each row i,

		max = INDEX(matrix, rowNum, i,0);				
		// The fist column is the maximum, 
		// unless there is someting greater than it in other columns.

		for (j=1; j< columnNum; j++) {			// Looking at all other columns,
			if (INDEX(matrix, rowNum, i,j) > max)
				max =INDEX(matrix, rowNum, i,j);
		}
		REAL(returnVal)[i] = max;
	}


	UNPROTECT(2);

	return returnVal;	
}






