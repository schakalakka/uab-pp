#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//Pi const taken from math.h
#define M_PI 3.14159265358979323846

// define rows and columns
#define rows 4096
#define columns 4096

// declare global variables: A and Anew
float A[rows][columns];
float Anew[rows][columns];

int main(int argc, char** argv)
{
	// declare local variables: error, iter_max ...
	float tol = 0.0001f;
	float error = 1;
	int iter_max = 50;
	int iter=0;

	// get iter_max from command line at execution time
	if (argc == 2) iter_max = atoi(argv[1]);



	// set all values in matrix as zero 
	for (int i=1;i<columns-1;i++){
		for (int j=1;j<rows-1;j++){
			A[j][i] = 0;
		}
	}
			
	// set boundary conditions
	//Initialization
	for (int i=0;i<columns;i++){
		A[0][i] = 0;
		A[rows-1][i] = 0;
	}
	for (int i=0;i<rows;i++){
		A[i][0] = sin(i*M_PI/(rows-1));
		A[i][columns-1] = sin(i*M_PI/(rows-1))*exp(-1*M_PI);
	}

    // Main loop: iterate until error <= tol or a maximum of iter_max iterations
	while ( error > tol && iter < iter_max ) {
		// Compute new values using main matrix and writing into auxiliary matrix
		for (int i=1;i<columns-1;i++){
			for (int j=1;j<rows-1;j++){
				Anew[j][i] = (A[j-1][i]+A[j+1][i]+A[j][i-1]+A[j][i+1])/4;
			}
		}	
	
		// Compute error = maximum of the square root of the absolute differences
		error = 0;
		for (int i=1;i<columns-1;i++){
				for (int j=1;j<rows-1;j++){
					if (error < sqrt(fabs(A[j][i]-Anew[j][i])))
						error = sqrt(fabs(A[j][i]-Anew[j][i]));
				}
		}
		// Copy from auxiliary matrix to main matrix
		for (int i=1;i<columns-1;i++){
			for (int j=1;j<rows-1;j++){
				A[j][i] = Anew[j][i];
			}
		}	
	// if number of iterations is multiple of 10 then print error on the screen
        iter++;
	    if (iter % 10 == 0) printf("Iteration: %d, Error: %f\n", error, iter);
	}
    printf("Number of iterations: %d\n, last error: %f", error, iter);
}
