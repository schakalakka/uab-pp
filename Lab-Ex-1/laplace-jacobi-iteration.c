// Usage: ./laplace-jacobi-iteration iterations, err_tol, rows, cols


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

#define M_PI 3.14159265358979323846


void print_matrix(float **A, int rows, int columns) {

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            printf("%f, ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

int main(int argc, char **argv) {


    // declare local variables: error, iter_max ...
    int iter = 0;
    float error = 1;

    //default values, can be changed via command line
    float tol = 0.01f;
    int iter_max = 50;
    int rows = 4096;
    int columns = 4096;


    // get iter_max from command line at execution time
    // get tolerance from command line at execution time
    // get rows and columns from command line at execution time
    // Usage: ./laplace-jacobi-iteration iterations, err_tol, rows, cols
    if (argc > 1) iter_max = atoi(argv[1]);
    if (argc > 2) {
        tol = atof(argv[2]);
    }
    if (argc == 4) {
        rows = atoi(argv[3]);
        columns = atoi(argv[3]);
    }
    if (argc == 5) {
        rows = atoi(argv[3]);
        columns = atoi(argv[4]);
    }

    printf("%ix%i matrix\nError tolerance: %f\nMaximum number of iterations: %i\n", rows, columns, tol, iter_max);


    float **A = (float **) calloc(rows, sizeof(float *));
    float **Anew = (float **) calloc(rows, sizeof(float *));

    //initialize each row/line with zeros
    for (int k = 0; k < rows; ++k) {
        A[k] = (float *) calloc(columns, sizeof(float));
        Anew[k] = (float *) calloc(columns, sizeof(float));
    }


    //initialize boundary conditions
    for (int i = 0; i < rows; i++) {
        A[i][0] = sin(i * M_PI / (rows - 1));
        A[i][columns - 1] = sin(i * M_PI / (rows - 1)) * exp(-1 * M_PI);
        Anew[i][0] = sin(i * M_PI / (rows - 1));
        Anew[i][columns - 1] = sin(i * M_PI / (rows - 1)) * exp(-1 * M_PI);
    }


    // Main loop: iterate until error <= tol or a maximum of iter_max iterations
    while (error > tol && iter < iter_max) {

        // Compute new values using main matrix and writing into auxiliary matrix
        for (int j = 1; j < rows - 1; j++) {
            for (int i = 1; i < columns - 1; i++) {
                Anew[j][i] = (A[j - 1][i] + A[j + 1][i] + A[j][i - 1] + A[j][i + 1]) * 0.25f;
            }
        }

        // Compute error = maximum of the square root of the absolute differences
        error = 0.0f;
        for (int j = 1; j < rows - 1; j++) {
            for (int i = 1; i < columns - 1; i++) {
                error = fmaxf(error, fabsf(A[j][i] - Anew[j][i]));
            }
        }
        error = sqrtf(error);

        float *temp_pointer = A;
        A = Anew;
        Anew = temp_pointer;

        iter++;
        // if number of iterations is multiple of 10 then print error on the screen
        if (iter % 10 == 0) {
            printf("Iteration: %i, Error: %f\n", iter, error);
        }

    }
    printf("The total number of iterations: %d\nLast error: %f\n", iter, error);
}
