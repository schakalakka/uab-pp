//
// Created by andy on 03.10.17.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


#define L 0.345678
# define M_PIl		3.141592653589793238462643383279502884L /* pi for double, taken from math.h*/




// define constants X, T
// declare matrix U
int main(int argc, char ** argv) {
    // declare local variables
    unsigned int X;
    unsigned int T;

    // get X from command line at execution time
    // get T from command line at execution time
    // Usage: ./finite-differences X T
    if (argc == 3) {
        X = atoi(argv[1]);
        T = atof(argv[2]);
    }
    else{
        printf("Error: Please provide values for X and T\nUsage: ./finite-differences X T");
        exit(1);
    }

    float **U = (float **) calloc(X+1, sizeof(float *));

    //initialize each row/line with zeros
    for (int k = 0; k < X+1; ++k) {
        U[k] = (float *) calloc(T+2, sizeof(float));
    }

    // initialize positions of matrix U
//    double U[X+1][T+2];
//    for (int i = 0; i <= T+1; ++i) {
//        U[0][i] = 0;
//        U[X][i] = 0;
//    }
    printf("%d, %d", X, T);

    double temp_cos_val = cos(M_PIl/T);
    for (int j = 1; j < X; ++j) {
        U[j][0] = sin(j*M_PIl/X);
        U[j][1] = U[j][0]*temp_cos_val;
    }

    // simulation program: body
    for (int t = 1; t < T+1; ++t) {
        for (int x = 1; x < X; ++x) {
            U[x][t+1] = 2*(1-L)*U[x][t] + L*U[x+1][t] + L*U[x-1][t] - U[x][t-1];
        }
    }
    // obtain checksum of final state and print on screen
    double sum = 0;
    for (int k = 0; k <= X; ++k) {
        sum += U[k][T+1];
    }
    printf("Checksum: %f", sum);
}
