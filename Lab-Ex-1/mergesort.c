#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void merge(float *X1, float *X2, float *Y, int N1, int N2) {
    int i1 = 0;
    int i2 = 0;
    while (i1 < N1 && i2 < N2) { //copy the elements of the arrays as long no end of one array is reached
        if (X1[i1] <= X2[i2]) {
            Y[i1 + i2] = X1[i1];
            i1++;
        } else {
            Y[i1 + i2] = X2[i2];
            i2++;
        }
    }
    while (i1 < N1) {
        Y[i1 + i2] = X1[i1];
        i1++;
    }
    while (i2 < N2) {
        Y[i1 + i2] = X2[i2];
        i2++;
    }
}

void mergeSort(float *X, float *Y, int N) {
    if (N == 1) {
        Y[0] = X[0];
    } else {
        float T[N];
        int N_half = N / 2; //discarding the decimal points
        mergeSort(X, T, N_half);
        mergeSort(&X[N_half], &T[N_half], N - N_half);
        merge(T, &T[N_half], Y, N_half, N - N_half);
    }
}



int main() {

    int N;
    printf("Length of the dynamically allocated arrays: ");
    scanf("%i", &N);

    //allocate to arrays
    float rand_arr[N];
    float result_arr[N];

    //initialize random with a seed depending on the time
    srand((unsigned int) time(NULL));

    //fill array with random floats (in [0,1]) and print it to the terminal
    printf("Initial random array: \n");
    for (int i = 0; i < N; ++i) {
        rand_arr[i] = (float) rand() / (float) RAND_MAX;
        printf("%f ", rand_arr[i]);
    }

    //output sorted array
    printf("\nSorted array: \n");
    mergeSort(rand_arr, result_arr, N);
    for (int i = 0; i < N; ++i) {
        printf("%f ", result_arr[i]);
    }
}
