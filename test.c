#include <stdio.h>
#include "matrix.h"

int main(int argc, char *argv[])
{
    mat *A = mat_new(3, 3);
    mat *B = mat_new(3, 1);
    A->data[0][0] = 2;
    A->data[0][1] = 5;
    A->data[0][2] = 1;

    A->data[1][0] = 4;
    A->data[1][1] = 5;
    A->data[1][2] = 1;

    A->data[2][0] = 6;
    A->data[2][1] = 5;
    A->data[2][2] = 0;

    B->data[0][0] = 20;
    B->data[1][0] = 10;
    B->data[2][0] = 0;

    mat *x = solve_linear_LU(A, B);
//    mat *x = mat_mul_naive(A,B);
    mat_print(x);
    mat_free(x);
    mat_free(A);
    mat_free(B);
    return 0;
}
