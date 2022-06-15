#ifndef MATRIX_H
#define MATRIX_H

#include <stdbool.h>

typedef struct mat_s
{
	unsigned int num_rows;
	unsigned int num_cols;
	double **data;
	int isSquare;
} mat;

typedef struct mat_lup
{
	mat *L;
	mat *U;
	mat *P;
	unsigned int num_permutations;
} mat_lup;

#define MIN_COEFF 0.00001

// constructing a matrix frame(allocating memory -> all elements are 0)
mat *mat_new(unsigned int num_rows, unsigned int num_cols);

// destructor (destroys memory -> each malloc should have its own free)
void mat_free(mat *matrix);

// constructing the LUP frame(allocating memory and all the elements are 0)
mat_lup *mat_lup_new(mat *L, mat *U, mat *P, unsigned int num_permutations);

// destructor (destroys memory->each malloc should have its own free(RULE))
void mat_lup_free(mat_lup *LUP);

// generate random number in interval
double rand_interval(int min, int max);

// method to create a random matrix
mat *mat_rnd(unsigned int num_rows, unsigned int num_cols, int min, int max);

// method to create a square matrix
mat *mat_sqr(unsigned int rowCol);

// method to create an identity matrix
mat *mat_iden(unsigned int rowCol);

// method to copy matrix to another matrix
mat *mat_cp(mat *matrix);

// read matrix from file[FIX FOR WINDOWS]
mat *mat_fromfile(FILE *f);
// FILE* pointer to the file which we are using

// method to create a random square matrix
mat *mat_sqr_rnd(unsigned int rowCol, int min, int max);

// check if matrix has same dimensions
int mat_equaldim(mat *m1, mat *m2);

// check matrix equality
int mat_equal(mat *m1, mat *m2);

// retreiving the column from the matrixx
mat *mat_get_column(mat *matrix, unsigned int col);

// retreiving the row from  the matrixx
mat *mat_get_row(mat *matrix, unsigned int row);

// setting all elements in matrix to a single number
mat *set_all_elements(mat *matrix, double num);

// setting all main diagonal elements to a single number
mat *set_diagonal_elements(mat *matrix, double num);

// multiply a row with a scalar
mat *row_multipy_scalar(mat *matrix, unsigned int row, double scalar);

// multiply a col with a scalar
mat *col_multiply_scalar(mat *matrix, unsigned int col, double scalar);

// adding rows in a matrix(for gaussian elimination or something like that
mat *rows_add(mat *matrix, unsigned int addendum, unsigned int original, double multiplier);

// printing a matrix
void mat_print(mat *matrix);

// multipy matrix with a scalar
mat *mat_multiply_scalar(mat *matrix, int scalar);

// removing a column in a matrix
mat *mat_remove_column(mat *matrix, unsigned int column);

// removing a row in a matrix
mat *mat_remove_row(mat *matrix, unsigned int row);

// swapping rows
mat *mat_swap_row(mat *matrix, unsigned int row1, unsigned int row2);

// swapping columns
mat *mat_swap_column(mat *matrix, unsigned int col1, unsigned int col2);

// horizontal concatenation of two matrices
mat *mat_horizontal_concat(int mnum, mat **marr);

// vertical concatenation
mat *mat_vertical_concat(int mnum, mat **marr);

// add matrices
mat *mat_add(mat *matrix1, mat *matrix2);

// subtracting matrices
mat *mat_subtract(mat *matrix1, mat *matrix2);

// multiplying matrices(using O(n^3))
mat *mat_mul_naive(mat *matrix1, mat *matrix2);

// multiplyting matrices(using O(n^2.87))
mat *mat_mul_strassen(mat *matrix1, mat *matrix2);

// get the pivot in a column
int mat_get_col_pivot(mat *matrix, unsigned int row1, unsigned int row2);

// row echelon form
mat *mat_ref(mat *matrix);

// reduced row echelon form
mat *mat_rref(mat *matrix);

// LU(P) decomposition
mat_lup *mat_LU(mat *matrix);

// Check if matrix is lower triangular
bool mat_isLowerTriangular(mat *L);

// Check if matrix is upper triangular
bool mat_isUpperTriangular(mat *U);

// forward substituition linear system
mat *solve_linear_forward(mat *L, mat *b);

// backward substitution linear system
mat *solve_linear_backward(mat *U, mat *b);

// linear system using LU(P)
mat *solve_linear_LU(mat *A, mat *b);

// inverse of matrix using LU(P)
mat *mat_inverse();

// determinant of matrix using LU(P)
double *mat_determinant();

// QR decomposition
mat *mat_QR();

#endif
