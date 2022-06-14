#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include "matrix.h"
#include <string.h>

#define uint unsigned int
#define MIN_COEFF 0.00001

mat *mat_new(uint num_rows, uint num_cols)
{
	/* STUDY MULTIDIMENSIONAL ARRAYS AND POINTERS*/
	// Check if rows are 0
	if (num_rows == 0)
	{
		printf("Invalid rows");
		return NULL;
	}

	// Check if columns are 0
	if (num_cols == 0)
	{
		printf("InvaliGd columns");
		return NULL;
	}

	// giving mat its memory
	mat *m = (mat *)malloc(sizeof(mat));
	m->num_rows = num_rows;
	m->num_cols = num_cols;
	m->isSquare = (num_rows == num_cols) ? 1 : 0;

	// data is pointer to pointer to double(set of arrays which are a set of doubles)
	// setting the outer layer(row layer)

	m->data = calloc(m->num_rows, sizeof(*m->data));
	// NP_CHECK(m->data);
	uint i;

	// structure established complete structure by adding the columns/values

	for (i = 0; i < m->num_rows; ++i)
	{

		// double de-referencing for individual numbers iterating over rows

		m->data[i] = calloc(m->num_cols, sizeof(double));

		// NP_CHECK(m->data[i]);
	}
	return m;
}
void mat_free(mat *matrix)
{
	uint i;
	// each rows memory is freed
	for (i = 0; i < matrix->num_rows; ++i)
	{
		free(matrix->data[i]);
	}
	// freeing the outer structure
	free(matrix->data);
	// freeing the whole matrix
	free(matrix);
}

double rand_interval(int min, int max)
{
	double d = (double)((rand() % (max - min + 1)) + min);
	return d;
}

mat *mat_rnd(uint num_rows, uint num_cols, int min, int max)
{
	mat *res = mat_new(num_rows, num_cols);
	uint i, j;

	for (i = 0; i < num_rows; ++i)
	{
		for (j = 0; j < num_cols; ++j)
		{
			res->data[i][j] = rand_interval(min, max);
		}
	}
	return res;
}

mat *mat_sqr(uint rowCol)
{
	return mat_new(rowCol, rowCol);
}

mat *mat_sqr_rnd(uint rowCol, int min, int max)
{
	return mat_rnd(rowCol, rowCol, min, max);
}

mat *mat_iden(uint rowCol)
{
	mat *iden = mat_sqr(rowCol);
	int i;
	for (i = 0; i < rowCol; ++i)
	{
		iden->data[i][i] = 1.0;
	}
	return iden;
}

mat *mat_fromfile(FILE *f)
{
	int i, j;
	uint num_rows = 0;
	uint num_cols = 0;
	fscanf(f, "%u", &num_rows);
	fscanf(f, "%u", &num_cols);
	mat *res = mat_new(num_rows, num_cols);
	for (i = 0; i < res->num_rows; ++i)
	{
		for (j = 0; j < res->num_cols; ++j)
		{
			fscanf(f, "%lf\t", &res->data[i][j]);
		}
	}
	return res;
}
int mat_equaldim(mat *m1, mat *m2)
{
	if (m1->num_cols == m2->num_cols && m1->num_rows == m2->num_rows)
	{
		return 1;
	}
	else
	{
		return 0;
	}
}

int mat_equal(mat *m1, mat *m2)
{
	if (mat_equaldim(m1, m2) == 0)
	{
		printf("Dimensions of matrix unequal");
		return 0;
	}
	else
	{
		int i, j;
		int flag = 1;
		for (i = 0; i < m1->num_rows; ++i)
		{
			for (j = 0; j < m1->num_rows; ++j)
			{
				if (m1->data[i][j] == m2->data[i][j])
				{
					continue;
				}
				else
				{
					flag = 0;
					break;
				}
			}
		}
		return flag;
	}
}

mat *mat_get_column(mat *matrix, uint col)
{
	int i;
	mat *column = mat_new(matrix->num_rows, 1);
	for (i = 0; i < matrix->num_rows; ++i)
	{
		column->data[i][0] = matrix->data[i][col - 1];
	}
	return column;
}

mat *mat_get_row(mat *matrix, uint row)
{
	mat *rowMat = mat_new(1, matrix->num_cols);
	// row is contigous set of memory
	// syntax of memcpy is void* mecmcpy(void* to,const void* from,size_t )
	memcpy(rowMat->data[0], matrix->data[row - 1], matrix->num_cols * sizeof(double));
	return rowMat;
}

mat *set_all_elements(mat *matrix, double num)
{
	int i, j;
	mat *copy = mat_cp(matrix);
	for (i = 0; i < copy->num_rows; ++i)
	{
		for (j = 0; j < copy->num_cols; ++j)
		{
			copy->data[i][j] = num;
		}
	}
	return copy;
}

mat *set_diagonal_elements(mat *matrix, double num)
{
	int i;
	mat *copy = mat_cp(matrix);
	for (i = 0; i < copy->num_rows; ++i)
	{
		copy->data[i][i] = num;
	}
	return copy;
}

mat *row_multipy_scalar(mat *matrix, uint row, double scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[row][i] *= scalar;
	}
	return matrix;
}
mat *col_multiply_scalar(mat *matrix, uint col, double scalar)
{
	int i;
	for (i = 0; i < matrix->num_cols; ++i)
	{
		matrix->data[i][col] *= scalar;
	}
	return matrix;
}

mat *rows_add(mat *matrix, uint addendum, uint original, double multiplier)
{
	int i;
	mat *copy = mat_cp(matrix);
	for (i = 0; i < copy->num_cols; ++i)
	{
		copy->data[original][i] = copy->data[original][i] + multiplier * copy->data[addendum][i];
	}
	return copy;
}

void mat_print(mat *matrix)
{
	int i, j;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			printf(" %lf ", matrix->data[i][j]);
		}
		printf("\n");
	}
	printf("\n\n");
}
mat *mat_multiply_scalar(mat *matrix, int scalar)
{
	mat *copy = mat_cp(matrix);
	int i, j;
	for (i = 0; i < copy->num_rows; ++i)
	{
		for (j = 0; j < copy->num_cols; ++j)
		{
			copy->data[i][j] = copy->data[i][j] * scalar;
		}
	}
	return copy;
}

mat *mat_remove_column(mat *matrix, uint column)
{
	if (column > matrix->num_cols)
	{
		perror("Invalid columns");
		return NULL;
	}
	mat *r = mat_new(matrix->num_rows, matrix->num_cols - 1);
	int i, j, k;
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0, k = 0; j < matrix->num_cols; ++j)
		{
			if (column != j)
			{
				// k resets back to 0 after i++ [one j is missing and that's going to lead to an incomplete array]
				// but k has to be contiguous
				r->data[i][k++] = matrix->data[i][j];
			}
		}
	}
	return r;
}

mat *mat_remove_row(mat *matrix, uint row)
{
	if (row > matrix->num_rows)
	{
		perror("Invalid rows");
		return NULL;
	}
	mat *r = mat_new(matrix->num_rows - 1, matrix->num_cols);
	int i, j, k;
	for (i = 0, k = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			if (row != i)
			{
				r->data[k][j] = matrix->data[i][j];
			}
		}
		k++;
	}
	return r;
}

mat *mat_swap_row(mat *matrix, uint row1, uint row2)
{
	if (row1 > matrix->num_rows || row2 > matrix->num_rows)
	{
		perror("Invalid row");
	}
	mat *copy = mat_cp(matrix);
	// swap the rows(contiguous memory locations)
	double *temp = copy->data[row1];
	copy->data[row1] = copy->data[row2];
	copy->data[row2] = temp;
	return copy;
}

mat *mat_swap_column(mat *matrix, uint col1, uint col2)
{
	if (col1 > matrix->num_cols || col2 > matrix->num_cols)
	{
		perror("Invalid column");
	}

	mat *copy = mat_cp(matrix);
	// swap the numbers between the rows
	int i;
	for (i = 0; i < copy->num_rows; ++i)
	{
		double temp = copy->data[i][col1];
		copy->data[i][col1] = copy->data[i][col2];
		copy->data[i][col2] = temp;
	}

	return copy;
}

mat *mat_cp(mat *matrix)
{
	int i, j;
	mat *copy = mat_new(matrix->num_rows, matrix->num_cols);
	for (i = 0; i < matrix->num_rows; ++i)
	{
		for (j = 0; j < matrix->num_cols; ++j)
		{
			copy->data[i][j] = matrix->data[i][j];
		}
	}
	return copy;
}

mat *mat_add(mat *matrix1, mat *matrix2)
{
	if (matrix1->num_rows != matrix2->num_rows || matrix1->num_cols != matrix2->num_cols)
	{
		perror("Cannot add matrices of different dimensions");
		return NULL;
	}
	mat *sum = mat_new(matrix1->num_rows, matrix1->num_cols);
	int i, j;
	for (i = 0; i < matrix1->num_rows; ++i)
	{
		for (j = 0; j < matrix1->num_cols; ++j)
		{
			sum->data[i][j] = matrix1->data[i][j] + matrix2->data[i][j];
		}
	}
	return sum;
}

mat *mat_subtract(mat *matrix1, mat *matrix2)
{
	if (matrix1->num_rows != matrix2->num_rows || matrix1->num_cols != matrix2->num_cols)
	{
		perror("Cannot subtract matrices of different dimensions");
		return NULL;
	}
	mat *diff = mat_new(matrix1->num_rows, matrix1->num_cols);
	int i, j;
	for (i = 0; i < matrix1->num_rows; ++i)
	{
		for (j = 0; j < matrix1->num_cols; ++j)
		{
			diff->data[i][j] = matrix1->data[i][j] - matrix1->data[i][j];
		}
	}
	return diff;
}

mat *mat_mul_naive(mat *matrix1, mat *matrix2)
{
	if (matrix1->num_cols != matrix2->num_rows)
	{
		perror("Cannot multiply matrices, incompatible dimensions");
		return NULL;
	}
	mat *mul = mat_new(matrix1->num_rows, matrix2->num_cols);
	uint i, j, k;
	for (i = 0; i < matrix1->num_rows; ++i)
	{
		for (j = 0; j < matrix2->num_cols; ++j)
		{
			for (k = 0; k < matrix1->num_cols; ++k)
			{
				mul->data[i][j] += matrix1->data[i][k] * matrix2->data[k][j];
			}
		}
	}
	return mul;
}

mat *mat_mul_strassen(mat *matrix1, mat *matrix2)
{
	if (matrix1->num_cols != matrix2->num_rows)
	{
		perror("Cannot multiply matrices, incompatible dimensions");
		return NULL;
	}
	return mat_iden(3);
	// Filled by strassen algorithm
}

int mat_get_col_pivot(mat *matrix, uint col, uint row)
{
	int i;
	for (i = row; i < matrix->num_rows; ++i)
	{
		if (fabs(matrix->data[i][col]) > MIN_COEFF)
		{
			return i;
		}
	}
	return -1;
}

mat *mat_ref(mat *matrix)
{

	// Row echelon form(staircase form)
	//  -- all zeroes at bottom
	//  -- first nonzero entry from left is a 1(leading 1)
	//  -- each leading 1 is to the right of the leading 1s in rows above

	// Gaussian algorithm - convert matrices to REF()
	//  -- If matrix all zeroes in REF
	//  -- For each column
	//  -- Find first row from top containing nonzero entry , move to top
	//  -- Multiplly new top row by 1/a to get leading 1
	//  -- subtract multiples of that row from rows below it to get zeros under leading 1
	mat *copy = mat_cp(matrix);
	int i, j;
	for (i = 0; i < copy->num_cols; ++i)
	{
		int pivot = mat_get_col_pivot(copy, i, i);
		if (pivot < 0)
			continue;
		if (pivot != i)
			copy = mat_swap_row(copy, i, pivot);
		for (j = i + 1; j < copy->num_rows; ++j)
		{
			copy = rows_add(copy, i, j, -(copy->data[j][i]) / copy->data[i][i]);
		}
	}
	return copy;
	// !REFACTOR change into a for loop

	//	num_rows,num_cols = arr.shape
	// if(num_rows != num_cols): return "error"
	// for k in range(num_rows-1):
	//		for i in range(k+1,num_rows):
	//			if arr[i,k] == 0: continue
	//			factor = arr[k,k]/arr[i,k]
	//			for j in range(k,num_rows):
	//				arr[i,j] = arr[k,j] - arr[i,j]*factor
}

mat *mat_rref(mat *matrix)
{
	mat *copy = mat_cp(matrix);
	int i, j;
	for (i = 0; i < copy->num_cols; ++i)
	{
		int pivot = mat_get_col_pivot(copy, i, i);
		if (pivot < 0)
			continue;
		if (pivot != i)
			copy = mat_swap_row(copy, i, pivot);
		copy = row_multipy_scalar(copy, i, 1 / copy->data[i][i]);
		for (j = 0; j < copy->num_rows; ++j)
		{
			if (j != i)
			{
				copy = rows_add(copy, i, j, -(copy->data[j][i]));
			}
			continue;
		}
	}
	return copy;
}

mat_lup *mat_lup_new(mat *L, mat *U, mat *P, uint num_permutations)
{
	mat_lup *r = malloc(sizeof(*r));
	r->L = L;
	r->U = U;
	r->P = P;
	r->num_permutations = num_permutations;
	return r;
}

void mat_lup_free(mat_lup *LUP)
{
	mat_free(LUP->L);
	mat_free(LUP->U);
	mat_free(LUP->P);
	free(LUP);
}

mat_lup *mat_LU(mat *matrix)
{
	// Any square matrix can be written as LU
	// After decomposing A, easy to solve Ax = b
	// LUx = b, where Ux = y
	// Ly = b, which is solved using forward substitution
	// solving x via back substitution(Ux = y)
	if (!matrix->isSquare)
	{
		perror("Invalid dimensions: Please enter a square matrix");
	}
	mat *U = mat_cp(matrix);
	uint i, j;
	uint num_permutations = 0;
	mat *P = mat_iden(U->num_rows);
	mat *L = mat_sqr(U->num_rows);
	for (i = 0; i < U->num_cols; ++i)
	{
		int pivot = mat_get_col_pivot(U, i, i);
		if (pivot < 0)
			continue;
		if (pivot != i)
		{
			U = mat_swap_row(U, i, pivot);
			P = mat_swap_row(P, i, pivot);
			L = mat_swap_row(L, i, pivot);
			num_permutations++;
		}
		for (j = i + 1; j < U->num_rows; ++j)
		{
			L->data[j][i] = (U->data[j][i] / U->data[i][i]);
			U = rows_add(U, i, j, -(U->data[j][i]) / U->data[i][i]);
		}
	}
	L = set_diagonal_elements(L, 1);

	mat_lup *res = mat_lup_new(L, U, P, num_permutations);

	return res;
}

bool mat_isLowerTriangular(mat *L)
{
	uint i, j;
	bool isLowerTriangular = 1;
	for (i = 0; i < L->num_rows; ++i)
	{
		for (j = i + 1; j < L->num_cols; ++j)
		{
			if (L->data[i][j] != 0)
			{
				isLowerTriangular = 0;
				break;
			}
		}
	}
	return isLowerTriangular;
}

bool mat_isUpperTriangular(mat *U)
{
	uint i, j;
	bool isUpperTriangular = 1;
	for (i = 0; i < U->num_rows; ++i)
	{
		for (j = 0; j < i; ++j)
		{
			if (U->data[i][j] != 0)
			{
				isUpperTriangular = 0;
				break;
			}
		}
	}
	return isUpperTriangular;
}

mat *solve_linear_forward(mat *L, mat *b)
{
	if (!mat_isLowerTriangular(L) || b->num_cols != 1)
	{
		perror("Error not lower triangular matrix or b is not vector");
		return mat_new(1, 1);
	}
	if (L->num_rows != b->num_rows)
	{
		perror("Error: Matrix and vector is of incompatible dimensions");
		return mat_new(1, 1);
	}
	uint i, j;
	mat *a = mat_new(L->num_rows, 1);
	for (i = 0; i < L->num_rows; ++i)
	{
		double accum = b->data[i][0];
		for (j = 0; j < i; ++j)
		{
			accum -= L->data[i][j] * a->data[j][0];
		}
		a->data[i][0] = accum / L->data[i][i];
	}
	return a;
}

mat *solve_linear_backward(mat *U, mat *b)
{
	if (!mat_isUpperTriangular(U) || b->num_cols != 1)
	{
		perror("Error not upper triangular matrix or b is not vector");
		return mat_new(1, 1);
	}
	if (U->num_rows != b->num_rows)
	{
		perror("Error: Matrix and vector is of incompatible dimensions");
		return mat_new(1, 1);
	}
	int i, j;
	mat *a = mat_new(U->num_rows, 1);
	for (i = U->num_rows - 1; i >= 0; --i)
	{
		double accum = b->data[i][0];
		for (j = i + 1; j < U->num_cols; ++j)
		{
			accum -= U->data[i][j] * a->data[j][0];
		}
		a->data[i][0] = accum / U->data[i][i];
	}
	return a;
}

mat *solve_linear_LU(mat *A, mat *b)
{
	mat_lup *LU = mat_LU(A);
	mat *Pb = mat_mul_naive(LU->P, b);
	mat *y = solve_linear_forward(LU->L, Pb);
	mat *x = solve_linear_backward(LU->U, y);
    mat_lup_free(LU);

	return x;
}
