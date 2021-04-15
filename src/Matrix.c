/*
 ============================================================================
 Name        : Matrix.c
 Author      : Jaime Feller
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Function definitions of Matrix.h
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <Math.h>

#include "Matrix.h"

//Internal
static int Matrix_Internal_checkDimensionsExact(Matrix_t* In_1, Matrix_t* In_2){ /*checks if In_1 and In_2 have the exact same dimensions*/
	if(In_1 == NULLMTR || In_1 == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->row != In_2->row || In_1->col != In_2->col){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkIndex(Matrix_t* In, unsigned int row, unsigned int col){/* checks if row and col are valid indexes for matrix M */
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(row >= In->row || col >= In->col || (int)row < 0 || (int)col < 0){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkSquare(Matrix_t* In){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row != In->col){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkTriangularUpper(Matrix_t* In){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row == 1){
		return MATRIX_OK;
	}

	char flag = 0x00;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = 0; j < i; ++j){
			if(In->numbers[i*In->col + j] != 0){
				flag = 0x01;
				break;
			}
		}
	}
	return (flag == 0x00) ? MATRIX_OK : MATRIX_ERROR;
}

static int Matrix_Internal_checkAddress(Matrix_t* In_1, Matrix_t* In_2){

	if(In_1 == In_2){
		return MATRIX_SAME;
	}
	if(In_1 == NULLMTR || In_2 == NULLMTR){
			return MATRIX_ERROR | MATRIX_NULL;
	}
	return MATRIX_DIFFERENT;
}

static int Matrix_Internal_checkTriangularLower(Matrix_t* In){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row == 1){
		return MATRIX_OK;
	}
	char flag = 0x00;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = i + 1; j < In->col; ++j){
			if(In->numbers[i*In->col + j] != 0){
				flag = 0x01;
				break;
			}
		}
	}
	return (flag == 0x00) ? MATRIX_OK : MATRIX_ERROR;
}

/*
static int Matrix_Internal_checkDiagonal(Matrix_t* In){
	if(In == NULLMTR){
	return MATRIX_ERROR | MATRIX_NULL;
	}
	return (CheckTriangularUpper(In) == MATRIX_OK && CheckTriangularLower(In) == MATRIX_OK) ? MATRIX_OK : MATRIX_ERROR;
}
*/

/*
static Matrix_t* Matrix_Internal_hConcat(Matrix_t* In_1, Matrix_t* In_2){
	if(In_1 == NULLMTR || In_2 == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	Matrix_t* Out = Matrix_new(In_1->row, In_1->col + In_2->col, 0);
	int k = 0;
	for(int i = 0; i < In_1->row; ++i){
		for(int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
		for(int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return Out;
}
*/

static Matrix_t* Matrix_Internal_vConcat(Matrix_t* In_1, Matrix_t* In_2){
	Matrix_t* Out = Matrix_new(In_1->row + In_2->row, In_1->col, 0);
	int k = 0;
	for(unsigned int i = 0; i < In_1->row; ++i){
		for(unsigned int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
	}
	for(unsigned int i = 0; i < In_2->row; ++i){
		for(unsigned int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return Out;
}

static Matrix_t* Matrix_Internal_submatrix(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right){
	if(lower < upper || right <  left || Matrix_Internal_checkIndex(In, lower, upper) != MATRIX_OK){
		return Matrix_new(0,0,0);
	}

	Matrix_t* Out = Matrix_new(lower - upper + 1, right - left + 1, 0);
	int k = 0;
	for(unsigned int i = upper; i <= lower; ++i){
		for(unsigned int j = left; j <= right; ++j){
			Out->numbers[k++] = In->numbers[i * (In->col) + j];
		}
	}
	return Out;
}

static Matrix_t* Matrix_Internal_removeRow(Matrix_t* In, unsigned int row){
	Matrix_t* Out = Matrix_new(In->row-1, In->col, 0);
	int k = 0;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = 0; j < In->col; ++j){
			if(i != row){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}
	return Out;
}

/*
static Matrix_t* Matrix_Internal_removeColumn(Matrix_t* In, unsigned int col){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	Matrix_t* Out = Matrix_new(In->row-1, In->col, 0);
	int k = 0;
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < In->col; ++j){
			if(j != col){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}
	return Out;
}
*/

/*
//Matrix_Internal_transpose is defined abode Matrix_transpose
static Matrix_t* Matrix_Internal_transpose_2(Matrix_t* In){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	Matrix_t* Out = Matrix_new(In->col, In->row, 0);
	int k = 0;
	for(int i = 0; i < In->col; ++i){
		for(int j = 0; j < In->row; ++j){
			Out->numbers[k++] = In->numbers[j * In->row + i];
		}
	}
	return Out;
}
*/

//Constructors
Matrix_t* Matrix_new(size_t row, size_t col, double val){
	if((int)row < 0 || (int)col < 0){
		return (Matrix_t*)NULL;
	}

	Matrix_t* Out = (Matrix_t *)malloc(sizeof(Matrix_t));
	Out->row = row;
	Out->col = col;
	Out->numbers = (double*)malloc(sizeof(double)*row*col);
	int k = 0;
	for(unsigned int i = 0; i < (Out->row * Out->col); ++i){
		Out->numbers[k++] = val;
	}
	return Out;
}

void Matrix_free(Matrix_t* In){
	if(In == NULLMTR){
		free((Matrix_t*) In);
		return;
	}
	free((double*)In->numbers);
	free((Matrix_t*)In);
}

Matrix_t* Matrix_identity(size_t n){
	if((int)n < 0){
		return NULL;
	}

	Matrix_t* I = Matrix_new(n, n, 0);
	for(unsigned int i = 0; i < n; ++i){
		I->numbers[i*I->col + i] = 1;
	}
	return I;
}

Matrix_t* Matrix_zeroes(size_t row, size_t col){
	if((int)row < 0 || (int)col < 0){
		return NULL;
	}

	Matrix_t* Z = Matrix_new(row, col, 0);
	return Z;
}

Matrix_t* Matrix_ones(size_t row, size_t col){
	if((int)row < 0 || (int)col < 0){
		return NULL;
	}
	Matrix_t* A = Matrix_new(row, col, 1);
	return A;
}

Matrix_t* Matrix_random(size_t row, size_t col, double from, double to){
	if((int)row < 0 || (int)col < 0){
		return NULL;
	}

	Matrix_t* R = Matrix_new(row, col, 0);
	int k = 0;
	for(unsigned int i = 0; i < row * col ; ++i){
		double r = ((double)rand()) / ((double)RAND_MAX);
		R->numbers[k++] = from + (to-from)*r;
	}
	return R;
}

Matrix_t* Matrix_random_int(size_t row, size_t col, int from, int to){
	if((int)row < 0 || (int)col < 0){
		return NULL;
	}

	Matrix_t* R = Matrix_new(row, col, 0);
	int k = 0;
	for(unsigned int i = 0; i < row * col ; ++i){
		double r = ((double)rand()) / ((double)RAND_MAX);
		R->numbers[k++] = (int)(from + (to-from)*r);
	}
	return R;
}

int Matrix_get(Matrix_t* In, unsigned int row, unsigned int col, double *out){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkIndex(In, row, col) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_INDEX;
	}
	*out = In->numbers[row * In->col + col];
	return MATRIX_OK;
}

int Matrix_set(Matrix_t* In, unsigned int row, unsigned int col, double val){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkIndex(In, row, col) != MATRIX_OK){
			return MATRIX_ERROR | MATRIX_INDEX;
		}
	In->numbers[ row * In->col + col] = val;
	return MATRIX_OK;
}

//Algebra
int Matrix_scalarMultiplication(Matrix_t* In, double factor, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkDimensionsExact(In, Out)){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	for(unsigned int k = 0; k < In->row * In->col; ++k){
		Out->numbers[k] = In->numbers[k] * factor;
	}
	return MATRIX_OK;
}

int Matrix_addition(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkDimensionsExact(In_1, In_2) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In_1, Out) != MATRIX_OK){
			return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	for(unsigned int k = 0; k < In_1->row * In_1->col; ++k){
		Out->numbers[k] = In_1->numbers[k] + In_2->numbers[k];
	}
	return MATRIX_OK;
}

int Matrix_subtraction(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkDimensionsExact(In_1, In_2) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In_1, Out) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	for(unsigned int k = 0; k < In_1->row * In_1->col; ++k){
		Out->numbers[k] = In_1->numbers[k] - In_2->numbers[k];
	}
	return MATRIX_OK;
}

static void Matrix_Internal_multiplication(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	// [n x m] * [m x p] = [n x q].

	if(In_1->row == 1 && In_1->col == 1){
		Matrix_scalarMultiplication(In_2, In_1->numbers[0], Out);
		return;
	}
	if(In_2->row == 1 && In_2->col == 1){
		Matrix_scalarMultiplication(In_1, In_2->numbers[0], Out);
		return;
	}

	for(unsigned int i = 0; i < In_1->row; ++i){
		for(unsigned int j = 0; j < In_2->col; ++j){
			for(unsigned int k = 0; k < In_2->row;++k){
				Out->numbers[i * Out->col + j] += In_1->numbers[i * (In_1->col) + k] * In_2->numbers[k * (In_2->col) + j];
			}
		}
	}

	return;
}

int Matrix_multiplication(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->col != In_2->row || In_1->row != Out->row || In_2->col != Out->col){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	//case A = AxB || B = AxB || A=AxA
	if( (Matrix_Internal_checkAddress(In_1, Out) == MATRIX_SAME) || (Matrix_Internal_checkAddress(In_2, Out) == MATRIX_SAME)){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);
		Matrix_Internal_multiplication(In_1, In_2, temp);
		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}

	Matrix_setMatrixtoZero(Out);
	Matrix_Internal_multiplication(In_1, In_2, Out);

	return MATRIX_OK;
}

static void Matrix_Internal_multiplication_recursive(
		Matrix_t* A, unsigned int a_upper, unsigned int a_lower, unsigned int a_left, unsigned int a_right,
		Matrix_t* B, unsigned int b_upper, unsigned int b_lower, unsigned int b_left, unsigned int b_right,
		Matrix_t* C){

	//printf("A: %d %d %d %d\n", a_upper, a_lower, a_left, a_right);
	//printf("B: %d %d %d %d\n", b_upper, b_lower, b_left, b_right);

	//          1          x         1                     1          x         1
	if( (a_upper == a_lower && a_left == a_right) && (b_upper == b_lower && b_left == b_right)){
		C->numbers[a_upper * C->col + b_left] +=
				A->numbers[ a_upper * A->col + a_right] * B->numbers[b_upper * B->col + b_left];
		return;
	}

	//           1          x          n                      n         x          1
	if( (a_upper == a_lower && a_right != a_left) && (b_upper != b_lower && b_right == b_left) ){
		// if A = [A11 A12]  and B = [B11  then C = [A11xB11+A12xB21]
		//	 	   			          B21]

		//C11_a = A11 x B11
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
				B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
					(b_left)													, (b_right),
				C);

		//C11_b = A12 x B21
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
				B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
					(b_left)													, (b_right),
				C);
		return;
	}

	//          1 x n                   n x p
	if( (a_upper == a_lower)  && (b_right != b_left) ){
		if(a_left == a_right){
			// if A = [A11]  and B = [B11 B12]  then C = [A11xB11 A11xB12]
			//

			//C11 = A11 x B11
			Matrix_Internal_multiplication_recursive(
					A,	(a_upper)													, (a_lower),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
					C);

			//C12 = A12 x B21
			Matrix_Internal_multiplication_recursive(
					A,	(a_upper)													, (a_lower),
						(a_left)													, (a_right),
					B,	(b_upper) 													, (b_lower),
						((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
					C);
			return;
		}
		// if A = [A11 A12]  and B = [B11 B12  then C = [A11xB11+A12xB21 A11xB12+A12xB22]
		//	 	   			          B21 B22]

		//C11_a = A11 x B11
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
				B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
					(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
				C);

		//C11_b = A12 x B21
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
				B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
					(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
				C);

		//C12_a = A11 x B12
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
				B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
					((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
				C);

		//C12_b = A12 x B22
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, (a_lower),
					((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
				B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
					((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
				C);
		return;
	}

	//          m x n                   n x 1
	if( (a_upper != a_lower) && (b_right == b_left)){
		if(a_left == a_right){
			// if A = [A11  and B = [B11]  then C = [A11xB11
			//         A21]							 A21xB11]

			//C11 = A11 x B11
			Matrix_Internal_multiplication_recursive(
					A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						(b_left)													, (b_right),
					C);

			//C21 = A21 x B11
			Matrix_Internal_multiplication_recursive(
					A,	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						(b_left)													, (b_right),
					C);

			return;
		}
		// if A = [A11 A12  and B = [B11  then C = [A11xB11+A12xB21
		//	 	   A21 A22]          B21]           A21xB11+A22xB21]

		//C11_a = A11 x B11
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
					(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
				B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
					(b_left)													, (b_right),
				C);

		//C11_b = A12 x B21
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
					((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
				B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
					(b_left)													, (b_right),
				C);

		//C21_a = A21 x B11
		Matrix_Internal_multiplication_recursive(
				A,	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
					(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
				B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
					(b_left)													, (b_right),
				C);

		//C21_b = A22 x B21
		Matrix_Internal_multiplication_recursive(
				A, 	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
					((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
				B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
					(b_left)													, (b_right),
				C);
		return;
	}

	//           n         x         1                       1         x           n
	if((a_upper != a_lower && a_right == a_left) && (b_upper == b_lower && b_right != b_left)){
		// if A = [A11  and B = [B11 B12]  then C = [A11xB11 A11xB12
		//	 	   A21]                     		 A21xB11 A21xB12]

		//C11 = A11 x B11
		Matrix_Internal_multiplication_recursive(
				A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
					(a_left)													, (a_right),
				B,	(b_upper)													, (b_lower),
					(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
				C);
		//C12 = A11 x B12
		Matrix_Internal_multiplication_recursive(
					A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
					C);
		//C21 = A21 x B11
		Matrix_Internal_multiplication_recursive(
					A,	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
					C);
		//C22 = A21 x B12
		Matrix_Internal_multiplication_recursive(
					A, 	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
						(a_left)													, (a_right),
					B,	(b_upper)													, (b_lower),
						((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
					C);
		return;
	}

	//INDEXES
	/*
	//int index_uppe = (a_upper);
	//int index_miup = ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1);
	//int index_milo = ((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1));
	//int index_lowe = (a_lower);

	//int index_left = (a_left);
	//int index_mile = ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1);
	//int index_miri = ((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1));
	//int index_righ = (a_right);
	*/

	//multiply the matrixes

	// if A = [A11 A12  and B = [B11 B12  then C = [A11xB11+A12xB21 A11xB12+A12xB22
	//	 	   A21 A22]          B21 B22]           A21xB11+A22xB21 A21xB12+A22xB22]

	//C11_a = A11 x B11
	Matrix_Internal_multiplication_recursive(
			A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
				(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
			B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
				(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
			C);
	//C11_b = A12 x B21
	Matrix_Internal_multiplication_recursive(
			A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
				((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
			B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
				(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
			C);
	//C12_a = A11 x B12
	Matrix_Internal_multiplication_recursive(
			A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
				(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
			B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
				((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
			C);
	//C12_b = A12 x B22
	Matrix_Internal_multiplication_recursive(
			A,	(a_upper)													, ((a_upper) + (((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)) - 1),
				((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
			B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
				((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
			C);
	//C21_a = A21 x B11
	Matrix_Internal_multiplication_recursive(
			A,	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
				(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
			B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
				(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
			C);
	//C21_b = A22 x B21
	Matrix_Internal_multiplication_recursive(
			A, 	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
				((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
			B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
				(b_left)													, ((b_left) + (((b_right-b_left)>>1) + ((b_right-b_left)&1)) - 1),
			C);
	//C22_a = A21 x B12
	Matrix_Internal_multiplication_recursive(
			A, 	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
				(a_left)													, ((a_left) + (((a_right-a_left)>>1) + ((a_right-a_left)&1)) - 1),
			B,	(b_upper)													, ((b_upper) + (((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)) - 1),
				((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
			C);
	//C22_b = A22 x B22
	Matrix_Internal_multiplication_recursive(
			A, 	((a_upper) + ((a_lower-a_upper)>>1) + ((a_lower-a_upper)&1)), (a_lower),
				((a_left) + ((a_right-a_left)>>1) + ((a_right-a_left)&1))	, (a_right),
			B,	((b_upper) + ((b_lower-b_upper)>>1) + ((b_lower-b_upper)&1)), (b_lower),
				((b_left) + ((b_right-b_left)>>1) + ((b_right-b_left)&1))	, (b_right),
			C);

	return;
}

int Matrix_multiplication_recursive(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	// [n x m] * [m x p] = [n x q].
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->col != In_2->row || In_1->row != Out->row || In_2->col != Out->col){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	//case A = AxB || B = AxB || A=AxA
	if( (Matrix_Internal_checkAddress(In_1, Out) == MATRIX_SAME) || (Matrix_Internal_checkAddress(In_2, Out) == MATRIX_SAME)){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);
		Matrix_Internal_multiplication_recursive(In_1, 0, In_1->row-1, 0, In_1->col-1,
								 	 	 	 	 In_2, 0, In_2->row-1, 0, In_2->col-1,
												 temp);
		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}

	Matrix_setMatrixtoZero(Out);
	Matrix_Internal_multiplication_recursive(In_1, 0, In_1->row-1, 0, In_1->col-1,
							 	 	 	 	 In_2, 0, In_2->row-1, 0, In_2->col-1,
											 Out);
	return MATRIX_OK;
}

//definitions

static void Matrix_Internal_transpose(Matrix_t* In, Matrix_t* Out){
	int k = 0;
	for(unsigned int i = 0; i < In->col; ++i){
		for(unsigned int j = 0; j < In->row; ++j){
			Out->numbers[k++] = In->numbers[j * In->row + i];
		}
	}
}

int Matrix_transpose(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row != Out->col || In->col != Out->row){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	//case square matrix, transpose in itself
	if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

		Matrix_Internal_transpose(In, temp);

		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}
	Matrix_Internal_transpose(In, Out);
	return MATRIX_OK;

}

int Matrix_determinant(Matrix_t* In, double* out){
	if( In == NULLMTR || out == (double*)0 ){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}

	if(In->row == 1){
		*out = In->numbers[0];
		return MATRIX_OK;
	}

	Matrix_t* M1 = Matrix_Internal_removeRow(In, 0);
	Matrix_t* M2 = Matrix_zeroes(In->row-1, In->col-1);
	double d = 0, sign = +1;
	for(unsigned int j = 0; j < In->col; ++j){
		double c = In->numbers[j];
		Matrix_removeColumn(M1, j, M2);
		double temp = 0;
		Matrix_determinant(M2, &temp);
		d += sign * temp * c;
		sign*=-1;
	}
	Matrix_free(M1);
	Matrix_free(M2);

	*out = d;

	return MATRIX_OK;
}

int Matrix_trace(Matrix_t* In, double* out){
	if(In == NULLMTR || out == (double*)0){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}
	*out = 0;
	for(unsigned int i = 0; i < In->col; ++i){
		*out += In->numbers[i*In->row+i];
	}
	return MATRIX_OK;
}

static void Matrix_Internal_adjoint(Matrix_t* In, Matrix_t* Out){
	Matrix_t* A1 = Matrix_new(In->row - 1, In->col, 0);
	Matrix_t* A2 = Matrix_new(In->row - 1, In->col - 1, 0);

	for(unsigned int i = 0; i < In->row;++i){
		Matrix_removeRow(In,i,A1);
		for(unsigned int j = 0; j < In->col;++j){
			Matrix_removeColumn(A1, j, A2);
			double si = (i+j)%2 == 0 ? 1: -1;
			double determ = 0;
			Matrix_determinant(A2, &determ);
			Out->numbers[i*Out->col+j] = (double)si*(double)determ;
		}
	}
	Matrix_transpose(Out, Out);

	Matrix_free(A1);
	Matrix_free(A2);

	return;
}

int Matrix_adjoint(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if( (Matrix_Internal_checkSquare(In) != MATRIX_OK) || Matrix_Internal_checkDimensionsExact(In, Out) != MATRIX_OK || (Matrix_Internal_checkSquare(Out) != MATRIX_OK) ){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}
	if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

		Matrix_Internal_adjoint(In, temp);

		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}
	Matrix_Internal_adjoint(In, Out);
	return MATRIX_OK;
}

static int Matrix_Internal_inverse(Matrix_t* In, Matrix_t* Out){
	Matrix_t* B = Matrix_new(In->row, In->col, 0);

    Matrix_adjoint(In, B);
	double det = 0;
	Matrix_determinant(In, &det);
	if(det == 0){
		Matrix_free(B);
		return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
	}
	Matrix_scalarMultiplication(B, 1/det, Out);
	Matrix_free(B);

	return MATRIX_OK;
}

int Matrix_inverse(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if( (Matrix_Internal_checkSquare(In) != MATRIX_OK) || Matrix_Internal_checkDimensionsExact(In, Out) != MATRIX_OK || (Matrix_Internal_checkSquare(Out) != MATRIX_OK) ){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}
	if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

		if( Matrix_Internal_inverse(In, temp) != MATRIX_OK){
			Matrix_free(temp);
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}

		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}
	if( Matrix_Internal_inverse(In, Out) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_inverse_triangular_upper(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	for(unsigned int i = 0; i < Out->row; ++i){
		//if element in diagonal is 0, then non invertible
		if(In->numbers[i * (In->col) + i] == 0){
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}
		for(unsigned int j = 0; j < Out->col; j++){
			if(j < i){
				Out->numbers[i*Out->col+j] = 0;
			} else if(j == i){
				Out->numbers[i*Out->col+j] = 1/In->numbers[i * (In->col) + j];
			} else {
				Out->numbers[i*Out->col+j] = In->numbers[i*(In->col) + j]/In->numbers[j*In->col+j];
			}
		}
	}
	return MATRIX_OK;
}

static int Matrix_Internal_inverse_triangular_lower(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	for(unsigned int i = 0; i < Out->row; ++i){
		//if element in diagonal is 0, then non invertible
		if(In->numbers[i * (In->col) + i] == 0){
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}
		for(unsigned int j = 0; j < Out->col; j++){
			if(j < i){
				Out->numbers[i*Out->col+j] = In->numbers[i*(In->col) + j]/In->numbers[j*In->col+j];
			} else if(j == i){
				Out->numbers[i*Out->col+j] = 1/In->numbers[i * (In->col) + j];
			} else {
				Out->numbers[i*Out->col+j] = 0;
			}
		}
	}
	return MATRIX_OK;
}

int Matrix_inverse_triangular(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}

	if( (Matrix_Internal_checkSquare(In) != MATRIX_OK) || Matrix_Internal_checkDimensionsExact(In, Out) != MATRIX_OK || (Matrix_Internal_checkSquare(Out) != MATRIX_OK) ){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}

	if(Matrix_Internal_checkTriangularUpper(In) == MATRIX_OK){
		if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
			Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

			if( Matrix_Internal_inverse_triangular_upper(In, temp) != MATRIX_OK){
				Matrix_free(temp);
				return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
			}

			Matrix_copyTo(temp, Out);
			Matrix_free(temp);
			return MATRIX_OK;
		}

		return Matrix_Internal_inverse_triangular_upper(In, Out);
	}

	if(Matrix_Internal_checkTriangularLower(In) == MATRIX_OK){
		if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
			Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

			if( Matrix_Internal_inverse_triangular_lower(In, temp) != MATRIX_OK){
				Matrix_free(temp);
				return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
			}

			Matrix_copyTo(temp, Out);
			Matrix_free(temp);
			return MATRIX_OK;
		}

		return Matrix_Internal_inverse_triangular_lower(In, Out);
	}

	return MATRIX_ERROR | MATRIX_TRIANGULAR;
}

// Utilities
void Matrix_display(Matrix_t* In){
	if(In == NULLMTR){
		printf("NULL MATRIX\n");
		return;
	}
	if(In->row > 0 && In->col > 0){
		int k = 0;
		printf("[");
		for(unsigned int i = 0; i < In->row; ++i){
			for(unsigned int j = 0; j < In->col; ++j){
				if(j < In->col - 1){
					printf("\t%f", In->numbers[k++]);
				} else {
					printf("\t%f", In->numbers[k++]);
				}
			}
			if(i < In->row - 1){
				printf("\n");
			} else {
				printf("\t]\n");
			}
		}
	} else {
		printf("[\t]\n");
	}
	printf("\n");
}

Matrix_t* Matrix_copy(Matrix_t* In){
	if(In == NULLMTR){
		return NULLMTR;
	}
	Matrix_t* Out = Matrix_new(In->row, In->col, 0);
	for(unsigned int k = 0; k < In->row * In->col; k++){
		Out->numbers[k]=In->numbers[k];
	}
	return Out;
}

int Matrix_copyTo(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkDimensionsExact(In, Out)){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	for(unsigned int k = 0; k < In->row * In->col; k++){
			Out->numbers[k]=In->numbers[k];
		}
	return MATRIX_OK;
}

int Matrix_submatrix(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(lower < upper || right <  left || Matrix_Internal_checkIndex(In, lower, right) != MATRIX_OK){
		return MATRIX_ERROR || MATRIX_INDEX;
	}
	if(Out->row !=  lower-upper+1 || Out->col != right-left+1){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}

	int k = 0;
	for(unsigned int i = upper; i <= lower; ++i){
		for(unsigned int j = left; j <= right; ++j){
			Out->numbers[k++] = In->numbers[i*In->col + j ];
		}
	}
	return MATRIX_OK;
}

int Matrix_removeRow(Matrix_t* In, unsigned int row, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row != Out->row + 1 || In->col != Out->col ){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = 0; j < In->col; ++j){
			if(i != row){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}
	return MATRIX_OK;
}

int Matrix_removeColumn(Matrix_t* In, unsigned int col, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row != Out->row || In->col != Out->col + 1 ){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = 0; j < In->col; ++j){
			if(j != col){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}

	return MATRIX_OK;
}

Matrix_t* Matrix_rowEchelon(Matrix_t* In){
	//__TODO__ Matrix_rowEchelon(Matrix_t* In, Matrix_t* Out)
	if(In->row == 1){
		for(unsigned int j = 0; j < In->col;++j){
			//search first element different to 0 and normalize
			if( !(fabs(In->numbers[j]-0.0) <= FLT_EPSILON) ){
				Matrix_t* B = Matrix_new(In->row, In->col, 0);
				Matrix_scalarMultiplication(In, 1/In->numbers[j], B);
				return B;
			}
		}
		Matrix_t* B = Matrix_new(In->row, In->col, 0);
		return B;
	}
	Matrix_t* B = Matrix_copy(In);
	unsigned int ind1 = B->col;
	unsigned int ind2 = 0;
	//search first not zero entry
	for(unsigned int i = 0; i < B->row; ++i){
		for(unsigned int j = 0; j < B->col; ++j){
			if(B->numbers[i*B->col+j]!= 0 && j < ind1){
				ind1 = j;
				ind2 = i;
				break;
			}
		}
	}
	if(ind2 > 1){
		//if I am not in the first row, swap it with the first
		for(unsigned int j = 0; j < B->col; ++j){
			double temp = B->numbers[j];
			B->numbers[j] = B->numbers[ind2*B->col+j];
			B->numbers[ind2*B->col+j] = temp;
		}
	}
	if(B->numbers[0]!=0){
		//normalize first row so the first digit is 1 (or 0)
		double coef = B->numbers[0];
		for(unsigned int j = 0; j < B->col; ++j){
			B->numbers[j] /= coef;
		}
		//Row = row - row_coef * row_0
		for(unsigned int i = 1; i < B->row; ++i){
			coef = B->numbers[i*B->col];
			for(unsigned int j = 0; j < B->col; ++j){
				B->numbers[i*B->col+j] -= coef * B->numbers[j];
			}
		}
	} else {
		//normalize first row
		double coef = 0;
		for(unsigned int j = 0; j < B->col;++j){
			if(B->numbers[j] != 0 && coef == 0){
				coef = B->numbers[j];
				B->numbers[j] = 1;
			} else if(B->numbers[j]!= 0){
				B->numbers[j] /= coef;
			}
		}
		for(unsigned int i = 1; i < B->row; ++i){
			coef = B->numbers[i*B->col];
			for(unsigned int j = 0; j < B->col; ++j){
				if(B->numbers[i*B->col+j] != 0 && coef == 0){
					coef = B->numbers[i*B->col+j];
					B->numbers[i*B->col+j] = 0;
				} else if(B->numbers[j]!= 0){
					B->numbers[i*B->col+j] -= coef * B->numbers[j];
				}
			}
		}
	}
	//remove first row, repeat process
	Matrix_t* B1 = Matrix_Internal_removeRow(B,0);
	Matrix_t* Be = Matrix_rowEchelon(B1);

	//copy results from Be to B
	for(unsigned int i = 0; i < Be->row; ++i){
		for(unsigned int j = 0; j < Be->col; ++j ){
			B->numbers[(i+1)*B->col+j] = Be->numbers[(i)*Be->col+j];
		}
	}
	Matrix_free(B1);
	Matrix_free(Be);
	return B;
}

static Matrix_t* Matrix_Internal_reducedRowEchelon(Matrix_t* In){
	if(In->row == 1){
		Matrix_t* R = Matrix_copy(In);
		return R;
	}
	Matrix_t* R = Matrix_copy(In);
	//find first non-zero digit in the last
	int idx1 = In->row -1, idx2 = In->col - 1;
	char flag = 0x00;

	for(int i = R->row - 1; i >= 0 && flag == 0x00; --i){
		for(unsigned int j = 0; j < R->col; ++j){
			if(R->numbers[i*R->col+j] != 0.0){
				idx1 = i; idx2 = j;
				flag = 0x01;
				break;
			}
		}
	}
	//subtract in the rows abode
	for(int i = idx1-1; i >= 0; --i){
		double coef = R->numbers[i*R->col+idx2];
		for(unsigned int j = 0; j < R->col; ++j){
			R->numbers[i*R->col + j] -= coef*R->numbers[idx1*R->col + j];
		}
	}
	//eliminate the lower n rows and calculate recursively
	Matrix_t* R1 = Matrix_Internal_submatrix(R, 0, idx1-1, 0, R->col-1);
	Matrix_t* R2 = Matrix_Internal_reducedRowEchelon(R1);
	Matrix_free(R1);
	//copy values
	for(int i = 0; i < idx1; ++i){
		for(unsigned int j = 0; j < R->col; ++j){
			R->numbers[i*R->col + j] = R2->numbers[i*R->col + j];
		}
	}
	Matrix_free(R2);

	return R;
}

Matrix_t* Matrix_reducedRowEchelon(Matrix_t* In){
	Matrix_t* R1 = Matrix_rowEchelon(In);
	Matrix_t* R2 = Matrix_Internal_reducedRowEchelon(R1);
	Matrix_free(R1);
	return R2;
}

int Matrix_hConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->col + In_2->col != Out->col || In_1->row != In_2->row || In_1->row != Out->row){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(unsigned int i = 0; i < In_1->row; ++i){
		for(unsigned int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
		for(unsigned int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return MATRIX_OK;
}

int Matrix_vConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->row + In_2->row != Out->row || In_1->col != In_2->col || In_1->col != Out->col){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(unsigned int i = 0; i < In_1->row; ++i){
		for(unsigned int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
	}
	for(unsigned int i = 0; i < In_2->row; ++i){
		for(unsigned int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return MATRIX_OK;
}

double Matrix_norm(Matrix_t* In){
	double d = 0;
	for(unsigned int k = 0; k < In->row*In->col; ++k){
		d += In->numbers[k]*In->numbers[k];
	}
	return sqrt(d);
}

Matrix_t* Matrix_kernel(Matrix_t* In){
	return Matrix_nullSpace(In);
}

Matrix_t* Matrix_nullSpace(Matrix_t* In){
	Matrix_t* RM = Matrix_reducedRowEchelon(In);

	int k = RM->row-1;
	for(int i = RM->row-1; i >= 0; --i){
		char flag = 0x00;
		for(unsigned int j =0; j < RM->col; ++j){
			if(RM->numbers[i*RM->col+j]!=0){
				flag = 0x01;
				break;
			}
		}
		if(flag == 0x01){
			k = i;
			break;
		}
	}

	Matrix_t* RRM = Matrix_Internal_submatrix(RM, 0, k, 0, RM->col-1);
	Matrix_free(RM);

	unsigned int nn = RRM->col - RRM->row;

	if(nn == 0){
		Matrix_t* N = Matrix_new(0,0,0);
		return N;
	}

	//obtain free columns R2
	Matrix_t* F = Matrix_Internal_submatrix(RRM, 0, RRM->row - 1, RRM->row, RRM->col-1);
	Matrix_free(RRM);


	Matrix_t* I = Matrix_identity(nn);

	//reciprocate F
	Matrix_scalarMultiplication(F, -1, F);

	//concat with identity
	Matrix_t* N = Matrix_Internal_vConcat(F,I);
	Matrix_free(F);
	Matrix_free(I);

	return N;
}

int Matrix_descomposition_LU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_U){
	if(In == NULLMTR || Out_L == NULLMTR || Out_U == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	//check squareness
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In, Out_L) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In, Out_U) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}

	if(Matrix_Internal_checkAddress(In, Out_L) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(In, Out_U) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(Out_L, Out_U) == MATRIX_SAME){
		return MATRIX_ERROR | MATRIX_SAME;
	}

	if(In->row == 1){
		Out_U->numbers[0] = In->numbers[0];
		Out_L->numbers[0] = 1;
		return MATRIX_OK;
	}

	double a = In->numbers[0];
	double c = 0;
	if(a != 0){
		c = 1/a;
	}

	Matrix_t* w  = Matrix_Internal_submatrix(In, 0, 0, 1, In->col-1);
	Matrix_t* v  = Matrix_Internal_submatrix(In, 1, In->row-1, 0,0);
	Matrix_t* Ab = Matrix_Internal_submatrix(In, 1, In->row-1, 1, In->col-1);
	Matrix_t* T1 = Matrix_new(v->row, w->col, 0);

	Matrix_multiplication(v, w, T1);
	Matrix_scalarMultiplication(T1, -c, T1);
	Matrix_addition(Ab, T1, T1);
	Matrix_free(Ab);

	Matrix_t* l = Matrix_new(T1->row, T1->col, 0);
	Matrix_t* u = Matrix_new(T1->row, T1->col, 0);
	Matrix_descomposition_LU(T1, l, u);
	Matrix_free(T1);

	//put the descomposition inside our matrixes
	int k = 0;
	for(unsigned int i = 0; i < In->row; ++i){
		for(unsigned int j = 0; j < In->col; ++j){
			if(i == 0 && j == 0){
				Out_L->numbers[k] = 1;
				Out_U->numbers[k] = a;
				++k;
			} else if(i == 0 && j > 0) {
				Out_U->numbers[k] = w->numbers[j-1];
				Out_L->numbers[k] = 0;
				++k;
			} else if(i > 0 && j == 0){
				Out_L->numbers[k] = c*v->numbers[i-1];
				Out_U->numbers[k] = 0;
				++k;
			} else{
				Out_L->numbers[k] = l->numbers[(i-1)*l->col + (j - 1)];
				Out_U->numbers[k] = u->numbers[(i-1)*u->col + (j - 1)];
				++k;
			}
		}
	}

	Matrix_free(w);
	Matrix_free(v);
	Matrix_free(l);
	Matrix_free(u);

	return MATRIX_OK;
}

int Matrix_descomposition_LDU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_D, Matrix_t* Out_U){
	if(In == NULLMTR || Out_L == NULLMTR || Out_D == NULLMTR || Out_U == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_L) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_D) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_U) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}

	if(Matrix_Internal_checkAddress(In, Out_L) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(In, Out_D) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(In, Out_U) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(Out_L, Out_D) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(Out_D, Out_U) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(Out_L, Out_U) == MATRIX_SAME){
		return MATRIX_ERROR | MATRIX_SAME;
	}

	Matrix_descomposition_LU(In, Out_L, Out_U);
	//Matrix U has the coeficients of Matrix D
	//must assure that matrix D has zeros in non diagonal index;
	for(unsigned int i = 0; i < In->col; ++i){
		for(unsigned int j = 0; j < In->row; ++j){
			if(i == j){
				Out_D->numbers[i*In->col+j] = Out_U->numbers[i*In->col+j];
				Out_U->numbers[i*In->col+j] = 1;
			}
		}
	}
	return MATRIX_OK;
}

static void Matrix_Internal_orthogonalization(Matrix_t* In, Matrix_t* Out){
	Matrix_t* In_col = Matrix_new(In->row, 1, 0);
	Matrix_t* Out_col = Matrix_new(In->row, 1, 0);
	Matrix_t* temp = Matrix_new(In->row, 1, 0);
	//u1 = v1
	//e1 = 1 / ||u1|| * u1
	//u_n = v_n - sum( [<V_n, U_k>/<u_k, u_k>]*u_k )
	//e_n = u_n / ||u_n||
	for(unsigned int j = 0; j < In->col; ++j){
		Matrix_submatrix(In, 0, In->row-1, j, j, In_col); // V_n
		Matrix_submatrix(In, 0, In->row-1, j, j, temp);
		double euclideannorm = 0;
		for(unsigned int k = 0; k < j; ++k){
			Matrix_submatrix(Out, 0, Out->row-1, k, k, Out_col); // U_k
			double d1 = 0, d2 = 0;
			Matrix_dotProduct(In_col, Out_col, &d1); // <V_n, U_k>
			Matrix_dotProduct(Out_col, Out_col, &d2); // <U_k, U_k>
			d2 = 1/d2;
			for(unsigned int i = 0; i < In->row; ++i){
				temp->numbers[i] -= d1 * d2 * Out_col->numbers[i]; // sum( [<V_n, U_k>/<u_k, u_k>]*u_k )
			}
		}
		Matrix_magnitude(temp, &euclideannorm);
		euclideannorm = 1/euclideannorm;

		for(unsigned int i = 0; i < In->row; ++i){
			Out->numbers[i*In->col + j] = temp->numbers[i] * euclideannorm;
		}

	}
	Matrix_free(In_col);
	Matrix_free(Out_col);
	Matrix_free(temp);

}

int Matrix_orthonormalization(Matrix_t* In, Matrix_t* Out){
	if(In == NULLMTR || Out == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(Matrix_Internal_checkDimensionsExact(In, Out)!= MATRIX_OK){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);
		Matrix_Internal_orthogonalization(In, temp);
		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}
	Matrix_Internal_orthogonalization(In, Out);
	return MATRIX_OK;
}

int Matrix_descomposition_QR(Matrix_t* In, Matrix_t* Out_Q, Matrix_t* Out_R){
	//In->row >= In->col
	//Q is In->row x In->row;
	//R is In->row x In->col;
	if(In == NULLMTR|| Out_Q == NULLMTR|| Out_R == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->row < In->col || Matrix_Internal_checkDimensionsExact(In, Out_Q) ||
			Out_R->row != In->col || Out_R->col != In->col ){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	if(Matrix_Internal_checkAddress(In, Out_Q) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(In, Out_R) == MATRIX_SAME ||
			Matrix_Internal_checkAddress(Out_Q, Out_R) == MATRIX_SAME){
		return MATRIX_ERROR | MATRIX_SAME;
	}

	Matrix_t* In_col = Matrix_new(In->row, 1, 0);
	Matrix_t* Out_col = Matrix_new(In->row, 1, 0);
	Matrix_t* temp = Matrix_new(In->row, 1, 0);

	for(unsigned int j = 0; j < In->col; ++j){
		Matrix_submatrix(In, 0, In->row-1, j, j, In_col); // V_n
		Matrix_submatrix(In, 0, In->row-1, j, j, temp);
		double euclideannorm = 0;
		for(unsigned int k = 0; k < j; ++k){
			Matrix_submatrix(Out_Q, 0, Out_Q->row-1, k, k, Out_col); // U_k
			double d1 = 0, d2 = 0;
			Matrix_dotProduct(In_col, Out_col, &d1); // <V_n, U_k>
			Matrix_dotProduct(Out_col, Out_col, &d2); // <U_k, U_k>
			d2 = 1/d2;
			for(unsigned int i = 0; i < In->row; ++i){
				temp->numbers[i] -= d1 * d2 * Out_col->numbers[i]; // sum( [<V_n, U_k>/<u_k, u_k>]*u_k )
			}
		}
		Matrix_magnitude(temp, &euclideannorm);
		euclideannorm = 1/euclideannorm;

		for(unsigned int i = 0; i < In->row; ++i){
			Out_Q->numbers[i*In->col + j] = temp->numbers[i] * euclideannorm;
		}
		//calculate R
		for(unsigned int l = j; l < In->col; ++l){
			//as <k*U, V> = k*<U,V>
			double d = 0;
			Matrix_submatrix(In, 0, In->row-1, l, l, In_col);
			Matrix_dotProduct(temp, In_col , &d);
			Out_R->numbers[j*Out_R->col + l] = d*euclideannorm;
		}
	}
	Matrix_free(In_col);
	Matrix_free(Out_col);
	Matrix_free(temp);

	return MATRIX_OK;
}

int Matrix_setMatrixtoZero(Matrix_t* In){
	if(In == NULLMTR){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	for(unsigned int k = 0; k < In->row*In->col; ++k){
		In->numbers[k] = 0;
	}
	return MATRIX_OK;
}

int Matrix_compare(Matrix_t* In_1, Matrix_t* In_2){
	int i1 = Matrix_norm(In_1);
	int i2 = Matrix_norm(In_2);
	return (i1 > i2) - (i1 < i2);
}

/*

/int Matrix_eigenValuesVectors(Matrix_t* In, Matrix_t* Out_Eigenvalues_real, Matrix_t* Out_Eigenvalues_imaginary, Matrix_t* Out_Eigenvectors_real, Matrix_t* Out_Eigenvectors_imaginary){
	//In should be square
	//Eigenvectors should have the same dimension of In
	//Eigenvalues should be 1x
	return MATRIX_OK;

}
*/

//VECTOR

int Matrix_dotProduct(Matrix_t* In_1, Matrix_t* In_2, double* out){
	if(In_1 == NULLMTR || In_2 == NULLMTR || out == (double*)0){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In_1->col != 1 || In_2->col != 1 || In_1->row != In_2->row){
		return MATRIX_ERROR | MATRIX_VECTOR;
	}
	*out = 0;
	for(unsigned int k = 0; k<In_1->row/* *In_1->col */; ++k){
		*out += In_1->numbers[k] * In_2->numbers[k];
	}


	return MATRIX_OK;
}

int Matrix_magnitude(Matrix_t* In, double* out){
	if(In == NULLMTR || out == (double*)0){
		return MATRIX_ERROR | MATRIX_NULL;
	}
	if(In->col != 1){
		return MATRIX_ERROR | MATRIX_VECTOR;
	}
	Matrix_dotProduct(In, In, out);
	*out = sqrt(*out);
	return MATRIX_OK;
}

//PARALLEL
#ifdef __PARALLEL__
static void multiplication_thread(Matrix*t In_1, double* Out){

}

int Matrix_multiplication_forked(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	// [n x m] * [m x p] = [m x p]
	if(In_1->col != In_2->row || In_1->row != Out->row || In_2->col != Out->col){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	/* Each matrix should divide in 4 submatrixes [ A B
	 * 												C D ]
	 *
	 * the result of the multiplication is [A1*A2+B1*C2 A1*B2++B1*D2
	 * 										C1*A2+D1*C2 C1*B2+D1*D2]
	 *
	 * Matrix In_1 is n x m => A1 is
	 *
	 */

	//TODO: Forked Method
	return MATRIX_ERROR;
}
#endif // __PARALLEL__
