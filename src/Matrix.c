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
	if(In_1->row != In_2->row || In_1->col != In_2->col){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkIndex(Matrix_t* In, int row, int col){/* checks if row and col are valid indexes for matrix M */
	if(row > In->row || col > In->col){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkSquare(Matrix_t* In){
	if(In->row != In->col){
		return MATRIX_ERROR;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_checkTriangularUpper(Matrix_t* In){
	if(In->row == 1){
		return MATRIX_OK;
	}

	char flag = 0x00;
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < i; ++j){
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
	return MATRIX_DIFFERENT;
}

static int Matrix_Internal_checkTriangularLower(Matrix_t* In){
	if(In->row == 1){
		return MATRIX_OK;
	}
	char flag = 0x00;
	for(int i = 0; i < In->row; ++i){
		for(int j = i + 1; j < In->col; ++j){
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
	return (CheckTriangularUpper(In) == MATRIX_OK && CheckTriangularLower(In) == MATRIX_OK) ? MATRIX_OK : MATRIX_ERROR;
}
*/

//Constructors
Matrix_t* Matrix_new(unsigned int row, unsigned int col, double val){
	Matrix_t* Out = (Matrix_t *)malloc(sizeof(Matrix_t));
	Out->row = row;
	Out->col = col;
	Out->numbers = (double*)malloc(sizeof(double)*row*col);
	int k = 0;
	for(int i = 0; i < (Out->row * Out->col); ++i){
		Out->numbers[k++] = val;
	}
	return Out;
}

void Matrix_free(Matrix_t* In){
	free(In->numbers);
	free(In);
}

Matrix_t* Matrix_identity(unsigned int n){
	Matrix_t* I = Matrix_new(n, n, 0);
	for(int i = 0; i < n; ++i){
		I->numbers[i*I->col + i] = 1;
	}
	return I;
}

Matrix_t* Matrix_zeroes(unsigned int row, unsigned int col){
	Matrix_t* A = Matrix_new(row, col, 0);
	return A;
}

Matrix_t* Matrix_ones(unsigned int row, unsigned int col){
	Matrix_t* A = Matrix_new(row, col, 1);
	return A;
}

Matrix_t* Matrix_random(unsigned int row, unsigned int col, double from, double to){
	Matrix_t* R = Matrix_new(row, col, 0);
	int k = 0;
	for(int i = 0; i < row * col ; ++i){
		double r = ((double)rand()) / ((double)RAND_MAX);
		R->numbers[k++] = from + (to-from)*r;
	}
	return R;
}

Matrix_t* Matrix_random_int(unsigned int row, unsigned int col, int from, int to){
	Matrix_t* R = Matrix_new(row, col, 0);
	int k = 0;
	for(int i = 0; i < row * col ; ++i){
		double r = ((double)rand()) / ((double)RAND_MAX);
		R->numbers[k++] = (int)(from + (to-from)*r);
	}
	return R;
}

int Matrix_get(Matrix_t* In, unsigned int row, unsigned int col, double *out){
	if(Matrix_Internal_checkIndex(In, row, col) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_INDEX;
	}
	*out = In->numbers[row * In->col + col];
	return MATRIX_OK;
}

int Matrix_set(Matrix_t* In, unsigned int row, unsigned int col, double val){
	if(Matrix_Internal_checkIndex(In, row, col) != MATRIX_OK){
			return MATRIX_ERROR | MATRIX_INDEX;
		}
	In->numbers[ row * In->col + col] = val;
	return MATRIX_OK;
}

//Algebra
int Matrix_scalarMultiplication(Matrix_t* In, double factor, Matrix_t* Out){
	if(Matrix_Internal_checkDimensionsExact(In, Out)){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	for(int k = 0; k < In->row * In->col; ++k){
		Out->numbers[k] = In->numbers[k] * factor;
	}
	return MATRIX_OK;
}

int Matrix_addition(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(Matrix_Internal_checkDimensionsExact(In_1, In_2) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In_1, Out) != MATRIX_OK){
			return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	for(int k = 0; k < In_1->row * In_1->col; ++k){
		Out->numbers[k] = In_1->numbers[k] + In_2->numbers[k];
	}
	return MATRIX_OK;
}

int Matrix_subtraction(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(Matrix_Internal_checkDimensionsExact(In_1, In_2) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In_1, Out) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}

	for(int k = 0; k < In_1->row * In_1->col; ++k){
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

	for(int i = 0; i < In_1->row; ++i){
		for(int j = 0; j < In_2->col; ++j){
			for(int k = 0; k < In_2->row;++k){
				Out->numbers[i * Out->col + j] += In_1->numbers[i * (In_1->col) + k] * In_2->numbers[k * (In_2->col) + j];
			}
		}
	}

	return;
}

int Matrix_multiplication(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
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

	// base cases => n = 1 (p = 1) => multiply vector
	//		else partition
	if(a_left == a_right || a_upper == a_lower || b_left == b_right || b_upper == b_lower){

		for(int i = 0; i < (a_lower - a_upper + 1); ++i){
			for(int j = 0; j < (b_right-b_left + 1); ++j){
				for(int k = 0; k < (b_lower - b_upper + 1);++k){
					C->numbers[(a_upper + i)*C->col + (b_left + j)] +=
							A->numbers[(a_upper + i) * (A->col) + (a_right + k)] * B->numbers[(b_upper + k) * (B->col) + (b_right + j)];
				}
			}
		}
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
	for(int i = 0; i < In->col; ++i){
		for(int j = 0; j < In->row; ++j){
			Out->numbers[k++] = In->numbers[j * In->row + i];
		}
	}
}

int Matrix_transpose(Matrix_t* In, Matrix_t* Out){
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
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}

	if(In->row == 1){
		*out = In->numbers[0];
		return MATRIX_OK;
	}

	Matrix_t* M1 = Matrix_removeRow_2(In, 0);
	Matrix_t* M2 = Matrix_zeroes(In->row-1, In->col-1);
	double d = 0, sign = +1;
	for(int j = 0; j < In->col; ++j){
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
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}
	*out = 0;
	for(int i = 0; i < In->col; ++i){
		*out += In->numbers[i*In->row+i];
	}
	return MATRIX_OK;
}

static void Matrix_Internal_adjoint(Matrix_t* In, Matrix_t* Out){
	Matrix_t* A1 = Matrix_new(In->row - 1, In->col, 0);
	Matrix_t* A2 = Matrix_new(In->row - 1, In->col - 1, 0);

	for(int i = 0; i < In->row;++i){
		Matrix_removeRow(In,i,A1);
		for(int j = 0; j < In->col;++j){
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

static int Matrix_internal_inverse(Matrix_t* In, Matrix_t* Out){
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
	if( (Matrix_Internal_checkSquare(In) != MATRIX_OK) || Matrix_Internal_checkDimensionsExact(In, Out) != MATRIX_OK || (Matrix_Internal_checkSquare(Out) != MATRIX_OK) ){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}
	if(Matrix_Internal_checkAddress(In, Out) == MATRIX_SAME){
		Matrix_t* temp = Matrix_zeroes(Out->row, Out->col);

		if( Matrix_internal_inverse(In, temp) != MATRIX_OK){
			Matrix_free(temp);
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}

		Matrix_copyTo(temp, Out);
		Matrix_free(temp);
		return MATRIX_OK;
	}
	if( Matrix_internal_inverse(In, Out) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
	}
	return MATRIX_OK;
}

static int Matrix_Internal_inverse_triangular_upper(Matrix_t* In, Matrix_t* Out){

	for(int i = 0; i < Out->row; ++i){
		//if element in diagonal is 0, then non invertible
		if(In->numbers[i * (In->col) + i] == 0){
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}
		for(int j = 0; j < Out->col; j++){
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

	for(int i = 0; i < Out->row; ++i){
		//if element in diagonal is 0, then non invertible
		if(In->numbers[i * (In->col) + i] == 0){
			return MATRIX_ERROR | MATRIX_NON_INVERTIBLE;
		}
		for(int j = 0; j < Out->col; j++){
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
	//__TODO__ case lower triangular

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
	if(In->row > 0 && In->col > 0){
		int k = 0;
		printf("[");
		for(int i = 0; i < In->row; ++i){
			for(int j = 0; j < In->col; ++j){
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
	Matrix_t* Out = Matrix_new(In->row, In->col, 0);
	for(int k = 0; k < In->row * In->col; k++){
		Out->numbers[k]=In->numbers[k];
	}
	return Out;
}

int Matrix_copyTo(Matrix_t* In, Matrix_t* Out){
	if(Matrix_Internal_checkDimensionsExact(In, Out)){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	for(int k = 0; k < In->row * In->col; k++){
			Out->numbers[k]=In->numbers[k];
		}
	return MATRIX_OK;
}

int Matrix_submatrix(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right, Matrix_t* Out){
	if(lower < upper || right <  left || Matrix_Internal_checkIndex(In, lower, right) != MATRIX_OK){
		return MATRIX_ERROR || MATRIX_INDEX;
	}
	if(Out->row !=  lower-upper+1 || Out->col != right-left+1){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}

	int k = 0;
	for(int i = upper; i <= lower; ++i){
		for(int j = left; j <= right; ++j){
			Out->numbers[k++] = In->numbers[i * (In->col) + j];
		}
	}
	return MATRIX_OK;
}

int Matrix_removeRow(Matrix_t* In, unsigned int row, Matrix_t* Out){
	if(In->row != Out->row + 1 || In->col != Out->col ){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < In->col; ++j){
			if(i != row){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}
	return MATRIX_OK;
}

int Matrix_removeColumn(Matrix_t* In, unsigned int col, Matrix_t* Out){
	if(In->row != Out->row || In->col != Out->col + 1 ){
		return MATRIX_ERROR | MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < In->col; ++j){
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
		for(int j = 0; j < In->col;++j){
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
	int ind1 = B->col;
	int ind2 = 0;
	//search first not zero entry
	for(int i = 0; i < B->row; ++i){
		for(int j = 0; j < B->col; ++j){
			if(B->numbers[i*B->col+j]!=0 && j < ind1){
				ind1 = j;
				ind2 = i;
				break;
			}
		}
	}
	if(ind2 > 1){
		//if I am not in the first row, swap it with the first
		for(int j = 0; j < B->col; ++j){
			double temp = B->numbers[j];
			B->numbers[j] = B->numbers[ind2*B->col+j];
			B->numbers[ind2*B->col+j] = temp;
		}
	}
	if(B->numbers[0]!=0){
		//normalize first row so the first digit is 1 (or 0)
		double coef = B->numbers[0];
		for(int j = 0; j < B->col; ++j){
			B->numbers[j] /= coef;
		}
		//Row = row - row_coef * row_0
		for(int i = 1; i < B->row; ++i){
			coef = B->numbers[i*B->col];
			for(int j = 0; j < B->col; ++j){
				B->numbers[i*B->col+j] -= coef * B->numbers[j];
			}
		}
	} else {
		//normalize first row
		double coef = 0;
		for(int j = 0; j < B->col;++j){
			if(B->numbers[j] != 0 && coef == 0){
				coef = B->numbers[j];
				B->numbers[j] = 1;
			} else if(B->numbers[j]!= 0){
				B->numbers[j] /= coef;
			}
		}
		for(int i = 1; i < B->row; ++i){
			coef = B->numbers[i*B->col];
			for(int j = 0; j < B->col; ++j){
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
	Matrix_t* B1 = Matrix_removeRow_2(B,0);
	Matrix_t* Be = Matrix_rowEchelon(B1);

	//copy results from Be to B
	for(int i = 0; i < Be->row; ++i){
		for(int j = 0; j < Be->col; ++j ){
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
		for(int j = 0; j < R->col; ++j){
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
		for(int j = 0; j < R->col; ++j){
			R->numbers[i*R->col + j] -= coef*R->numbers[idx1*R->col + j];
		}
	}
	//eliminate the lower n rows and calculate recursively
	Matrix_t* R1 = Matrix_submatrix_2(R, 0, idx1-1, 0, R->col-1);
	Matrix_t* R2 = Matrix_Internal_reducedRowEchelon(R1);
	Matrix_free(R1);
	//copy values
	for(int i = 0; i < idx1; ++i){
		for(int j = 0; j < R->col; ++j){
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
	if(In_1->col + In_2->col != Out->col || In_1->row != In_2->row || In_1->row != Out->row){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(int i = 0; i < In_1->row; ++i){
		for(int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
		for(int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return MATRIX_OK;
}

int Matrix_vConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out){
	if(In_1->row + In_2->row != Out->row || In_1->col != In_2->col || In_1->col != Out->col){
		return MATRIX_ERROR || MATRIX_DIMENSIONS;
	}
	int k = 0;
	for(int i = 0; i < In_1->row; ++i){
		for(int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
	}
	for(int i = 0; i < In_2->row; ++i){
		for(int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return MATRIX_OK;
}

double Matrix_norm(Matrix_t* In){
	double d = 0;
	for(int k = 0; k < In->row*In->col; ++k){
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
		for(int j =0; j < RM->col; ++j){
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

	Matrix_t* RRM = Matrix_submatrix_2(RM, 0, k, 0, RM->col-1);
	Matrix_free(RM);

	unsigned int nn = RRM->col - RRM->row;

	if(nn == 0){
		Matrix_t* N = Matrix_new(0,0,0);
		return N;
	}

	//obtain free columns R2
	Matrix_t* F = Matrix_submatrix_2(RRM, 0, RRM->row - 1, RRM->row, RRM->col-1);
	Matrix_free(RRM);


	Matrix_t* I = Matrix_identity(nn);

	//reciprocate F
	Matrix_scalarMultiplication(F, -1, F);

	//concat with identity
	Matrix_t* N = Matrix_vConcat_2(F,I);
	Matrix_free(F);
	Matrix_free(I);

	return N;
}

int Matrix_descomposition_LU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_U){
	//check squareness
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In, Out_L) != MATRIX_OK || Matrix_Internal_checkDimensionsExact(In, Out_U) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
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

	Matrix_t* w  = Matrix_submatrix_2(In, 0, 0, 1, In->col-1);
	Matrix_t* v  = Matrix_submatrix_2(In, 1, In->row-1, 0,0);
	Matrix_t* Ab = Matrix_submatrix_2(In, 1, In->row-1, 1, In->col-1);
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
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < In->col; ++j){
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
	if(Matrix_Internal_checkSquare(In) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_L) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_D) != MATRIX_OK ||
			Matrix_Internal_checkDimensionsExact(In, Out_U) != MATRIX_OK){
		return MATRIX_ERROR | MATRIX_SQUARE;
	}



	Matrix_descomposition_LU(In, Out_L, Out_U);
	//Matrix U has the coeficients of Matrix D
	//must assure that matrix D has zeros in non diagonal index;
	for(int i = 0; i < In->col; ++i){
		for(int j = 0; j < In->row; ++j){
			if(i == j){
				Out_D->numbers[i*In->col+j] = Out_U->numbers[i*In->col+j];
				Out_U->numbers[i*In->col+j] = 1;
			}
		}
	}
	return MATRIX_OK;
}

//static double Matrix_Internal_innermultiply(Matrix_t* A, Matrix_t* B);
//int Matrix_descomposition_QR(Matrix_t* M, Matrix_t* Q, Matrix_t* R);

int Matrix_setMatrixtoZero(Matrix_t* In){
	for(int k = 0; k < In->row*In->col; ++k){
		In->numbers[k] = 0;
	}
	return MATRIX_OK;
}


int Matrix_compare(Matrix_t* In_1, Matrix_t* In_2){
	int i1 = Matrix_norm(In_1);
	int i2 = Matrix_norm(In_2);
	return (i1 > i2) - (i1 < i2);
}

//Please check input size
Matrix_t* Matrix_transpose_2(Matrix_t* In){
	Matrix_t* Out = Matrix_new(In->col, In->row, 0);
	int k = 0;
	for(int i = 0; i < In->col; ++i){
		for(int j = 0; j < In->row; ++j){
			Out->numbers[k++] = In->numbers[j * In->row + i];
		}
	}
	return Out;
}

Matrix_t* Matrix_submatrix_2(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right){
	if(lower < upper || right <  left || Matrix_Internal_checkIndex(In, lower, upper) != MATRIX_OK){
		return Matrix_new(0,0,0);
	}

	Matrix_t* Out = Matrix_new(lower - upper + 1, right - left + 1, 0);
	int k = 0;
	for(int i = upper; i <= lower; ++i){
		for(int j = left; j <= right; ++j){
			Out->numbers[k++] = In->numbers[i * (In->col) + j];
		}
	}
	return Out;
}

Matrix_t* Matrix_removeRow_2(Matrix_t* In, unsigned int row){
	Matrix_t* Out = Matrix_new(In->row-1, In->col, 0);
	int k = 0;
	for(int i = 0; i < In->row; ++i){
		for(int j = 0; j < In->col; ++j){
			if(i != row){
				Out->numbers[k++] = In->numbers[i * In->col + j];
			}
		}
	}
	return Out;
}

Matrix_t* Matrix_removeColumn_2(Matrix_t* In, unsigned int col){
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

Matrix_t* Matrix_hConcat_2(Matrix_t* In_1, Matrix_t* In_2){
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

Matrix_t* Matrix_vConcat_2(Matrix_t* In_1, Matrix_t* In_2){
	Matrix_t* Out = Matrix_new(In_1->row + In_2->row, In_1->col, 0);
	int k = 0;
	for(int i = 0; i < In_1->row; ++i){
		for(int j = 0; j < In_1->col; ++j){
			Out->numbers[k++] = In_1->numbers[i*In_1->col+j];
		}
	}
	for(int i = 0; i < In_2->row; ++i){
		for(int j = 0; j < In_2->row; ++j){
			Out->numbers[k++] = In_2->numbers[i*In_2->col+j];
		}
	}
	return Out;
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
