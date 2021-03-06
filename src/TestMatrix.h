/*
 * TestMatrix.h
 *
 *  Created on: Apr 13, 2021
 *      Author: jjfll
 */

#ifndef TESTMATRIX_H_
#define TESTMATRIX_H_

/*
	//Not being tested
	void Matrix_free(Matrix_t* In);
	Matrix_t* Matrix_random(unsigned int row, unsigned int col, double from, double to);
	Matrix_t* Matrix_random_int(unsigned int row, unsigned int col, int lower, int upper);



	int Matrix_scalarMultiplication(Matrix_t* In, double factor, Matrix_t* Out);
	int Matrix_addition(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);
	int Matrix_subtraction(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);
	int Matrix_multiplication(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);
	int Matrix_multiplication_recursive(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);

	//vectors
	int Matrix_dotProduct(Matrix_t* In_1, Matrix_t* In_2, double* out);
	int Matrix_magnitude(Matrix_t* In_1, double* out);

	// Utilities
	void Matrix_display(Matrix_t* In);

	Matrix_t* Matrix_copy(Matrix_t* In);
	int Matrix_copyTo(Matrix_t* In, Matrix_t* Out);

	int Matrix_submatrix(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right, Matrix_t* Out);

	int Matrix_removeRow(Matrix_t* In, unsigned int row, Matrix_t* Out);
	int Matrix_removeColumn(Matrix_t* In, unsigned int col, Matrix_t* Out);
	int Matrix_transpose(Matrix_t* In, Matrix_t* Out);
	int Matrix_determinant(Matrix_t* In, double* out);
	int Matrix_trace(Matrix_t* In, double* out);

	int Matrix_adjoint(Matrix_t* In, Matrix_t* Out);
	int Matrix_inverse(Matrix_t* In, Matrix_t* Out);
	int Matrix_inverse_triangular(Matrix_t* In, Matrix_t* Out);

	Matrix_t* Matrix_rowEchelon(Matrix_t* In);
	Matrix_t* Matrix_reducedRowEchelon(Matrix_t* In);

	int Matrix_hConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);
	int Matrix_vConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out);
	double Matrix_norm(Matrix_t* In);
	Matrix_t* Matrix_kernel(Matrix_t* In);
	Matrix_t* Matrix_nullSpace(Matrix_t* In);

	int Matrix_compare(Matrix_t* In_1, Matrix_t* In_2);

	int Matrix_orthonormalization(Matrix_t* In, Matrix_t* Out);
	int Matrix_descomposition_LU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_U);
	int Matrix_descomposition_LDU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_D, Matrix_t* Out_U);

	int Matrix_setMatrixtoZero(Matrix_t* In);
	//TODO
	int Matrix_descomposition_QR(Matrix_t* In, Matrix_t* Out_Q, Matrix_t* Out_R);
	int Matrix_eigenValuesVectors(Matrix_t* In, Matrix_t* Out_Eigenvalues_real, Matrix_t* Out_Eigenvalues_imaginary, Matrix_t* Out_Eigenvectors_real, Matrix_t* Out_Eigenvectors_imaginary);

 */
#include <stdlib.h>
#include <stdio.h>
#include "Matrix.h"


#define TEST_PASS printf("PASSED\n");
#define TEST_FAIL printf("FAILED (Line:%d)\n", __LINE__);
#define INTRO printf("%s \t", __FUNCTION__ );
#define __TODO__  printf("\tTo do:%s\n", __FUNCTION__);


#include <math.h>
//1 if a is more, -1 if b is more, 0 if they are FLT_EPSILON appart
static int comparefloat(double a, double b){
	//fabs(In->numbers[j]-0.0) <= FLT_EPSILON)
	if(fabs(a - b) <= FLT_EPSILON){
		return 0;
	}
	return (a > b) - (a < b);
}

void Test_Matrix_new(){
	//New Matrix
	INTRO

	unsigned int row = 3;
	unsigned int col = 4;
	double val = 1.5;

	Matrix_t* M = Matrix_new(row,col, val);
	//expected new matrix with  elements equal to 1
	if(M->row != row || M->col != col){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	for(unsigned int k = 0; k < row*col; ++k){
		if(M->numbers[k] != val){
			Matrix_free(M);
			TEST_FAIL;
			return;
		}
	}

	Matrix_t* N = Matrix_new(-1, -1, 1);
	if(N != NULL){
		Matrix_free(N);
		TEST_FAIL;
		return;
	}

	Matrix_free(M);
	Matrix_free(N);
	TEST_PASS;
	return;
}

void Test_Matrix_identity(){
	INTRO
	unsigned int n = 4;
	Matrix_t* M = Matrix_identity(n);
	if(M->row != n || M->col != n){
		Matrix_free(M);
		TEST_FAIL
		return;
	}
	for(unsigned int i = 0; i < n; ++i){
		for(unsigned int j = 0; j < n; ++j){
			if(i != j){
				if(M->numbers[i*n+j] != 0){
					Matrix_free(M);
					TEST_FAIL;
					return;
				}
			} else {
				if(M->numbers[i*n+j] != 1){
					Matrix_free(M);
					TEST_FAIL;
					return;
				}
			}
		}
	}

	Matrix_t* N = Matrix_identity(-1);
	if(N != NULL){
		Matrix_free(N);
		TEST_FAIL
		return;
	}

	Matrix_free(M);
	Matrix_free(N);
	TEST_PASS;
	return;

}

void Test_Matrix_zeroes(){
		//New Matrix
		INTRO

		unsigned int row = 4;
		unsigned int col = 5;
		Matrix_t* M = Matrix_zeroes(row,col);
		//expected new matrix with  elements equal to 1
		if(M->row != row || M->col != col){
			Matrix_free(M);
			TEST_FAIL
			return;
		}
		for(unsigned int k = 0; k < row*col; ++k){
			if(M->numbers[k] != 0){
				Matrix_free(M);
				TEST_FAIL
				return;
			}
		}


		Matrix_t* N = Matrix_zeroes(-1, 1);
		if(N != NULL){
			Matrix_free(N);
			TEST_FAIL
			return;
		}

		Matrix_free(M);
		Matrix_free(N);
		TEST_PASS;
		return;
}

void Test_Matrix_ones(){
		//New Matrix
		INTRO;

		unsigned int row = 3;
		unsigned int col = 2;
		Matrix_t* M = Matrix_ones(row,col);
		//expected new matrix with  elements equal to 1
		if(M->row != row || M->col != col){
			Matrix_free(M);
			TEST_FAIL
			return;
		}

		for(unsigned int k = 0; k < row*col; ++k){
			if(M->numbers[k] != 1){
				Matrix_free(M);
				TEST_FAIL
				return;
			}
		}

		Matrix_t* N = Matrix_ones(-1, 2);
		if(N != NULL){
			Matrix_free(N);
			TEST_FAIL
			return;
		}

		Matrix_free(M);
		Matrix_free(N);
		TEST_PASS
		return;
}

void Test_Matrix_get(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val = 1.5;
	Matrix_t* M = Matrix_new(row,col, val);


	double get = 0;

	//Out of bounds get
	if(Matrix_get(M, row, col-1, &get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_get(M,row-1, col, &get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_get(M,row-1, -2, &get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_get(M, -2, col-1, &get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	//in bounds get
	for(unsigned int i = 0; i < row; ++i){
		for(unsigned int j = 0; j < col; ++j){
			if(Matrix_get(M, i, j, &get) != MATRIX_OK){
				Matrix_free(M);
				TEST_FAIL;
				return;
			}
			if(get != val){
				Matrix_free(M);
				TEST_FAIL;
				return;
			}
		}
	}

	Matrix_free(M);
	TEST_PASS;

}

void Test_Matrix_set(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val = 3.5;
	Matrix_t* M = Matrix_new(row,col, val);

	double get = 15;

	//out of bounds set
	if(Matrix_set(M, row, col-1, get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_set(M,row-1, col, get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_set(M,row-1, -2, get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	if(Matrix_set(M, -2, col-1, get) == MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	// in Bound set

	if(Matrix_set(M, (unsigned int)((row - 1)/2), (unsigned int)((col - 1)/2), get) != MATRIX_OK){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}
	if(M->numbers[(row - 1)/2 * M->col + (col - 1)/2] != get){
		Matrix_free(M);
		TEST_FAIL;
		return;
	}

	Matrix_free(M);
	TEST_PASS;

}

void Test_Matrix_scalarMultiplication(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val = 3.5;
	Matrix_t* M = Matrix_new(row,col, val);
	Matrix_t* N = Matrix_zeroes(row,col);
	Matrix_t* O = Matrix_zeroes(col, row);


	//wrong output size
	if(Matrix_scalarMultiplication(M,val,O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		TEST_FAIL;
		return;
	}

	//Right output size
	if(Matrix_scalarMultiplication(M,val,N) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		TEST_FAIL;
		return;
	}

	for(unsigned int k = 0; k < N->row*N->col; ++k){
		if(comparefloat(N->numbers[k], val*M->numbers[k]) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			TEST_FAIL;
			return;
		}
	}

	//Input = output
	//copy N as M
	if(Matrix_scalarMultiplication(M,1,N) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		TEST_FAIL;
		return;
	}

	//N = val * N
	if(Matrix_scalarMultiplication(N,val,N) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		TEST_FAIL;
		return;
	}

	//check if done correctly
	for(unsigned int k = 0; k < N->row*N->col; ++k){
		if(comparefloat(N->numbers[k], val*M->numbers[k]) != 0){



			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			TEST_FAIL;
			return;
		}
	}

	Matrix_free(M);
	Matrix_free(N);
	Matrix_free(O);
	TEST_PASS;
	return;
}

void Test_Matrix_addition(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val1 = 3.5, val2 = 5.25;
	Matrix_t* M = Matrix_new(row, col, val1);
	Matrix_t* N = Matrix_new(row, col, val2);
	Matrix_t* O = Matrix_zeroes(row, col);
	Matrix_t* Q = Matrix_new(row, col + 1, val1);

	//Wreong size input
	if(Matrix_addition(Q, N, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	if(Matrix_addition(M, Q, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//wrong size output
	if(Matrix_addition(M, N, Q) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//right op
	if(Matrix_addition(M, N, O) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < O->row*O->col; ++k){
		if(comparefloat(O->numbers[k], val1 + val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	//input same as output
	if(Matrix_addition(M, N, N) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < N->row*N->col; ++k){
		if(comparefloat(N->numbers[k], val1 + val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	if(Matrix_addition(M, N, M) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < M->row*M->col; ++k){
		//N is val1 + val2, M is val1
		if(comparefloat(M->numbers[k], val1 + val1 + val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	Matrix_free(M);
	Matrix_free(N);
	Matrix_free(O);
	Matrix_free(Q);
	TEST_PASS;
	return;
}

void Test_Matrix_subtraction(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val1 = 3.5, val2 = 5.25;
	Matrix_t* M = Matrix_new(row, col, val1);
	Matrix_t* N = Matrix_new(row, col, val2);
	Matrix_t* O = Matrix_zeroes(row, col);
	Matrix_t* Q = Matrix_new(row, col + 1, val1);

	//Wreong size input
	if(Matrix_subtraction(Q, N, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	if(Matrix_subtraction(M, Q, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//wrong size output
	if(Matrix_subtraction(M, N, Q) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//right op
	if(Matrix_subtraction(M, N, O) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < O->row*O->col; ++k){
		if(comparefloat(O->numbers[k], val1 - val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	//input same as output
	if(Matrix_subtraction(M, N, N) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < N->row*N->col; ++k){
		if(comparefloat(N->numbers[k], val1 - val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	if(Matrix_subtraction(M, N, M) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < M->row*M->col; ++k){
		//N is (val1 - val2), M is val1
		if(comparefloat(M->numbers[k], val1 - (val1 - val2)) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	Matrix_free(M);
	Matrix_free(N);
	Matrix_free(O);
	Matrix_free(Q);
	TEST_PASS;
	return;
}

void Test_Matrix_multiplication(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val1 = 3.5, val2 = 5.25;
	Matrix_t* M = Matrix_new(row, col, val1);
	Matrix_t* N = Matrix_new(col, row, val2);
	Matrix_t* O = Matrix_zeroes(row, row);
	Matrix_t* Q = Matrix_new(row+1, col+1, val1);

	//Wreong size input
	if(Matrix_multiplication(Q, N, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	if(Matrix_multiplication(M, Q, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//wrong size output
	if(Matrix_multiplication(M, N, Q) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//right op
	if(Matrix_multiplication(M, N, O) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < O->row*O->col; ++k){
		if(comparefloat(O->numbers[k], col * val1 * val2) != 0){
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	Matrix_free(M);
	Matrix_free(N);
	Matrix_free(O);
	Matrix_free(Q);

	//square matrixes
	Matrix_t* P = Matrix_new(col, col, val1);
	Matrix_t* R = Matrix_new(col, col, val2);

	if(Matrix_multiplication(P, R, R) != MATRIX_OK){
		Matrix_free(P);
		Matrix_free(R);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < R->row*R->col; ++k){
		if(comparefloat(R->numbers[k], col * val1 * val2) != 0){
			Matrix_free(P);
			Matrix_free(R);
			TEST_FAIL;
			return;
		}
	}

	if(Matrix_multiplication(P, P, P) != MATRIX_OK){
		Matrix_free(P);
		Matrix_free(R);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < P->row*P->col; ++k){
		if(comparefloat(P->numbers[k], col * val1 * val1) != 0){

			Matrix_free(P);
			Matrix_free(R);
			TEST_FAIL;
			return;
		}
	}
	Matrix_free(P);
	Matrix_free(R);
	TEST_PASS;
	return;
}

void Test_Matrix_multiplication_recursive(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 5;
	double val1 = 3.5, val2 = 5.25;
	Matrix_t* M = Matrix_new(row, col, val1);
	Matrix_t* N = Matrix_new(col, row, val2);
	Matrix_t* O = Matrix_zeroes(row, row);
	Matrix_t* Q = Matrix_new(row+1, col+1, val1);

	//Wreong size input
	if(Matrix_multiplication_recursive(Q, N, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	if(Matrix_multiplication_recursive(M, Q, O) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//wrong size output
	if(Matrix_multiplication_recursive(M, N, Q) == MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	//right op
	if(Matrix_multiplication_recursive(M, N, O) != MATRIX_OK){
		Matrix_free(M);
		Matrix_free(N);
		Matrix_free(O);
		Matrix_free(Q);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < O->row*O->col; ++k){
		if(comparefloat(O->numbers[k], col * val1 * val2) != 0){
			Matrix_display(O);
			Matrix_free(M);
			Matrix_free(N);
			Matrix_free(O);
			Matrix_free(Q);
			TEST_FAIL;
			return;
		}
	}

	Matrix_free(M);
	Matrix_free(N);
	Matrix_free(O);
	Matrix_free(Q);

	//square matrixes
	Matrix_t* P = Matrix_new(col, col, val1);
	Matrix_t* R = Matrix_new(col, col, val2);

	if(Matrix_multiplication_recursive(P, R, R) != MATRIX_OK){
		Matrix_free(P);
		Matrix_free(R);
		TEST_FAIL;
		return;
	}

	for(unsigned int k = 0; k < R->row*R->col; ++k){

		if(comparefloat(R->numbers[k], col * val1 * val2) != 0){
			Matrix_display(R);
			Matrix_free(P);
			Matrix_free(R);
			TEST_FAIL;
			return;
		}
	}

	if(Matrix_multiplication_recursive(P, P, P) != MATRIX_OK){
		Matrix_free(P);
		Matrix_free(R);
		TEST_FAIL;
		return;
	}

	for(unsigned int k = 0; k < P->row*P->col; ++k){
		if(comparefloat(P->numbers[k], col * val1 * val1) != 0 ){
			Matrix_free(P);
			Matrix_free(R);
			TEST_FAIL;
			return;
		}
	}
	Matrix_free(P);
	Matrix_free(R);
	TEST_PASS;
	return;
}

void Test_Matrix_dotProduct(){
	INTRO
	__TODO__
}
void Test_Matrix_magnitude(){
	INTRO
	__TODO__
}
void Test_Matrix_display(){
	INTRO
	__TODO__
}
void Test_Matrix_copy(){
	INTRO
	__TODO__
}
void Test_Matrix_copyTo(){
	INTRO
	__TODO__
}
void Test_Matrix_submatrix(){
	INTRO
	__TODO__
}
void Test_Matrix_removeRow(){
	INTRO
	__TODO__
}
void Test_Matrix_removeColumn(){
	INTRO
	__TODO__
}
void Test_Matrix_transpose(){
	INTRO
	__TODO__
}
void Test_Matrix_determant(){
	INTRO
	__TODO__
}
void Test_Matrix_trace(){
	INTRO
	__TODO__
}
void Test_Matrix_adjoint(){
	INTRO
	__TODO__
}
void Test_Matrix_inverse(){
	INTRO
	__TODO__
}
void Test_Matrix_inverse_triangular(){
	INTRO
	__TODO__
}
void Test_Matrix_rowEchelon(){
	INTRO
	__TODO__
}
void Test_Matrix_reducedRowEchelon(){
	INTRO
	__TODO__
}
void Test_Matrix_hConcat(){
	INTRO
	__TODO__
}
void Test_Matrix_vConcat(){
	INTRO
	__TODO__
}
void Test_Matrix_norm(){
	INTRO
	__TODO__
}
void Test_Matrix_kernel(){
	INTRO
	__TODO__
}
void Test_Matrix_nullSpace(){
	INTRO
	__TODO__
}
void Test_Matrix_compare(){
	INTRO
	__TODO__
}
void Test_Matrix_orthonormalization(){
	INTRO
	__TODO__
}
void Test_Matrix_descomposition_LU(){
	INTRO
	__TODO__
}
void Test_Matrix_descomposition_LDU(){
	INTRO
	__TODO__
}
void Test_Matrix_setMatrixtoZero(){
	INTRO
	unsigned int row = 3;
	unsigned int col = 4;
	double val = 1.5;
	Matrix_t* A = Matrix_new(row, col, val);
	if(Matrix_setMatrixtoZero(A) != MATRIX_OK){
		Matrix_free(A);
		TEST_FAIL;
		return;
	}
	for(unsigned int k = 0; k < A->row*A->col; ++k){
		if(comparefloat(A->numbers[k], 0) != 0){
			Matrix_free(A);
			TEST_FAIL;
			return;
		}
	}
	TEST_PASS;
	return;
}



#endif /* TESTMATRIX_H_ */
