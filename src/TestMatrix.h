/*
 * TestMatrix.h
 *
 *  Created on: Apr 13, 2021
 *      Author: jjfll
 */

#ifndef TESTMATRIX_H_
#define TESTMATRIX_H_

/*
	Matrix_t* Matrix_random(unsigned int row, unsigned int col, double from, double to);
	Matrix_t* Matrix_random_int(unsigned int row, unsigned int col, int lower, int upper);
	//De-constructor
	void Matrix_free(Matrix_t* In);
	//Getter
	int Matrix_get(Matrix_t* In, unsigned int row, unsigned int col, double* out);
	//Setter
	int Matrix_set(Matrix_t* In, unsigned int row, unsigned int col, double val);
	//operations
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

#define TEST_PASS printf("%s PASSED\n", __FUNCTION__);
#define TEST_FAIL printf("%s FAILED (%d)\n", __FUNCTION__, __LINE__);
#define INTRO printf("%s ... ", __FUNCTION__ );

void Test_Matrix_new(){
	//New Matrix
	INTRO

	unsigned int row = 3;
	unsigned int col = 4;
	double val = 1.5;

	Matrix_t* M = Matrix_new(row,col, val);
	//expected new matrix with  elements equal to 1
	if(M->row != row || M->col != col){
		TEST_FAIL
		return;
	}

	for(unsigned int k = 0; k < row*col; ++k){
		if(M->numbers[k] != val){
			TEST_FAIL
			return;
		}
	}

	Matrix_t* N = Matrix_new(-1, -1, 1);
	if(N != NULL){
		TEST_FAIL
		return;
	}

	TEST_PASS;
	//Matrix_free(M);
	return;
}

void Test_Matrix_identity(){
	INTRO
	unsigned int n = 4;
	Matrix_t* M = Matrix_identity(n);
	if(M->row != n || M->col != n){
		TEST_FAIL
		return;
	}
	for(unsigned int i = 0; i < n; ++i){
		for(unsigned int j = 0; j < n; ++j){
			if(i != j){
				if(M->numbers[i*n+j] != 0){
					TEST_FAIL;
					return;
				}
			} else {
				if(M->numbers[i*n+j] != 1){
					TEST_FAIL;
					return;
				}
			}
		}
	}
	TEST_PASS;
	Matrix_free(M);
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
			TEST_FAIL
			return;
		}

		for(unsigned int k = 0; k < row*col; ++k){
			if(M->numbers[k] != 0){
				TEST_FAIL
				return;
			}
		}
		TEST_PASS;
		Matrix_free(M);
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
			TEST_FAIL
			return;
		}

		for(unsigned int k = 0; k < row*col; ++k){
			if(M->numbers[k] != 1){
				TEST_FAIL
				return;
			}
		}
		TEST_PASS
		Matrix_free(M);
		return;
}



#endif /* TESTMATRIX_H_ */
