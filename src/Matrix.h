/*
 ============================================================================
 Name        :	Matrix.h
 Author      :	Jaime Feller
 Version     :	1.0
 Copyright   : 	Your copyright notice
 Description :	based in https://www.codeproject.com/Articles/5283245/Matrix-Library-in-C
 	 	 	 	Incorporates 	- error return
 	 	 	 	 	 	 	 	- Checks in Operations
 	 	 	 	 	 	 	 	- Matrixes zero indexed (as oposed to the original)
 	 	 	 	 	 	 	 	- Recursive algorithm to multiply matrices of any order
 	 	 	 	 	 	 	 	- Reduced Row Echelon form
 	 	 	 	 	 	 	 	- Inverse of any triangular Matrix (original only upper triangular)

To do		: 	Pthread Matix multiplication
				QR descomposition
				EigenValues
				EigenVectors

 ============================================================================
 */

#ifndef MATRIX_H_
#define MATRIX_H_
	//struct definition

	typedef struct Matrix_t{
		unsigned int row;
		unsigned int col;
		double* numbers;
	} Matrix_t;

/*
	typedef struct matrix_list_t{
		struct matrix_list_t* Next;
		Matrix_t* Matrix;
	} Matrix_List_t;
*/

	//Constructor
	Matrix_t* Matrix_new(unsigned int row, unsigned int col, double val); /* Allocates new matrix of row x col dimensions, filled with val */
	Matrix_t* Matrix_identity(unsigned int n); /* Allocates new Identity matrix of dimension n x n*/
	Matrix_t* Matrix_zeroes(unsigned int row, unsigned int col); /* Allocates new matrix of dimensions row * col, filled with zeroes */
	Matrix_t* Matrix_ones(unsigned int row, unsigned int col); /* Allocates new matrix of dimensions row * col, filled with ones */
	Matrix_t* Matrix_random(unsigned int row, unsigned int col, double from, double to); /* Allocates new matrix of dimensions row * col, filled with values between lower and upper inclusives */
	Matrix_t* Matrix_random_int(unsigned int row, unsigned int col, int lower, int upper); /* Allocates new matrix of dimensions row * col, filled with values between lower and upper inclusives */
	//destructor
	void Matrix_free(Matrix_t* In); /* free allocation of matrix M */
	//Getter
	int Matrix_get(Matrix_t* In, unsigned int row, unsigned int col, double* out); /* Get value of element in row row and column col, cero indexed*/
	//Setter
	int Matrix_set(Matrix_t* In, unsigned int row, unsigned int col, double val); /* Set value of element in row row and column col, cero indexed*/
	//operations
	int Matrix_scalarMultiplication(Matrix_t* In, double factor, Matrix_t* Out); /* Performs scalar multiplication, Out = Factor x In */
	int Matrix_addition(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /* Performs Matrix Addition, Out = In_1 + In_2 */
	int Matrix_subtraction(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /* Performs Matrix subtraction, Out = In_1 - In_2 */
	int Matrix_multiplication(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /*Performs Matrix multiplication, Out = In_1 * In_2 */
	int Matrix_multiplication_recursive(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /*Recursive algorithm for multiplying matrixes*/

	// Utilities
	void Matrix_display(Matrix_t* In); /* Prints the Matrix M in the console */

	Matrix_t* Matrix_copy(Matrix_t* In); /* generates a new matrix that is a copy of In */
	int Matrix_copyTo(Matrix_t* In, Matrix_t* Out); /* Copies values of In to Out */

	int Matrix_submatrix(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right, Matrix_t* Out); /* Obtains, in Out, the subMatrix of In between Columns Left and Right, and between rows Upper and Lower */

	int Matrix_removeRow(Matrix_t* In, unsigned int row, Matrix_t* Out); /* Removes Row row of the matrix */
	int Matrix_removeColumn(Matrix_t* In, unsigned int col, Matrix_t* Out); /* Removes column col of the matrix */
	int Matrix_transpose(Matrix_t* In, Matrix_t* Out); /* Transposes In Matrix */
	int Matrix_determinant(Matrix_t* In, double* out); /* Obtains the determinant of the matrix In*/
	int Matrix_trace(Matrix_t* In, double* out); /* Obtains the trace of the matrix */

	int Matrix_adjoint(Matrix_t* In, Matrix_t* Out); /* Calculates the adjoint matrix of In */
	int Matrix_inverse(Matrix_t* In, Matrix_t* Out); /* Inverts the Matrix In */
	int Matrix_inverse_triangular(Matrix_t* In, Matrix_t* Out); /* Inverts the triangular matrix In */

	Matrix_t* Matrix_rowEchelon(Matrix_t* In); /* Echelon form of the Matrix M */
	Matrix_t* Matrix_reducedRowEchelon(Matrix_t* In); /* Reduced Echelon form of the Matrix M */


	int Matrix_hConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /* Concatenates the matrix In_2 to In_1 horizontally */
	int Matrix_vConcat(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /* Concatenates the matrix In_2 to In_1 vertically */
	double Matrix_norm(Matrix_t* In); /* Calculates the norm of the Matrix In */
	Matrix_t* Matrix_kernel(Matrix_t* In); /* gets the kernel/null Space of the matrix In, each column of the output is a null vector */
	Matrix_t* Matrix_nullSpace(Matrix_t* In); /* gets the kernel/null Space of the matrix In, each column of the output is a null vector */

	int Matrix_compare(Matrix_t* In_1, Matrix_t* In_2); /*returns 0 if norm(In1) == norm(In2), <0 if norm(In1) < norm(In2), >0 if norm(In1) > norm(In2) */

	int Matrix_descomposition_LU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_U); /* Calculates the LU descomposition of the Matrix */
	int Matrix_descomposition_LDU(Matrix_t* In, Matrix_t* Out_L, Matrix_t* Out_D, Matrix_t* Out_U); /* Calculates the LDU descomposition of the Matrix */
	int Matrix_descomposition_QR(Matrix_t* In, Matrix_t* Out_Q, Matrix_t* Out_R); /* TODO: Calculates the QR descomposition of the Matrix */

	int Matrix_setMatrixtoZero(Matrix_t* In); /* Sets all the values of the Matrix to zero.*/

	//MUST CHECK MANUALLY DIMENSIONS, TRY NOT TO USE
	Matrix_t* Matrix_submatrix_2(Matrix_t* In, unsigned int upper, unsigned int lower, unsigned int left, unsigned int right); /* generates a new subMatrix of In between Columns Left and Right, and between rows Upper and Lower */
	Matrix_t* Matrix_hConcat_2(Matrix_t* In_1, Matrix_t* In_2); /* Concatenates the matrix In_2 to In_1 horizontally */
	Matrix_t* Matrix_vConcat_2(Matrix_t* In_1, Matrix_t* In_2); /* Concatenates the matrix In_2 to In_1 vertically */
	Matrix_t* Matrix_removeRow_2(Matrix_t* In, unsigned int row); /* Removes Row row of the matrix */
	Matrix_t* Matrix_removeColumn_2(Matrix_t* In, unsigned int col); /* Removes column col of the matrix */
	Matrix_t* Matrix_transpose_2(Matrix_t* In); /* Transposes In Matrix */

	//forked operations
#ifdef __PARALLEL__
#include <pthread.h>
	int Matrix_multiplication_forked(Matrix_t* In_1, Matrix_t* In_2, Matrix_t* Out); /*TODO: Forked Method*/
#endif

	//Error Codes
#define MATRIX_OK				0b000000000000000000000000

#define MATRIX_DIMENSIONS		0b000000000000000000000001
#define MATRIX_INDEX			0b000000000000000000000010
#define MATRIX_SQUARE       	0b000000000000000000000100
#define MATRIX_NON_INVERTIBLE   0b000000000000000000001000
#define MATRIX_TRIANGULAR_UPPER	0b000000000000000000010000
#define MATRIX_TRIANGULAR_LOWER	0b000000000000000000100000
#define MATRIX_TRIANGULAR		0b000000000000000001000000

#define MATRIX_DIFFERENT	 	0b001000000000000000000000
#define MATRIX_SAME				0b010000000000000000000000
#define MATRIX_ERROR			0b100000000000000000000000

#ifndef FLT_EPSILON
// This is the difference between 1 and the smallest floating point number of type float that is greater than 1.
// https://frama-c.com/2013/05/09/Definition-of-FLT_EPSILON.html
#define FLT_EPSILON 1.19209290E-07F
#endif //MATRIX_H_

#endif /* SRC_UTILS_MATRIX_H_ */
