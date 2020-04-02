#include <iostream>
#include <math.h>
#include <ctime>
#include <string>
#include "Matrix.h"
#include "Matrix.cpp"	// mandatory if using TEMPLATING
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <cblas.h>
#include <stdio.h>
#include <fstream>
using namespace std;

// may need to use __restrict on some of the parameters
// may be good to compare BLAS/LAPACK for performance comparisons and potentially plot with gnuplot.
void  solver_runtime(string mattype, string solver)
{ /* Finds runtime for variety of sparse solvers*/

	// Contiguous memory is always faster. Where incy and incy are skipping
	// over array-items means that memory stored in cache is accessed sporadically,
	// thus we pay a performance penalty;

	//////////////////////////////////////////
	if (mattype == "dense") { // runs dense matrix solvers

		if (solver == "matMatMult") {
			//cout << "Running time tests on matMatMult\n";
			fstream outfile("matMatMult.txt", ios::out); //output file stream
			outfile.precision(10);
			for (int k = 24; k <= 5600; k *= 2) // starts off 10x10 
			{

				const int rows = k; // altering rows and cols with each loop 
				const int cols = k;

				auto* A = new Matrix<double>(rows, rows, rows, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
				auto* B = new Matrix<double>(rows, rows, rows, 2, 1, "fixed");
				auto* output = new Matrix<double>(rows, rows, true);


				auto* blas_A = new double[rows * cols];
				auto* blas_B = new double[rows * cols];
				auto* blas_output = new double[rows * cols];

				for (int i = 0; i < rows * cols; i++) {
					blas_A[i] = A->values[i];
					blas_B[i] = (B->values[i]);
					blas_output[i] = (output->values[i]);
				}

				cout << "Rows = " << k << endl;

				// calculate residuals 
				clock_t start = clock();	// timing calculation time 
				A->matMatMult(*B, *output); // outputting our exac5 b vector from a known x and A
				clock_t end = clock();

				clock_t start_blas = clock(); // timing calculation time 
				////cblas_dgemm(ORDER, TRANSPOSE, TRANSPOSE, int, int, int, double, const double*, int, const double*, int, double, double*, int) 
				////cblas_dgemm(order, transA, transB, M, N, K, alpha, A, LDA, B, LDB, beta, C, LDC)
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, rows, rows, rows, 1., blas_A, rows, blas_B, rows, 1., blas_output, rows);
				clock_t end_blas = clock();


				//cout << "Time spent (ms??): " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << "," << (double)(end_blas - start_blas) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

				delete A;
				delete B;
				delete output;
				delete[] blas_A;
				delete[] blas_B;
				delete[] blas_output;

			}

			outfile.close();

		}
		if (solver == "matVecMult") {
			cout << "Running time tests on matVecMult\n";
			fstream outfile("dgetrf.txt", ios::out); //output file stream
			outfile.precision(10);
			for (int k = 10; k <= 200; k *= 1.2) // starts off 10x10 - prev 2500
			{
				const int rows = k; // altering rows and cols with each loop 
				const int cols = k;

				auto* A = new Matrix<double>(rows, rows, rows, 2, 1, "fixed");	 // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
				auto* P = new Matrix<double>(*A);
				auto* L = new Matrix<double>(rows, "zeros");
				auto* U = new Matrix<double>(rows, "eye");

				auto* blas_A = new double[rows * cols];
				auto* IPIV = new int[rows];
				
				for (int i = 0; i < rows * cols; i++) {
					blas_A[i] = A->values[i];
				}
				for (int i = 0; i < rows; i++) {
					IPIV[i] = 0;
				}

				cout << "Rows = " << k << endl;


				// calculate residuals 
				clock_t start = clock();	// timing calculation time 
				//A->LU_Decomposition(Matrix<T> & iA, Matrix<T> & P_, Matrix<T> & L, Matrix<T> & U);	// outputting our exact b vector from a known x and A
				A->LU_Decomposition(*A, *P, *L, *U);
				clock_t end = clock();
				delete P;
				delete U;
				delete L;

				int INFO;
				clock_t start_blas = clock(); // timing calculation time 
				////cblas_dgemm(ORDER, TRANSPOSE, TRANSPOSE, int, int, int, double, const double*, int, const double*, int, double, double*, int) 
				////cblas_dgemm(order, transA, M, N, alpha, A, LDA, B, LDB, beta, C, LDC)
				// cblas_dgetrf(int M, int N, double A, int LDA, int IPIV, int INFO);
				cblas_dgetrf(400, 400, *blas_A, 1, IPIV, INFO);
				clock_t end_blas = clock();

				//cout << "Time spent (ms??): " << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << "," << (double)(end_blas - start_blas) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

				delete A;
				delete[] blas_A;
				delete[] IPIV;
			}

			outfile.close();

		}

	}

}



int main()
{

	solver_runtime("dense", "matMatMult");
	//solver_runtime("dense", "matVecMult");
	return 0;
}
