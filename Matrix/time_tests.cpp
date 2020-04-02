#include <iostream>
#include <math.h>
#include <ctime>
#include <string>
#include "Matrix.h"
#include "Matrix.cpp"	// mandatory if using TEMPLATING
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <fstream>

// function to find relationship between runtime and matrix size for each solver
// inputs: strings - (i) dense or sparse, (ii) the desired solver
void  solver_runtime(string mattype, string solver)
{ /* Calculates runtime for variety of sparse solvers with progressively increasing system sizes.
  Outputs results to a txt file for plotting.
  */

  // Timing tests for our sparse solvers
	if (mattype == "sparse") // runs sparse matrix solvers
	{

		if (solver == "conjgrad")
		{
			fstream outfile("conjgrad_sparse.txt", ios::out); //output file stream
			for (int k = 1000; k <= 16500; k = k + 250) // starts off 10x10 
			{

				int rows = k; // increasing rows and cols with each loop 
				int cols = k;
				double tol = 10e-10;
				int iter = 100000;


				auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // call SPD constructor to make a matrix
				auto* A = new CSRMatrix<double>(*A_dense); // convert dense matrix to sparse
				auto* x_exact = new double[cols]; // making a known test x-vector
				auto* x_guess = new double[cols]; // making an input x-vector
				auto* b = new double[rows]; // output vector


				for (int i = 0; i < cols; i++)
				{
					x_exact[i] = (i + 2) * 2;  // set exact x values
					x_guess[i] = 0; // set initial guess array to zeros

				}


				// Timing 
				clock_t start = clock();
				A->conjgrad(x_guess, b, tol, iter); // calculating our estimated x
				clock_t end = clock();
				float total = (float)(end - start) / (float)(CLOCKS_PER_SEC) * 1000.0; // convert time to ms
				outfile << rows << "," << (float)(end - start) / (float)(CLOCKS_PER_SEC) * 1000.0 << endl;
				//cout << rows << "," << total << endl;

				// freeing-up memory 
				delete A_dense;
				delete A;
				delete[] x_exact;
				delete[] x_guess;
				delete[] b;


			}
			outfile.close();

		}


		else if (solver == "gauss_seidel")
		{
			fstream outfile("gauss_seidel_sparse.txt", ios::out); //output file stream

			for (int k = 1000; k <= 16500; k = k + 250) // loop over increasing matrix sizes
			{
				int rows = k; // increasing rows and cols with each loop 
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;
				unique_ptr<Matrix<double>> A_dense(new Matrix<double>(rows, 500, 495, 2, 1, "fixed")); // call SPD constructor
				unique_ptr<CSRMatrix<double>> A(new CSRMatrix<double>(*A_dense)); // convert to a sparse matrix
				unique_ptr<Matrix<double>> x_guess(new Matrix<double>(rows, 1, true)); // making a guess x-vector
				unique_ptr<Matrix<double>> x_exact(new Matrix<double>(rows, 1, true)); // making a known test x-vector
				unique_ptr<Matrix<double>> b(new Matrix<double>(rows, 1, true));

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2; // set x vector values
					x_guess->values[i] = 0;

				}

				// timing
				clock_t start = clock(); // timeing calculation time 
				A->gauss_seidel(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl; // << "," << res << endl;

			}


			outfile.close();

		}

	}

	// Timing tests for our dense solvers
	else {  // runs dense matrix solvers

		if (solver == "gaussian_elimination_solve") {

			fstream outfile("gauss_elim_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // loop over increasing array sizes
			{

				int rows = k;
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;

				unique_ptr<Matrix<double>> A(new Matrix<double>(rows, 500, 495, 2, 1, "fixed")); // call SPD constructor
				unique_ptr<Matrix<double>> x_guess(new Matrix<double>(rows, 1, true)); // making a guess x-vector
				unique_ptr<Matrix<double>> x_exact(new Matrix<double>(rows, 1, true)); // making a known test x-vector
				unique_ptr<Matrix<double>> b(new Matrix<double>(rows, 1, true));

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2;
					x_guess->values[i] = 0;

				}

				A->matMatMult(*x_exact, *b); // outputting our exact b vector from a known x and A

				// Timing
				clock_t start = clock();
				A->gaussian_elimination_solve(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

			}

			outfile.close();

		}

		else if (solver == "cholesky_solve")

		{
			fstream outfile("cholesky_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // loop over increasing array sizes
			{

				int rows = k; // increasing rows and cols with each loop 
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;
				unique_ptr<Matrix<double>> A(new Matrix<double>(rows, 500, 495, 2, 1, "fixed")); // call SPD constructor
				unique_ptr<Matrix<double>> x_guess(new Matrix<double>(rows, 1, true)); // making a guess x-vector
				unique_ptr<Matrix<double>> x_exact(new Matrix<double>(rows, 1, true)); // making a known test x-vector
				unique_ptr<Matrix<double>> b(new Matrix<double>(rows, 1, true));

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2; // setting x vector values 
					x_guess->values[i] = 0;

				}


				A->matMatMult(*x_exact, *b); // outputting our exac5 b vector from a known x and A

				// timing
				clock_t start = clock();
				A->cholesky_solve(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

			}




		}

		else if (solver == "gauss_seidel")

		{
			fstream outfile("gauss_seidel_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // loop over increasing array sizes
			{

				int rows = k; // increasing rows and cols with each loop
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;
				unique_ptr<Matrix<double>> A(new Matrix<double>(rows, 500, 495, 2, 1, "fixed")); // call SPD constructor
				unique_ptr<Matrix<double>> x_guess(new Matrix<double>(rows, 1, true)); // making a guess x-vector
				unique_ptr<Matrix<double>> x_exact(new Matrix<double>(rows, 1, true)); // making a known test x-vector
				unique_ptr<Matrix<double>> b(new Matrix<double>(rows, 1, true));

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2; // setting x vector values
					x_guess->values[i] = 0;

				}


				A->matMatMult(*x_exact, *b); // outputting our exac5 b vector from a known x and A

				// timing 
				clock_t start = clock();
				A->gauss_seidel(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl; // << "," << res << endl;

			}




		}

		else if (solver == "jacobi")

		{
			fstream outfile("jacobi_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // increasing array sizes 
			{

				int rows = k;
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;
				unique_ptr<Matrix<double>> A(new Matrix<double>(rows, 500, 495, 2, 1, "fixed")); // call SPD constructor
				unique_ptr<Matrix<double>> x_guess(new Matrix<double>(rows, 1, true)); // making a guess x-vector
				unique_ptr<Matrix<double>> x_exact(new Matrix<double>(rows, 1, true)); // making a known test x-vector
				unique_ptr<Matrix<double>> b(new Matrix<double>(rows, 1, true));

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2; // setting x vector values
					x_guess->values[i] = 0;

				}

				A->matMatMult(*x_exact, *b); // outputting our exac5 b vector from a known x and A

				// timing 
				clock_t start = clock();
				A->jacobi(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl; // << "," << res << endl;

			}

		}

		else if (solver == "conjugate_gradient")

		{
			fstream outfile("conjugate_gradient_dense.txt", ios::out); //output file stream
			for (int k = 10; k <= 700; k = k + 10) // increasing array size with each loop
			{

				int rows = k;
				int cols = k;
				double tol = 10e-10;
				int iter = 100000;


				auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // making linear system
				auto* x_exact = new double[cols];
				auto* x_guess = new double[cols];
				auto* b = new double[rows];

				for (int i = 0; i < cols; i++)
				{
					x_exact[i] = (i + 2) * 2; // setting x vector values
					x_guess[i] = 0;

				}

				// timing 
				clock_t start = clock();
				A->conjgrad(x_guess, b, tol, iter);
				clock_t end = clock();
				float total = (float)(end - start) / (float)(CLOCKS_PER_SEC) * 1000.0;
				outfile << rows << "," << total << endl;
				total = 0;

				// free-up memory 
				delete A;
				delete[] x_exact;
				delete[] x_guess;
				delete[] b;


			}
			outfile.close();

		}

		else if (solver == "LU_solve") {

			fstream outfile("LU_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // loop over increasing array sizes
			{

				int rows = k;
				int cols = k;
				double tol = 10e-10;
				int iter = 10000;
				auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // make matrix system
				auto* x_guess = new Matrix<double>(rows, 1, true);
				auto* x_exact = new Matrix<double>(rows, 1, true);
				auto* b = new Matrix<double>(rows, 1, true);

				for (int i = 0; i < rows; i++)
				{
					x_exact->values[i] = (i + 1) * 2;
					x_guess->values[i] = 0;

				}

				A->matMatMult(*x_exact, *b); // outputting our exact b vector from a known x and A

				// Timing
				clock_t start = clock();
				A->LU_solve(*A, *b, *x_guess);
				clock_t end = clock();
				outfile << rows << "," << (double)(end - start) / (double)(CLOCKS_PER_SEC) * 1000.0 << endl;

				// free-up memory
				delete A;
				delete x_exact;
				delete x_guess;
				delete b;

			}

			outfile.close();

		}


	}

}

void time_tests()
{

	cout << "Timing .txt files written out to current directory\n";
	cout << endl;

	solver_runtime("dense", "gaussian_elimination_solve");
	cerr << "Gauss-elim dense finished" << endl;

	solver_runtime("dense", "cholesky_solve");
	cerr << "cholesk dense finished" << endl;

	solver_runtime("dense", "gauss_seidel");
	cerr << "Gauss-S dense finished" << endl;

	solver_runtime("dense", "jacobi");
	cerr << "jacobi dense finished" << endl;

	solver_runtime("dense", "conjugate_gradient");
	cerr << "conjgrad finished" << endl;

	solver_runtime("dense", "LU_solve");
	cerr << "LU dense finished" << endl;

	solver_runtime("sparse", "gauss_seidel");
	cerr << "Gauss-S sparse finished" << endl;

	solver_runtime("sparse", "conjgrad");
	cerr << "conjgrad sparse finished" << endl;

}