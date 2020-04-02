#include <iostream>
#include <math.h>
#include <ctime>
#include <string>
#include "Matrix.h"
#include "Matrix.cpp"	// mandatory if using TEMPLATING
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <fstream>

using namespace std;

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
				auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // call SPD constructor
				auto* A = new CSRMatrix<double>(*A_dense); // convert to a sparse matrix
				auto* x_guess = new Matrix<double>(rows, 1, true); // making a guess x-vector
				auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
				auto* b = new Matrix<double>(rows, 1, true);

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

				// freeing-up memory 
				delete A_dense;
				delete A;
				delete x_exact;
				delete x_guess;
				delete b;


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
				A->gaussian_elimination_solve(*A, *b, *x_guess);
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

		else if (solver == "cholesky_solve")

		{
			fstream outfile("cholesky_dense.txt", ios::out); //output file stream

			for (int k = 10; k < 700; k = k + 10) // loop over increasing array sizes
			{

				int rows = k; // increasing rows and cols with each loop 
				int cols = k; 
				double tol = 10e-10;
				int iter = 10000;
				auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // make linear system
				auto* x_guess = new Matrix<double>(rows, 1, true); 
				auto* x_exact = new Matrix<double>(rows, 1, true); 
				auto* b = new Matrix<double>(rows, 1, true);

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

				// free-up memory 
				delete A;
				delete x_exact;
				delete x_guess;
				delete b;
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
				auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // making linear system
				auto* x_guess = new Matrix<double>(rows, 1, true);
				auto* x_exact = new Matrix<double>(rows, 1, true); 
				auto* b = new Matrix<double>(rows, 1, true);

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

				// free-up memory 
				delete A;
				delete x_exact;
				delete x_guess;
				delete b;

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
				auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // making linear system
				auto* x_guess = new Matrix<double>(rows, 1, true); 
				auto* x_exact = new Matrix<double>(rows, 1, true); 
				auto* b = new Matrix<double>(rows, 1, true);

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

				// free-up memory 
				delete A;
				delete x_exact;
				delete x_guess;
				delete b;
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
				string keyword = "noprint";


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
				A->conjgrad(x_guess, b, tol, iter, keyword); 
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


// tests each solver for 
void test_CSR_solvers(string solver)
{
	// TESTS for a known exact x vector in system

	// 1. Define an A matrix and x vector
	// 2. Calculate b vector: A*x = b
	// 3. Input A and b into user-specified solver 
	// 4. calculate RMS error between 

	int rows = 100;
	// 1. 

	if (solver == "conjugate_gradient") {
		auto* Adense = new Matrix<double>(rows, 500, 490, 10, 1, "fixed"); // CALLING A PSEUDO-RANDOM SPD MATRIX CONSTRUCTOR
		auto* A = new CSRMatrix<double>(*Adense); // WE CONVERT OUR KNOWN DENSE A MATRIX TO SPARSE
		auto* x_exact = new double[rows];
		auto* x_guess = new double[rows];
		auto* b = new double[rows];
		double tol = 10e-6;
		int maxiter = 1000;
		double error = 0;

		A->printAsDense();

		cout << "------------------------------------------" << endl;
		cout << "Testing sparse solver: " << solver << endl;
		cout << "------------------------------------------" << endl;
		cout << "Input linear system is: " << endl << endl;
		cout << "1. Our CSR-formatted matrix A is: " << endl << endl;
		A->printMatrix();
		cout << endl;
		cout << "2. Our exact x vector is: [";

		for (int i = 0; i < rows; i++)
		{

			//x_exact[i] = rand() % 50 + 2;; // pseudo-random exact x vector
			x_exact[i] = (i * 2) + 2; // pseudo-random exact x vector
			x_guess[i] = 0; // setting initial guess to zero
			cerr << x_exact[i] << " ";

		}
		cerr << "]" << endl;/*
		cerr << "3. Our initial x guess is: ";
		for (int i = 0; i < rows; i++)
		{
			cerr << x_guess[i] << " ";
		}*/

		// 2. calculate b
		A->matVecMult(x_exact, b); /// A * x_exact = b, we output b vector here.
		cerr << endl << "3. Our b vector: [";
		for (int i = 0; i < rows; i++)
		{

			cout << b[i] << " ";

		}
		cout << "]" << endl;




		// 3. CALL SOLVER: for now just conjugate gradient
		A->conjgrad(x_guess, b, tol, maxiter, "noprint"); // now output our x guess

		// 4. calculate RMS error
		cerr << endl << "3. Our output x vector guess using " << solver << " method is: [";
		for (int i = 0; i < A->rows; i++)
		{
			error += pow((x_guess[i] - x_exact[i]), 2);  // square of difference between values 
			cout << x_guess[i] << " ";
		}
		cerr << "]" << endl;

		error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
		cout << "with an RMS error between exact and solver x vectors of: " << error << endl;



		delete Adense;
		delete A;
		delete[] x_exact;
		delete[] x_guess;
		delete[] b;

	}
	else if (solver == "gauss_seidel") {
		double error = 0;
		auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
		auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
		auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
		auto* b_dense = new Matrix<double>(rows, 1, true); // b vector
		
		cout << "------------------------------------------" << endl;
		cout << "Testing sparse solver: " << solver << endl;
		cout << "------------------------------------------" << endl;
		cout << "Input linear system is: " << endl << endl;
		cout << "1. Our exact x vector is: [";
		for (int i = 0; i < rows; i++)
		{
			//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
			x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
			x_guess->values[i] = 0; // setting initial guess to zero
			cout << x_exact->values[i] << " ";
		}
		cerr << "]" << endl << endl;

		// calculate b
		A_dense->matMatMult(*x_exact, *b_dense); /// A * x_exact = b, we output b vector here.

		// Now convert to sparse matrices so we can test our Gauss-seidel solver 
		auto* A = new CSRMatrix<double>(*A_dense); // WE CONVERT OUR KNOWN DENSE A MATRIX TO SPARS

		cerr << "2. Our b vector: [";
		for (int i = 0; i < rows; i++)
		{
			cout << b_dense->values[i] << " ";

		}
		cout << "]" << endl << endl;

		cout << "3. Our CSR-formatted matrix A is: " << endl << endl;
		A->printAsDense();

		// solve using Gauss-Seidel function
		A->gauss_seidel(*A, *b_dense, *x_guess);

		cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
		for (int i = 0; i < A->rows; i++)
		{
			error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
			cout << x_guess->values[i] << " ";
		}
		cerr << "]" << endl;

		error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
		cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

		delete A_dense;
		delete A;
		delete x_exact;
		delete x_guess;
		delete b_dense;
	}
	else
	{
		cerr << endl << "No valid solver specified." << endl;
	}
}

void test_dense_solvers(string solver) {


	// TESTS for a known exact x vector in system

	// 1. Define an A matrix and exact x vector
	// 2. Calculate b vector: A*x = b
	// 3. Input A and b into user-specified solver
	// 4. calculate RMS error between and output x guess



	int rows = 10;

	if (solver == "conjugate_gradient") {
		auto* A = new Matrix<double>(rows, 500, 490, 10, 1, "fixed"); // CALLING A PSEUDO-RANDOM SPD MATRIX CONSTRUCTOR
		auto* x_exact = new double[rows];
		auto* x_guess = new double[rows];
		auto* b = new double[rows];
		double tol = 10e-6;
		int maxiter = 1000;
		double error = 0;

		cout << "------------------------------------------" << endl;
		cout << "Testing dense solver: " << solver << endl; 
		cout << "------------------------------------------" << endl;
		cout << "Input linear sytem: " << endl << endl;
		cout << "1. Our matrix A is: " << endl << endl;
		A->printMatrix();
		cout << endl;
		cout << "2. Our exact x vector is: [";

		for (int i = 0; i < rows; i++)
		{

			//x_exact[i] = rand() % 50 + 2;; // pseudo-random exact x vector
			x_exact[i] = (i * 2) + 2; // pseudo-random exact x vector
			x_guess[i] = 0; // setting initial guess to zero
			cerr << x_exact[i] << " ";

		}
		cerr << "]" << endl;


		// 2. calculate b

		A->matVecMult(x_exact, b); /// A * x_exact = b, we output b vector here.
		cerr << endl << "3. Our b vector: [";
		for (int i = 0; i < rows; i++)
		{

			cout << b[i] << " ";

		}
		cout << "]" << endl;




		// 3. CALL SOLVER: for now just conjugate gradient
		A->conjgrad(x_guess, b, tol, maxiter, "noprint"); // now output our x guess

		// 4. calculate RMS error
		cerr << endl << "3. Our output x vector guess using " << solver << " method is: [";
		for (int i = 0; i < A->rows; i++)
		{
			error += pow((x_guess[i] - x_exact[i]), 2);  // square of difference between values 
			cout << x_guess[i] << " ";
		}
		cerr << "]" << endl;

		error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
		cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

	}

	else if (solver == "gauss_seidel") {
		double error = 0;
		auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
		auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
		auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
		auto* b = new Matrix<double>(rows, 1, true); // b vector

		cout << "------------------------------------------" << endl;
		cout << "Testing dense solver: " << solver << endl;
		cout << "------------------------------------------" << endl;
		cout << "Input linear system is: " << endl << endl;
		cout << "1. Our exact x vector is: [";
		for (int i = 0; i < rows; i++)
		{
			//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
			x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
			x_guess->values[i] = 0; // setting initial guess to zero
			cout << x_exact->values[i] << " ";
		}
		cerr << "]" << endl << endl;

		// calculate b
		A->matMatMult(*x_exact, *b); /// A * x_exact = b, we output b vector here.

		cerr << "2. Our b vector: [";
		for (int i = 0; i < rows; i++)
		{
			cout << b->values[i] << " ";

		}
		cout << "]" << endl << endl;

		cout << "3. Our matrix A is: " << endl << endl;
		A->printMatrix();



		// solve using Gauss-Seidel function
		A->gauss_seidel(*A, *b, *x_guess);

		cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
		for (int i = 0; i < A->rows; i++)
		{
			error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
			cout << x_guess->values[i] << " ";
		}
		cerr << "]" << endl;

		error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
		cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

		delete A;
		delete x_exact;
		delete x_guess;
		delete b;
	}

	else if (solver == "cholesky") {
	double error = 0;
	auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
	auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* b = new Matrix<double>(rows, 1, true); // b vector

	cout << "------------------------------------------" << endl;
	cout << "Testing dense solver: " << solver << endl;
	cout << "------------------------------------------" << endl;
	cout << "Input linear system is: " << endl << endl;
	cout << "1. Our exact x vector is: [";
	for (int i = 0; i < rows; i++)
	{
		//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
		x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
		x_guess->values[i] = 0; // setting initial guess to zero
		cout << x_exact->values[i] << " ";
	}
	cerr << "]" << endl << endl;

	// calculate b
	A->matMatMult(*x_exact, *b); /// A * x_exact = b, we output b vector here.

	cerr << "2. Our b vector: [";
	for (int i = 0; i < rows; i++)
	{
		cout << b->values[i] << " ";

	}
	cout << "]" << endl << endl;

	cout << "3. Our matrix A is: " << endl << endl;
	A->printMatrix();


	A->cholesky_solve(*A, *b, *x_guess);

	cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
	for (int i = 0; i < A->rows; i++)
	{
		error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
		cout << x_guess->values[i] << " ";
	}
	cerr << "]" << endl;

	error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
	cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

	delete A;
	delete x_exact;
	delete x_guess;
	delete b;
	cout << "Finished cholesky\n";
	}

	else if (solver == "jacobi") {
	double error = 0;
	auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
	auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* b = new Matrix<double>(rows, 1, true);		// b vector

	cout << "------------------------------------------" << endl;
	cout << "Testing dense solver: " << solver << endl;
	cout << "------------------------------------------" << endl;
	cout << "Input linear system is: " << endl << endl;
	cout << "1. Our exact x vector is: [";
	for (int i = 0; i < rows; i++)
	{
		//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
		x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
		x_guess->values[i] = 0; // setting initial guess to zero
		cout << x_exact->values[i] << " ";
	}
	cerr << "]" << endl << endl;

	// calculate b
	A->matMatMult(*x_exact, *b); /// A * x_exact = b, we output b vector here.

	cerr << "2. Our b vector: [";
	for (int i = 0; i < rows; i++)
	{
		cout << b->values[i] << " ";

	}
	cout << "]" << endl << endl;

	cout << "3. Our matrix A is: " << endl << endl;
	A->printMatrix();


	A->jacobi(*A, *b, *x_guess);

	cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
	for (int i = 0; i < A->rows; i++)
	{
		error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
		cout << x_guess->values[i] << " ";
	}
	cerr << "]," << endl;

	error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
	cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

	delete A;
	delete x_exact;
	delete x_guess;
	delete b;
	}

	else if (solver == "gaussian_elimination") {
	double error = 0;
	auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
	auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* b = new Matrix<double>(rows, 1, true); // b vector

	cout << "------------------------------------------" << endl;
	cout << "Testing dense solver: " << solver << endl;
	cout << "------------------------------------------" << endl;
	cout << "Input linear system is: " << endl << endl;
	cout << "1. Our exact x vector is: [";
	for (int i = 0; i < rows; i++)
	{
		//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
		x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
		x_guess->values[i] = 0; // setting initial guess to zero
		cout << x_exact->values[i] << " ";
	}
	cerr << "]" << endl << endl;

	// calculate b
	A->matMatMult(*x_exact, *b); /// A * x_exact = b, we output b vector here.

	cerr << "2. Our b vector: [";
	for (int i = 0; i < rows; i++)
	{
		cout << b->values[i] << " ";

	}
	cout << "]" << endl << endl;

	cout << "3. Our matrix A is: " << endl << endl;
	A->printMatrix();


	A->gaussian_elimination_solve(*A, *b, *x_guess);

	cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
	for (int i = 0; i < A->rows; i++)
	{
		error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
		cout << x_guess->values[i] << " ";
	}
	cerr << "]" << endl;

	error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
	cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

	delete A;
	delete x_exact;
	delete x_guess;
	delete b;
	}

	else if (solver == "LU_solve") {
	double error = 0;

	auto* A = new Matrix<double>(rows, 500, 495, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
	auto* x_guess = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* x_exact = new Matrix<double>(rows, 1, true); // making a known test x-vector
	auto* b = new Matrix<double>(rows, 1, true); // b vector

	cout << "------------------------------------------" << endl;
	cout << "Testing dense solver: " << solver << endl;
	cout << "------------------------------------------" << endl;
	cout << "Input linear system is: " << endl << endl;
	cout << "1. Our exact x vector is: [";
	for (int i = 0; i < rows; i++)
	{
		//x_exact->values[i] = rand() % 50 + 2;; // pseudo-random exact x vector
		x_exact->values[i] = (i * 2) + 2; // pseudo-random exact x vector
		x_guess->values[i] = 0; // setting initial guess to zero
		cout << x_exact->values[i] << " ";
	}
	cerr << "]" << endl << endl;

	// calculate b
	A->matMatMult(*x_exact, *b); /// A * x_exact = b, we output b vector here.

	cerr << "2. Our b vector: [";
	for (int i = 0; i < rows; i++)
	{
		cout << b->values[i] << " ";

	}
	cout << "]" << endl << endl;

	cout << "3. Our matrix A is: " << endl << endl;
	A->printMatrix();


	A->LU_solve(*A, *b, *x_guess);

	cerr << endl << "4. Our output x vector guess using " << solver << " method is:" << endl << "[";
	for (int i = 0; i < A->rows; i++)
	{
		error += pow((x_guess->values[i] - x_exact->values[i]), 2);  // square of difference between values 
		cout << x_guess->values[i] << " ";
	}
	cerr << "]" << endl;

	error = sqrt(error / A->rows); // divide sum of differences by number of values, then take the square root 
	cout << "with an RMS error between exact and solver x vectors of: " << error << endl;

	delete A;
	delete x_exact;
	delete x_guess;
	delete b;
	
	}

	else {
	// give error i
	cerr << endl << "Error - no valid solver specified. " << endl;
}


}

int main()
{

	int rows = 4;
	int cols = 4;

	double A_data[] = { 1., 0., 3., 7., 2., 1., 0., 4., 5., 4., 1., -2., 4., 1., 6., 2. };
	double* A_vals = A_data;
	auto* dense_mat = new Matrix<double>(rows, rows, A_vals);

	double b_data[] = { 1., 2., -3., 2. };
	double* b_vals = b_data;
	auto* b = new Matrix<double>(rows, 1, b_vals);

	auto* output = new Matrix<double>(rows, 1, true);

	std::cout << "Doing matrix multiplication between: \n";
	std::cout << "A: \n";
	dense_mat->printMatrix();
	std::cout << "and B: \n";
	b->printMatrix();
	// EXAMPLE FOR HOW TO RUN CODE NORMALLY (WITHOUT OVERLOADING)
	dense_mat->matMatMult(*b, *output);		// <==> *(dense_mat).matMatmult(...)
	std::cout << "Result: \n";
	output->printMatrix();

	std::cout << "\nOverwriting the * operator to do matrix multiplication between: \n";
	std::cout << "A: \n";
	dense_mat->printMatrix();
	std::cout << "and B: \n";
	b->printMatrix();
	// EXAMPLE FOR HOW TO RUN CODE USING OVERLOADING
	// You are now responsible for this new chunk of data - DOCUMENT THIS!!!
	output = *dense_mat * *b;
	std::cout << "Result: \n";
	output->printMatrix();

	std::cout << "Original Dense Matrix" << endl;
	dense_mat->printMatrix();
	output = dense_mat->transpose();
	std::cout << "Transposed Dense Matrix" << endl;
	output->printMatrix();
	delete output;
	output = nullptr;

	std::cout << "Print identity: \n";
	auto* identity = new Matrix<double>(5, "eye");
	identity->printMatrix();
	delete identity;

	std::cout << "Print zeros: \n";
	auto* zeros = new Matrix<double>(5, "zeros");
	zeros->printMatrix();
	delete zeros;

	std::cout << "Print diagonal: \n";
	auto* diag = new Matrix<double>(b_vals, 4, "diag");
	diag->printMatrix();
	delete diag;

	output = new Matrix<double>(rows, 1, true);

	std::cout << "\nSolving Ax = b system for x: \n";
	//DOCUMENT THAT THE OUTPUT NEEDS TO BE DELETED
	dense_mat->gaussian_elimination_solve(*dense_mat, *b, *output);
	std::cout << "GE solution for x: \n";
	output->printMatrix();
	
	//DOCUMENT THAT THE OUTPUT NEEDS TO BE DELETED
	dense_mat->LU_solve(*dense_mat, *b, *output);

	std::cout << "LU solution for x: \n";
	output->printMatrix();

	cout << "\nTrying to delete\n";
	delete dense_mat;
	delete b;
	delete output;

	//// TESTING SPARSE MATRICES
	// Define some sparse class matrices
	int rows_sparse = 4;
	int cols_sparse = 4;
	int nnzs = 6;

	double val_data[] = { 1., 2., 5., 3., 1., 3. };
	double* val_input = val_data;
	int col_data[] = { 2, 3, 1, 2, 0, 1 };
	int* col_input = col_data;
	int row_data[] = { 0, 2, 3, 4, 6 };
	int* row_input = row_data;

	// Define sparse class
	auto* sparse_mat = new CSRMatrix<double>(rows_sparse, cols_sparse, nnzs, val_input, row_input, col_input);
	std::cout << "Print Sparse Values: \n";
	sparse_mat->printValues();
	std::cout << endl;

	std::cout << "Print Sparse Matrix:\n";
	sparse_mat->printMatrix();
	std::cout << endl;

	std::cout << "Print Sparse Matrix as dense:\n";
	sparse_mat->printAsDense();
	std::cout << endl;

	auto* output_to_dense = new Matrix<double>();
	sparse_mat->sparse2dense(*output_to_dense);
	std::cout << "Print converted Matrix:\n";
	output_to_dense->printMatrix();
	std::cout << endl;
	delete output_to_dense;

	double* output_matVec = new double[rows];
	double input_vals_matVec[] = { 1., 3., 7., 2., };
	double* input_matVec = input_vals_matVec;
	std::cout << "Doing matVecMult\n";
	sparse_mat->matVecMult(input_matVec, output_matVec);
	std::cout << "Printing sparse_mat_poly: \n";
	for (int i = 0; i < rows; i++) {
		std::cout << output_matVec[i] << '\n';
	}
	std::cout << endl;
	delete[] output_matVec;

	double A_val_data[] = { 5., 8., 3., 6. };
	double* A_val_input = A_val_data;
	int A_col_data[] = { 0, 1, 2, 1 };
	int* A_col_input = A_col_data;
	int A_row_data[] = { 0, 0, 2, 3, 4 };
	int* A_row_input = A_row_data;

	// Define sparse class
	CSRMatrix<double>* A_mat = new CSRMatrix<double>(4, 5, 5, A_val_input, A_row_input, A_col_input);
	A_mat->printAsDense();

	double B_val_data[] = { 10., 20., 30., 40., 50., 60., 70., 80. };
	double* B_val_input = B_val_data;
	int B_col_data[] = { 0, 1, 1, 3, 2, 3, 4, 5 };
	int* B_col_input = B_col_data;
	int B_row_data[] = { 0, 2, 4, 7, 8, 8 };
	int* B_row_input = B_row_data;


	// Define sparse class
	CSRMatrix<double>* B_mat = new CSRMatrix<double>(5, 6, 8, B_val_input, B_row_input, B_col_input);
	//B_mat->printAsDense();
	CSRMatrix<double>* out_mat = new CSRMatrix<double>();
	cout << "Doing mat-mat mult\n";
	//A_mat->matMatMult(*B_mat, *out_mat);
	cout << "...using the overwritten * operator" << endl;
	out_mat = *A_mat * *B_mat;

	cout << "Printing output matrix: \n";
	out_mat->printMatrix();
	cout << endl;

	cout << "Printing output matrix as dense: \n";
	out_mat->printAsDense();
	cout << endl;

	delete sparse_mat;
	delete A_mat;
	delete B_mat;
	delete out_mat;

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


	solver_runtime("dense", "gauss_seidel");
	cerr << "Gauss-S dense finished" << endl;
	
	solver_runtime("dense", "jacobi");
	cerr << "jacobi dense finished" << endl;


	//////////////////////////////////////
	////////// TESTING SOLVERS ///////////
	//////////////////////////////////////


	// 1. Sparse matrix Solvers
	test_CSR_solvers("conjugate_gradient");
	test_CSR_solvers("gauss_seidel");

	// 2. Dense matrix Solvers
	test_dense_solvers("conjugate_gradient");
	test_dense_solvers("gauss_seidel");
	test_dense_solvers("cholesky");

	test_dense_solvers("jacobi");
	test_dense_solvers("gaussian_elimination");
	test_dense_solvers("LU_solve");

	return 0;
}
