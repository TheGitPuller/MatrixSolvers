#include "Matrix.h"
#include "Matrix.cpp"	// mandatory if using TEMPLATING
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"

// tests each solver 
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
		A->conjgrad(x_guess, b, tol, maxiter); // now output our x guess

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
void test_CSR_solvers(string solver)
{
	// TESTS for a known exact x vector in system

	// 1. Define an A matrix and x vector
	// 2. Calculate b vector: A*x = b
	// 3. Input A and b into user-specified solver 
	// 4. calculate RMS error between 

	int rows = 50;
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
		A->conjgrad(x_guess, b, tol, maxiter); // now output our x guess

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

void tests()
{
	//////////////////////////////////////
	////////// TESTING SOLVERS ///////////
	//////////////////////////////////////

	// 1. Dense matrix Solvers
	test_dense_solvers("conjugate_gradient");

	test_dense_solvers("gauss_seidel");

	test_dense_solvers("cholesky");

	test_dense_solvers("jacobi");

	test_dense_solvers("gaussian_elimination");

	test_dense_solvers("LU_solve");

	// 2. Sparse matrix Solvers
	test_CSR_solvers("conjugate_gradient");

	test_CSR_solvers("gauss_seidel");
}