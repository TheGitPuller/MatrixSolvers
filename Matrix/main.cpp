#include "Matrix.h"
#include "Matrix.cpp"	
#include "CSRMatrix.h"
#include "CSRMatrix.cpp"
#include <chrono>
using namespace std;

void GUI(){
	string mat_type;
	string method;
	int rows;
	char choice = 'Y';
	cout << "Welcome to this linear systems solver" << endl;
	cout << "This GUI enables you to generate pseudo-random Symmetric Positive Definite (SPD)" << endl;
	cout << "matrix system and solve it using different methods" << endl;
	cout << "Press enter a key to continue" << endl<< endl;
	cin.get();
	while (choice == 'Y') {
		while (mat_type != "dense" && mat_type != "sparse") {
			cout << "Would you like to generate a dense or a sparse matrix?" << endl;
			cout << "Input 'dense' or 'sparse'" << endl << endl;
			cin >> mat_type;
			if (mat_type == "dense" || mat_type == "sparse") {
				cout << "An NxN matrix will be generated." << endl;
				cout << "Please enter a value for N: " << endl << endl;
				cin >> rows;

				while (cin.fail()) {
					cout << "Error, input is not an integer" << std::endl;
					cin.clear();
					cin.ignore(256, '\n');
					cin >> rows;
				}
				if (rows < 0) {
					cout << "Error, input is negative: " << rows << " ... converting to positive: " << abs(rows) << endl;
					rows = abs(rows);
				}
			}
			if (mat_type == "dense") {

				auto* A = new Matrix<double>(rows, 500, 300, 10, 1, "random");

				auto* X = new Matrix<double>(rows, 1, true);

				for (int j = 0; j < rows; j++) { X->values[j] = 0; }

				auto* B = new Matrix<double>(rows, 1, true);

				for (int j = 0; j < rows; j++) { B->values[j] = rand() % (20 - 10 + 1) + 10; }

				cout << "The system Ax = b will be solved, where: " << endl << endl;
				cout << "A: " << endl;

				A->printMatrix();

				cout << "b: " << "[ ";
				for (int j = 0; j < rows; j++) { cerr << B->values[j] << " "; }
				cout << "]" << endl;

				//B->printMatrix();

				while (method != "jac" && method != "GS" && method != "CG" && method != "LU" && method != "cho") {
					cout << "The following dense solvers are available: " << endl << endl;
					cout << " - Jacobi (type 'jac')" << endl;
					cout << " - Gauss-Seidel (type 'GS')" << endl;
					cout << " - Conjugate Gradient (type 'CG')" << endl;
					cout << " - LU decomposition (type 'LU')" << endl;
					cout << " - Cholesky decomposition (type 'cho')" << endl;
					cin >> method;
					if (method == "jac") {
						auto start = std::chrono::high_resolution_clock::now();
						A->jacobi(*A, *B, *X);
						cout << "Solution 'x': " << endl;
						X->printMatrix();
						auto finish = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed = finish - start;
						std::cout << "Time taken: " << elapsed.count() << " s\n" << endl << endl;
						delete A;
						delete X;
						delete B;
					}
					else if (method == "LU") {
						//A->LU_solve(*A, *B, *X);
						//X->printMatrix();
					}
					else if (method == "GS") {
						auto start = std::chrono::high_resolution_clock::now();
						A->gauss_seidel(*A, *B, *X);
						cout << "Solution 'x': " << endl;
						X->printMatrix();
						auto finish = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed = finish - start;
						std::cout << "Time taken: " << elapsed.count() << " s\n" << endl << endl;
						delete A;
						delete X;
						delete B;
					}
					else if (method == "CG") {

					}
					else if (method == "cho") {
						auto start = std::chrono::high_resolution_clock::now();
						A->cholesky_solve(*A, *B, *X);
						cout << "Solution 'x': " << endl;
						X->printMatrix();
						auto finish = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed = finish - start;
						std::cout << "Time taken: " << elapsed.count() << " s\n" << endl << endl;
						delete A;
						delete X;
						delete B;
					}
				}
			}
			else if (mat_type == "sparse") {
				while (method != "GS" && mat_type != "CG") {
					auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "random"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
					auto* A = new CSRMatrix<double>(*A_dense); // WE CONVERT OUR KNOWN DENSE A MATRIX TO SPARSE
					auto* X = new Matrix<double>(rows, 1, true);
					for (int j = 0; j < rows * rows; j++) { X->values[j] = 0; }
					auto* B = new Matrix<double>(rows, 1, true);
					for (int j = 0; j < rows * rows; j++) { B->values[j] = 1; }
					cout << "The system Ax = b will be solved, where: " << endl << endl;
					cout << "A (sparse): " << endl;
					A->printMatrix();
					cout << "b: " << endl;
					B->printMatrix();
					cout << "The following sparse solvers are available: " << endl;
					cout << " - Gauss-Seidel (type 'GS')" << endl;
					cout << " - Conjugate Gradient(type 'CG')" << endl;
					cin >> method;
					if (method == "GS") {
						auto start = std::chrono::high_resolution_clock::now();
						A->gauss_seidel(*A, *B, *X);
						cout << "Solution 'x': " << endl;
						X->printMatrix();
						auto finish = std::chrono::high_resolution_clock::now();
						std::chrono::duration<double> elapsed = finish - start;
						std::cout << "Time taken: " << elapsed.count() << " s\n" << endl << endl;
						delete A;
						delete X;
						delete B;
					}
					else if (method == "CG") {
						//conjugate gradient stuff
					}



				}

			}
			else {
				cout << "Please enter 'dense' or 'sparse' or press Ctrl + C to exit" << endl;
			}


		}

	cout << "Would you like to solver another system?(Y/N): ";
	cin >> choice;
	}
}


int main_old()
{
	//// GUI enabled by default: Allows user to generate a pseudo-random matrix A
	//// and solve the system Ax = b (where b is a vector of 1's)
	//// Can be used to test all solvers (dense and sparse)
	GUI();

	//// !! WATCH OUT !! If any of the code below is uncommented, please make sure to uncomment the delete statements at the bottom
	//// Instead of using the GUI, the sections of code below can be uncommented to edit the 
	//// linear system by hand. Note that most solvers only work for Symmetric Positiv Definite Matrices (SPDs)


	//// INITIALISATION - MAKE SURE TO UNCOMMENT SECTION BELOW TO RUN DENSE/SPARSE SOLVERS
	//// This section is used to initialise objects to be passed to the solvers later on.
	//// The constructor below takes the following parameters:
	//// (number of rows, maximum diagonal value, minimum diagonal value, minimum value off diagonal, keyword (random or else))
	//// The matrix A is of size NxN. N is defined as 'rows' below:
	// int rows = 5;
	//// Now, for example, let's generate a dense matrix by:
	// auto* A = new Matrix<double>(rows, 500, 490, 2, 1, "not random");
	//// Let's initialise x and b (and fill x with zeros, b with 1's)
	// auto* X = new Matrix<double>(rows, 1, true);
	// for(int j = 0; j < rows * rows; j++){X->values[j] = 0;}
	// auto* B = new Matrix<double>(rows, 1, true);
	// for(int j = 0; j < rows * rows; j++){B->values[j] = 1;}
	//// Now the solvers can be called


	//// 1. DENSE SOLVERS

	//// 1.1. DENSE JACOBI
	//// Takes A, x and b as pointers and solves for x (modified directly)
	// A->jacobi(*A, B, X);
	//// The result can be printed by calling:
	// X->printMatrix();

	//// 1.2. DENSE GAUSS-SEIDEL
	//// Takes A, x and b as pointers and solves for x (modified directly)
	// A->gauss_seidel(*A, *B, *X);
	//// The result can be printed by calling:
	// X->printMatrix();

	//// 1.3. DENSE LU

	//// 1.4. DENSE CONJUGATE GRADIENT

	//// 1.5. DENSE CHOLESKY DECOMPOSITION ('A' must be symmetric)
	//// Takes A, x and b as pointers and solves for x (modified directly)
	// A->cholesky_solve(*A, *B, *X);
	//// The result can be printed by calling:
	// X->printMatrix();


	//// 2. SPARSE SOLVERS
	//// The dense matrix initialised above is converted into a sparse matrix by:
	// auto* A = new CSRMatrix<double>(*A_dense); // WE CONVERT OUR KNOWN DENSE A MATRIX TO SPARSE

	//// 2.1. SPARSE GAUSS-SEIDEL
	//// Takes A, x and b as pointers and solves for x (modified directly)
	// A->gauss_seidel(*A, *B, *X);
	//// The result can be printed by calling:
	// X->printMatrix();

	//// 2.2. SPARSE CONJUGATE GRADIENT


	// delete A;
	// delete X;
	// delete B;

	return 1;
}
