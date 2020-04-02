#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <string>
#include <memory>
#include <tuple>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <ctime>

using namespace std;
// We aim to template what stores the values of interest
template<class T>
class Matrix
{
public:
	// Public methods
	// Public variables
	///// Size of our matrix
	T* values = nullptr;	// explicitly using the c++11 nullptr convention
	int rows = -1;
	int cols = -1;
	int size_of_values = -1;

	// Create a default constructor such that you can create static instances of a class
	Matrix();
	// MUST USE ROW-MAJOR ORDER!!
	// OR at least include functionality to change it from col-major to row-major
	// Constructor - pass in variables
	// where we want to preallocate memory - own our own memory
	Matrix(int rows, int cols, bool preallocate);
	// where we have already preallocated memory outside
	Matrix(int rows, int cols, T* values_ptr);	// we know that there's no point in passing in the
												// preallocated boolean if there's already a pointer
	// Copy constructor
	Matrix(const Matrix& old_obj);

	// Identity and zeros matrix constructor
	Matrix(int len, std::string keyword);

	// Construct diagonal matrix with given diagonal values
	Matrix(T* diag_ptrs, int len, std::string keyword);

	// Construct SPD matrix
	Matrix(int len, int diag_max, int diag_min, int offdiag_max, int off_diag_min, string keyword);

	// Documentation
	virtual void printValues();
	// Documentation
	virtual void printMatrix();

	// Code to do matrix-matrix multiplication using cache optimisation
	virtual void matMatMult(const Matrix& mat_right, Matrix& output);		// pass by reference to do a *shallow copy*

	// Code to do vector-matrix multiplication
	virtual void matVecMult(T* input, T* output);

	// Code to do matrix-matrix addition
	virtual void matMatAdd(const Matrix& mat_right, Matrix& output);
	// Code to do matrix-matrix subtraction
	virtual void matMatSub(const Matrix& mat_right, Matrix& output);
	
	// Define a multiplication operator that overwrites the * operator and makes it do matrix multiplication
	// This creates a raw pointer on the heap which the user is now responsible for.
	virtual Matrix* operator* (const Matrix& obj);
	
	// Define a multiplication operator that overwrites the * operator and makes it do matrix addition
	// This creates a raw pointer on the heap which the user is now responsible for.
	virtual Matrix* operator+ (const Matrix& obj);

	// Define a subtraction operator that overwrites the - operator and makes it do matrix subtraction
	// This creates a raw pointer on the heap which the user is now responsible for.
	virtual Matrix* operator- (const Matrix& obj);

	// Function that transposes the matrix and returns a new matrix
	virtual Matrix* transpose();

	// Row swapping function for partial pivoting algorithm
	virtual void swap_row(Matrix<T>& A, int i, int j);

	// Upper triangulation methon to convert matrix and corresponding target vector into upper-triangle matrix
	virtual void upper_triangle_pp(Matrix<T>& A, Matrix<T>& b);

	// Deconstructs A into P, L, U matrices
	virtual void LU_Decomposition(Matrix<T>& iA, Matrix<T>& P_, Matrix<T>& L, Matrix<T>& U);

	// Forward substitution algorithm for LU decomposition
	virtual void forward_substitution(Matrix<T>& A, Matrix<T>& b, Matrix<T>& x);

	// Backward substitution algorithm for LU decomposition
	virtual void back_substitution(Matrix<T>& A, Matrix<T>& b, Matrix<T>& x);

	// Gaussian elimination algorithm (basic version of LU decomposition for less than one system to solve)
	virtual void gaussian_elimination_solve(Matrix<T>& A, Matrix<T>& b, Matrix<T>& output);
	
	// Gauss Seidel solver
	virtual void gauss_seidel(Matrix<T>& a, Matrix<T>& b, Matrix<T>& x);
	
	// Find if matrix is symmetric
	bool isSymmetric();

	// Wrapper to perform Cholesky decomposition
	Matrix<double>* choleskyDecomposition();

	// Solve system using Cholesky decomposition
	void cholesky_solve(Matrix<double>& a, Matrix<double>& b, Matrix<double>& x);
	
	// Jacobi solver
	void jacobi(Matrix<double>& a, Matrix<double>& b, Matrix<double>& x);

	// Wrapper to solve linear system using LU decomposition
	virtual void LU_solve(Matrix<T>& A, Matrix<T>& b, Matrix<T>& output);

	// Conjugate gradient descent algorithm
	void conjgrad(T* x_guess, T* b, double tol, int maxiter);

	// Destructor
	virtual ~Matrix();	// virtual means that an inherited class can go on and write its own destructor method
	// pure virtual function: virtual ~Matrix() = 0 // => inherited class HAS to write its own constructor.


protected:
	bool preallocated = false;	// Defaulted to False
	int offdiag_max = 1; // initialise values for SPD matrix constructor
   	int offdiag_min = 1;
   	int diag_min = 1;
   	int diag_max = 1;

};
#endif // MATRIX_H
