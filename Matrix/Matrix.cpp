#include "Matrix.h"
#include <string>
#include <cstdlib>

using namespace std;

// Default constructor
template<class T>
Matrix<T>::Matrix() {}

template<class T>
// :: means what class it's coming from
// Initialization list does everything in the block {}, but can do more (like set constants)
Matrix<T>::Matrix(int rows_in, int cols_in, bool preallocate)
: rows(rows_in), cols(cols_in), size_of_values(rows_in * cols_in), preallocated(preallocate)

{
	this->rows = rows;		// like `self`: access own variables (`rows = rows` also works, but not as safe!)
	this->cols = cols;		// like `self`: access own variables

	if (this->preallocated)	// `if (preallocated) {}`
		this->values = new T[this->size_of_values];
}

// constructor here we have to manage our own data
template<class T>
Matrix<T>::Matrix(int rows_in, int cols_in, T* values_ptr)
: rows(rows_in), cols(cols_in), size_of_values(rows_in * cols_in), values(values_ptr)
{}

// COPY CONSTRUCTOR
template<class T>
Matrix<T>::Matrix(const Matrix& old_obj) {
	*this = old_obj;
	// This `new` here is now detaching our `.values` attribute from the `old_obj`'s
	// to a new slab of memory and then writing in our old_obj values to this into
	// this new slab of memory. So the `delete` called during the constructor deletes
	// this slab of memory but keeps the (dangling) pointer.
	this->values = new T[old_obj.size_of_values];
	for (int i = 0; i < old_obj.size_of_values; ++i) {
		this->values[i] = old_obj.values[i];
	}
}

// Identitity matrix constructor
template<class T>
Matrix<T>::Matrix(int len, std::string keyword) : rows(len), cols(len)
{
	if (keyword == "eye") {
		this->preallocated = true;
		this->values = new T[(rows) * (cols)];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i == j) {						// fill with 1 if on diagonal (i==j)
					this->values[i * len + j] = 1.;
				}
				else {
					this->values[i * len + j] = 0.;
				}
			}
		}
	}
	else if (keyword == "zeros") {
		this->preallocated = true;
		this->values = new T[(rows) * (cols)];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				this->values[i * len + j] = 0.;		// fill entire matrix with zeros
			}
		}
	}
	else {
		cout << "Fatal: Invalid keyword argument\n";
	}
}

// makes an SPD matrix
template <class T>
Matrix<T>::Matrix(int len, int diag_max, int diag_min, int offdiag_max, int offdiag_min, string keyword) : rows(len), cols(len), diag_max(diag_max), diag_min(diag_min), offdiag_max(offdiag_max), offdiag_min(offdiag_min), size_of_values(len*len)
{
	if (keyword == "random")
	{
		this->preallocated = true;
		this->values = new T[this->rows * this->cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i == j) {						// fill with value if on diagonal (i==j)
					this->values[i * this->cols + j] = rand() % diag_max + diag_min;
				}
				else {
					this->values[i * this->cols + j] = rand() % offdiag_max + offdiag_min;
				}
			}
		}
	}
	else {
		this->preallocated = true;
		this->values = new T[this->rows * this->cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i == j) {						// fill with value if on diagonal (i==j)
					this->values[i * this->cols + j] = 4;
				}
				else if (i == j + 1 || i == j - 1) {
					this->values[i * this->cols + j] = 1;
				}
				else if (i == j + 2 || i == j - 2) {
					this->values[i * this->cols + j] = 0.25;
				}
				else {

					this->values[i * this->cols + j] = 0;
				}
			}
		}


	}
}

template <class T>
Matrix<T>::Matrix(T* diag_ptrs, int len, std::string keyword) : rows(len), cols(len)
{
	if (keyword == "diag") {
		this->preallocated = true;
		this->values = new T[this->rows * this->cols];
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				if (i == j) {						// fill with value if on diagonal (i==j)
					this->values[i * this->cols + j] = diag_ptrs[i];
				}
				else {
					this->values[i * this->cols + j] = 0.;
				}
			}
		}
	}
	else {
		cout << "Fatal Error: Incorrect keyword argument. Did you mean 'diag'?" << endl;
	}
}


template<class T>
void Matrix<T>::printValues()
{
	for (int i = 0; i < this->size_of_values; ++i)
		std::cout << this->values[i] << " ";
	std::cout << std::endl;
}

template<class T>
void Matrix<T>::printMatrix()
{
	// Print explicitly assumed row-major order
	for (int i = 0; i < this->rows; ++i) {
		for (int j = 0; j < this->cols; ++j) {
			std::cout << this->values[j + i * this->cols] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// Implicitly assumes that mat_right and output exist
// A @ B = C  <=>  input @ mat_right = output
template<class T>
void Matrix<T>::matMatMult(const Matrix& mat_right, Matrix& output)	// Matrix class passed by reference
{
	// Make robust and check whether the data is input correctly using tests
	if (this->cols != mat_right.rows) // check dimensions on matrices allows for multiplication
	{
		std::cerr << "\nInput dimensions don't match for multiplication\n"; exit(1);
	}


	// Initiate results array with zeros
	for (int i = 0; i < this->rows * mat_right.cols; ++i)
		output.values[i] = 0;

	/* Note on matrix storage from https://sites.cs.ucsb.edu/~tyang/class/240a17/slides/Cache3.pdf
	A matrix is a 2D array of alements by memory addresses are "1D", where C++ is row major
		=> ` A(i,j) at A + i*n + j `
	Row major:
	[[0  1  2  3 ]
	 [4  5  6  7 ]
	 [8  9  10 11]
	 [12 13 14 15]]
	 Let's assume we have two levels in the hierarchy, fast (cache memory levels in cores) and slow (ram, which
	 is accessed separately), and all data is initially in slow memory, and stored in fast cached memory when
	 immediately after accessing, in case it's needed later.
	*/
	// During optimisation, we assume that fast (cached) memort is fast enough to hold three vectors
	// and that he cost of a fast memory access is zero.

	//Set values to zero beforehand
	for (int i = 0; i < this->rows * mat_right.cols; i++)
	{
		output.values[i] = 0.;
	}

	T* A = new T;
	for (int i = 0; i < this->rows; ++i) {
		for (int j = 0; j < mat_right.cols; ++j) {
			*A = output.values[i * output.cols + j];
			for (int k = 0; k < this->cols; ++k) {
				// the ordering affects the performance of the cache-perfomance - if we change the loop orderings
				// we might expect different performance
				// the += is problematic as it's instantiated to God-knows what, we must instantiate it first (done above).
				*A += this->values[i * this->cols + k] * mat_right.values[k * mat_right.cols + j];
			}
			output.values[i * output.cols + j] = *A;
			
		}
	}
	delete A;

	//T* acc00 = new T(0.);
	//T* acc01 = new T(0.);
	//T* acc10 = new T(0.);
	//T* acc11 = new T(0.);

	//int N = this->rows;
	//int ib = 24;
	//int kb = 64;
	//for (int ii = 0; ii < N; ii += ib)
	//{

	//		for (int j = 0; j < N; j += 2)
	//		{
	//			for (int i = ii; i < ii + ib; i += 2)
	//			{

	//					for (int k = 0; k < N; k++)
	//					{
	//						*acc00 += mat_right.values[(k * N) + j + 0] * this->values[(i + 0) * N + k];
	//						*acc01 += mat_right.values[(k * N) + j + 1] * this->values[(i + 0) * N + k];
	//						*acc10 += mat_right.values[(k * N) + j + 0] * this->values[(i + 1) * N + k];
	//						*acc11 += mat_right.values[(k * N) + j + 1] * this->values[(i + 1) * N + k];
	//					}
	//					output.values[(i + 0) * N + j + 0] = *acc00;
	//					output.values[(i + 0) * N + j + 1] = *acc01;
	//					output.values[(i + 1) * N + j + 0] = *acc10;
	//					output.values[(i + 1) * N + j + 1] = *acc11;
	//			}
	//		}
	//}

	//delete acc00;
	//delete acc01;
	//delete acc10;
	//delete acc11;
}

template <class T>
void Matrix<T>::matVecMult(T* input, T* output)
{
	// Set values to zero beforehand
	for (int i = 0; i < this->cols; i++)
	{
		output[i] = 0;
	}


	for (int i = 0; i < this->rows; i++) // loop over matrix rows
	{
		for (int k = 0; k < this->cols; k++) // loop over matrix cols
		{
			// calculate value of each element in out vector array from matrix-vector mult.
			output[i] += this->values[i * this->cols + k] * input[k];


		}
	}

}


template<class T>
void Matrix<T>::matMatAdd(const Matrix& mat_right, Matrix& output)	// Matrix class passed by reference
{
	// Make robust and check whether the data is input correctly using tests
	if (this->cols != mat_right.cols && this->rows != mat_right.rows) // check dimensions on matrices allows for multiplication
	{
		cout << "Matrix dimensions do not match for addition" << endl;
		exit(1);
	}
	else { // if rows and columns match, do matrix addition for all elements
		// Initiate results array with zeros
		for (int i = 0; i < output.size_of_values; ++i)
			output.values[i] = 0;

		// For all rows and columns, access every index over each array and add them together to make the output array
		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->cols; ++j) {
				output.values[i * this->cols + j] += this->values[i * this->cols + j] + mat_right.values[i * this->cols + j];
			}
		}
	}
}

template<class T>
void Matrix<T>::matMatSub(const Matrix& mat_right, Matrix& output)	// Matrix class passed by reference
{
	// Make robust and check whether the data is input correctly using tests
	if (this->cols != mat_right.cols && this->rows != mat_right.rows) // check dimensions on matrices allows for multiplication
	{
		cout << "Matrix dimensions do not match for addition" << endl;
		exit(1);

	}
	else { // if rows and columns match, do matrix addition for all elements
		// Initiate results array with zeros
		for (int i = 0; i < output.size_of_values; ++i)
			output.values[i] = 0;

		// For all rows and columns, access every index over each array and add them together to make the output array
		for (int i = 0; i < this->rows; ++i) {
			for (int j = 0; j < this->cols; ++j) {
				output.values[i * this->cols + j] -= this->values[i * this->cols + j] + mat_right.values[i * this->cols + j];
			}
		}
	}
}

template<class T>
Matrix<T>* Matrix<T>::operator* (const Matrix<T>& obj) // matrix pointer passed by reference
// being passed by reference
{
	// Declare a new_entry object and assign the sum of the
	// member variables of the current object (using this) and the parameter object (obj)
	// to the new_entry object's regVariable member and return the object as a result.
	// Instantiate another object

	auto* output_prealloc = new Matrix<T>(this->rows, obj.cols, true);
	this->matMatMult(obj, *output_prealloc);

	return output_prealloc;
}

template<class T>
Matrix<T>* Matrix<T>::operator+ (const Matrix<T>& obj) // matrix pointer passed by reference
// being passed by reference
{
	// Declare a new_entry object and assign the sum of the
	// member variables of the current object (using this) and the parameter object (obj)
	// to the new_entry object's regVariable member and return the object as a result.
	// Instantiate another object

	auto* output_prealloc = new Matrix<T>(this->rows, obj.cols, true);

	this->matMatAdd(obj, *output_prealloc);
	return output_prealloc;
}

template<class T>
Matrix<T>* Matrix<T>::operator- (const Matrix<T>& obj) // matrix pointer passed by reference
// being passed by reference
{
	// Declare a new_entry object and assign the sum of the
	// member variables of the current object (using this) and the parameter object (obj)
	// to the new_entry object's regVariable member and return the object as a result.
	// Instantiate another object

	auto* output_prealloc = new Matrix<T>(this->rows, obj.cols, true);

	this->matMatSub(obj, *output_prealloc);
	return output_prealloc;
}


template<class T>
Matrix<T>* Matrix<T>::transpose()
{
	auto* output = new Matrix<T>(this->rows, this->cols, true);
	for (int i = 0; i < this->rows; ++i) {
		for (int j = 0; j < this->cols; ++j) {
			output->values[i * this->cols + j] = this->values[j * this->cols + i];
		}
	}

	return output;
}

template<class T>
Matrix<T>::~Matrix()		// Cannot pass parameters to constructors (already knows everything it needs to delete itself)
{
	if (this->preallocated) {
		delete[] this->values;
	}
}

template<class T>
void Matrix<T>::swap_row(Matrix<T>& A, int i, int j) {
	/* Swap rows i and j of the matrix and vector b.*/
	if (i == j) {
		return;
	}

	//cout << "swapping rows " <<  i << " and " << j << endl;
	// Declare a new set of matrix rows through a hard copy
	// to prevent us from overwriting our current data.
	// These new pointer arrays store values that are to be swapped.
	auto* iA = new T[A.cols];

	// Fill the new (replacement) iA and ib matrices with their 
	// corresponding rows from A and b
	for (int k = 0; k < A.cols; ++k) {
		iA[k] = A.values[i * A.cols + k];
	}

	// Overwrite the rows
	for (int k = 0; k < A.cols; ++k) {
		A.values[i * A.cols + k] = A.values[j * A.cols + k];
	}

	// Fill in the missing row
	for (int k = 0; k < A.cols; ++k) {
		A.values[j * A.cols + k] = iA[k];
	}

	/// Get rid of iA as it is no longer needed.
	delete[] iA;
}


template<class T>
void Matrix<T>::forward_substitution(Matrix<T>& A, Matrix<T>& b, Matrix<T>& x) {
	/* Function to perform forwards substitution as second major part of solving
	Ax = b linear system.
	Inputs:
		- A in lower triangular form, and;
		- b that results from upper-triangle operations
	Outputs:
		- new x solution to linear system.*/


	int n = b.size_of_values;
	// Initialize x with zeros
	for (int i = 0; i < x.size_of_values; ++i) {
		x.values[i] = 0;
	}


	T s;
	// Unlike Drake, let's start at top and work our way down
	for (int k = 0; k < n; ++k) {
		// can do all of this with vector - vector multiplication
		// Declare double that will act as our fraction
		s = 0;
		for (int j = 0; j < n; ++j) {
			// Create fraction to add to x solution at every pivot point
			s += A.values[k * A.cols + j] * x.values[j];
		}
		x.values[k] = (b.values[k] - s) / A.values[k * A.cols + k];
	}
}


template<class T>
void Matrix<T>::back_substitution(Matrix<T>& A, Matrix<T>& b, Matrix<T>& x) {
	/* Function to perform back substitution as second major part of solving
	Ax = b linear system.
	Inputs:
		- A in upper triangular form, and;
		- b that results from upper-triangle operations
	Outputs:
		- new x solution to linear system.*/

	int n = b.size_of_values;

	// Initialize x with zeros
	for (int i = 0; i < x.size_of_values; ++i) {
		x.values[i] = 0;
	}

	// Like Drake, let's start at bottom and work our way up
	for (int k = n - 1; k >= 0; --k) {
		// can do all of this with vector - vector multiplication
		// Declare double that will act as our fraction
		T s = 0.;
		for (int j = k + 1; j < n; ++j) {
			// Create fraction to add to x solution at every pivot point
			s += A.values[k * A.cols + j] * x.values[j];
		}
		x.values[k] = (b.values[k] - s) / A.values[k * A.cols + k];
	}
};

template<class T>
void Matrix<T>::LU_Decomposition(Matrix<T>& iA, Matrix<T>& P_, Matrix<T>& L, Matrix<T>& U)
{
	/*		Deconstructs A (design matrix) into P, L and U within LU solver.
	Inputs:
		- iA - matrix to deconstruct (does not get overwritten) (N x N)
		- P_ (reference) - ptr to ones matrix of size (N x N) input to be used as permutation matrix
		- L (reference) - ptr to zeros matrix of size (N x N) input to be used as lower triangular matrix
		- U (reference) - ptr to matrix that's a copy of iA used as upper triangular matrix
	Output:
		- void
	Result:
		- P_ (altered)
		- L (altered)
		- U (altered)

	Benchmarked with previous code from ACSE-3.3 Linear Solvers, Matthew Piggot. */

	// TRY UNIQUE POINTERS HERE

	const int m = iA.rows;
	const int n = iA.cols;

	// CREATE LOCAL COPY
	// We don't want to change A matrix, so let's do a deep copy
	// on the heap and delete the UPPER TRIANGLE after the operations have been done.
	// Initialise L to be a zeros array
	// L = Matrix<T>(m, "zeros");	// Make [L]ower triangle of A initially a matrix of zeros
	// Initialise P_ to be the permutation matrix (currently `I`)
	// P_ = Matrix<T>(m, "eye");

	int kmax = 0;	// Declare a max row magnitude bookmark
	// Loop over each pivot row, except the last
	for (int k = 0; k < m - 1; ++k) {
		// Implement partial pivoting
		// Order rows from highest to lowest magnitude

		// Loop over all entries below the pivot and select the k with largest abs value
		for (int i = k; i < n; ++i) {
			// Measure the absolute value of the array
			if (abs(iA.values[kmax * iA.cols + k]) < abs(iA.values[i * iA.cols + k])) {
				kmax = i;
			}
		}
		// j here will be the entry counting down the matrix from out current partial pivot element
		// so add j onto k to get the row index

		// Swap the current pivot row (k) with the row that contains the largest
		// magnitude row below that exact pivot element
		// kmax is the  b i g   b o i  element and k is the pivot element
		this->swap_row(iA, kmax, k);
		// record this swap in the permutation matrix
		this->swap_row(P_, kmax, k);
		// record this swap in lower triangle matrix
		this->swap_row(L, kmax, k);


		for (int i = k + 1; i < m; ++i) {
			// only one double value, so fine to store on stack for time being
			T s = (iA.values[i * iA.cols + k]) / (iA.values[k * iA.cols + k]);
			for (int j = k; j < n; ++j) {
				iA.values[i * iA.cols + j] -= s * iA.values[k * iA.cols + j];
				L.values[i * L.cols + k] = s;
			}
		}
		// Copy iA values into U
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				U.values[i * m + j] = iA.values[i * m + j];
			}
		}
	}

	// For some reason this is overwriting it here, but when we return it to the calling function, the eye matrix
	// isn't there anymore
	//L = *L + Matrix<T>(m, "eye");		// Add the identity matrix... haven't yet overwritten the += operator
	// ...must be a scoping thing...

	// Manually fill in identity matrix
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			if (i == j) {
				L.values[i * n + j] += 1.;
			}
		}
	}
};

template<class T>
void Matrix<T>::LU_solve(Matrix<T>& A, Matrix<T>& b, Matrix<T>& output)
{
	/* Wrapper for Linear Solver for Dense Matrices `Ax = b`, where we wish to solve for x*/
	// Before you build a shed, you have to check your toolbox
	//		- Mr. Ore's A-Level Mathematics class explaining which SUVAT equation to choose
	if (A.rows != A.cols) { cout << "\nFatal Error: Non-square matrix"; exit(-1); }
	if (A.cols != output.rows) {
		cout << "\nFatal Error: Matrix Dimensions do not match";
		cout << "Input matrix A of size (" << A.rows << ", " << A.cols;
		cout << ") cannot be solved for a vector of size (" << output.rows << ", " << output.cols << ").\n";
		exit(-1);
	}

	Matrix<T>* P = new Matrix<T>();		// matrix multiplication doesn't yet return unique pointers
	unique_ptr<Matrix<T>> P_(new Matrix<T>(A.rows, "eye"));
	unique_ptr<Matrix<T>> iA(new Matrix<T>(A));
	unique_ptr<Matrix<T>> ib(new Matrix<T>(b));
	unique_ptr<Matrix<T>> L(new Matrix<T>(A.rows, "zeros"));
	unique_ptr<Matrix<T>> U(new Matrix<T>(A));

	// Deconstruct the A matrix into upper and lower components
	this->LU_Decomposition(*iA, *P_, *L, *U);		// P_, L, U have been written into now
													// (A is unchanged)

	P = *P_ * *ib;								// Use overwritten * operator to do matVecMult;
	
	// Create new matrix to store answer from forward substitution
	unique_ptr<Matrix<T>> y(new Matrix<T>(b.rows, b.cols, true));

	// Substitution operators
	// First on a depository array which maps the changes from L to U triangle matrices
	// y vector initialized to zero within function
	this->forward_substitution(*L, *P, *y);
	// Secondly on the eventual output
	// output vector initialized to zero within function
	this->back_substitution(*U, *y, output);

	// Delete all scoped variables from the heap

	delete P;

}

template<class T>
void Matrix<T>::upper_triangle_pp(Matrix<T>& A, Matrix<T>& b)
{
	/* Function that converts A into upper triangular form using row operations,
	performing the same operations on vector b.

	This function uses partial pivoting to overcome the problems induced by zeros
	(especially in sparse matrices). This is when we swap the current pivot row with
	the row below that has the largest magnitude entry in the pivot column.*/

	int n = b.size_of_values;
	// Check that A is square
	if (A.rows != A.cols) {
		cout << "\nFatal Error: Matrix not square.\n";
		exit(-1);
	}
	// Check that A and b dimensions are compatible
	if (A.cols != b.rows) {
		cout << "\nFatal Error: Matrix and vector dimensions are not compatible.\n";
		exit(-1);
	}

	int kmax;	// Declare a max row magnitude bookmark
	// Loop over each pivot row, except the last
	for (int k = 0; k < n - 1; ++k) {
		// Implement partial pivoting
		// Order rows from highest to lowest magnitude
		kmax = k;
		// Loop over all entries below the pivot and select the k with largest abs value
		for (int i = k + 1; i < n; ++i) {
			// Measure the absolute value of the array
			if (abs(A.values[kmax * A.cols + k]) < abs(A.values[i * A.cols + k])) {
				kmax = i;
			}
		}

		// cout << "k: " << k << ", kmax: " << kmax << endl;
		// Swap the current pivot row (k) with the row that contains the largest
		// magnitude row below that exact pivot element

		this->swap_row(A, kmax, k);
		this->swap_row(b, kmax, k);

		// cout << "Completed swapping for pivot number " << k << endl;
		for (int i = k + 1; i < n; ++i) {
			T s = A.values[i * A.cols + k] / A.values[k * A.cols + k];
			for (int j = k; j < n; ++j) {
				// cout << "i: " << i << ", j: " << j << endl;
				A.values[i * A.cols + j] -= s * A.values[k * A.cols + j];
			}
			b.values[i] -= s * b.values[k];
		}
	}
}

template<class T>
void Matrix<T>::gaussian_elimination_solve(Matrix<T>& A, Matrix<T>& b, Matrix<T>& x)
{	/*		SOLVES `Ax = b` LINEAR SYSTEM FOR `x`.
	Wrapper for partial pivoting, upper triangulation and backward substitution functions.
	Requires a declared x as `new` class pointer on heap, user is responsible from
	hereon to delete this new pointer from their heap memory space.
	Inputs:
		- A: Square matrix of size (m x m)
		- b: Target vector of size (m x 1)
	Outputs:
		- x: Solution vector of size (m x 1)	*/
	unique_ptr<Matrix<T>> iA(new Matrix<T>(A));	// make copy of A so that we don't overwrite it
	unique_ptr<Matrix<T>> ib(new Matrix<T>(b));	// make copy of b so that we don't overwrite it

	this->upper_triangle_pp(*iA, *ib);
	this->back_substitution(*iA, *ib, x);

}

template<class T>
bool Matrix<T>::isSymmetric()
{
	//returns true if matrix is symmetric, false otherwise
	Matrix<T>* thisTranspose = this->transpose();
	for (int i = 0; i < this->rows; i++){ 
		for (int j = 0; j < this->cols; j++){
			if (this->values[i*this->cols + j] != thisTranspose->values[i*thisTranspose->cols + j]){
				return false; 
			}
		}
	}
	return true;
}

template<class T>
Matrix<double>* Matrix<T>::choleskyDecomposition() { 
    // Perform Cholesky Decomposition on matrix 
    // Return Lower triangular matrix
    // The Upper Triangular Matrix can be obtained by calling transpose()
    // Usage example:
    // int rows = 5, cols = 5;
    // auto* A = new Matrix<int>(rows, cols, true);                 // A
    // auto* result = A->choleskyDecomposition();                   // L 
    // Matrix<long double>* resultTranspose = result->transpose();  // L^t

    // Algorithm taken from: https://algowiki-project.org/en/Cholesky_decomposition

    //Initialise empty matrix => Lower triangular will be stored here
	Matrix<T>* lower = new Matrix<T>(this->rows, this->cols,  true); // change back to long double

    // Decomposing matrix into Lower Triangular 
    for (int i = 0; i < this->rows; i++) { 
        for (int j = 0; j <= i; j++) { 
            T sum = 0; 
            if (j == i){ // diagonals  L(j, j)
                for (int k = 0; k < j; k++) {
                    sum += pow(lower->values[j*lower->cols + k], 2); 
                }
                lower->values[j*lower->cols + j] = sqrt(this->values[j*this->cols + j] - sum); 
            } else { 
                // non-diagonals L(i, j) 
                for (int k = 0; k < j; k++){
                    sum += lower->values[i*lower->cols + k] * lower->values[j*lower->cols + k]; 
                }
                lower->values[i*lower->cols + j] = T((this->values[i*this->cols + j] - sum) / lower->values[j*this->cols + j]); 
            } 
        } 
    }
    return lower;
} 


template<class T> 
void Matrix<T>::cholesky_solve(Matrix<double>& a, Matrix<double>& b, Matrix<double>& x) { 
    // Check input dimensions are OK.
    if (a.rows != a.cols) { cout << "\nFatal Error: Non-square matrix"; exit(-1); }
	if (a.cols != x.rows) {
		cout << "\nFatal Error: Matrix Dimensions do not match";
		cout << "Input matrix A of size (" << a.rows << ", " << a.cols;
		cout << ") cannot be solved for a vector of size (" << b.rows << ", " << b.cols << ").\n";
		exit(-1);
	}

    // Cholesky decomposition turns A*x = b into L*L^t*x = b
    // Let L^t = y => L*y = b; AND L^t*x = y; Definition of 'y' is below:
	Matrix<T>* y = new Matrix<T>(a.rows, b.cols, true);

    for (int j = 0; j < (a.rows * b.cols) ; j++) {y->values[j] = 0;}

    // Find L and L^t
    Matrix<T>* L = a.choleskyDecomposition();
    Matrix<T>* Lt = L->transpose();

    //Solve L*y = b
    y->forward_substitution(*L, b, *y);
    //Solve L^t*x = y; vector 'b' is reused to store the solution
    y->back_substitution(*Lt, *y, x);

	delete L;
	delete Lt;
	delete y;
}

template <class T>
void Matrix<T>::gauss_seidel(Matrix<T>& a, Matrix<T>& b, Matrix<T>& x) {
//  Calculates GS iterations
//  Parameters:
//    Input:
//      <int> iterations: Number of GS iterations to compute
//      <int> rows: Number of rows of A
//      <int> cols: Number of cols of A 
//      For the system Ax = b:
//      Matrix<double>* a, x_0, b 
//    Output:
//      Matrix<double> x: Update x vector after GS iteration
//  Assumptions:
//    Ax = b has a unique solution
//    A has no zeros on its main diagonal
//    A is diagonally dominant
    double tolerance = 0.00001;
    T norm = 1.0; //initialise to value larger than tolerance
    int rows = a.rows; // could get rid of these variables
    int cols = a.cols; 
    auto n = rows;
    auto* y = new Matrix<T>(rows, cols, true);
    auto* x_prev = new Matrix<T>(rows, 1, true); // used to compute tolerance
    for (int j = 0; j < rows * cols; j++){y->values[j] = 0;}
    auto omega = 1; //change to implement SOR 
    //perform Gauss_Seidel iteration until tolerance is reached
    while (norm > tolerance){
        //design choice: copying array to compute tolerance -> overhead
        //possible better alternative: fixed number of iterations? OK for big matrices
        //not so good for small matrices
        //copy constructor
	    for (int i = 0; i < x.size_of_values; ++i) {
		    x_prev->values[i] = x.values[i];
	    }
        norm = 0;
        // Gauss-Seidel algorithm taken from: https://www.geeksforgeeks.org/gauss-seidel-method/
        for (int i = 0; i < n; i++){
            y->values[i] = omega*(b.values[i] / a.values[i+ i*rows]);
            for (int j = 0; j < n; j++){
                if (j == i) continue;
                y->values[i] = omega*(y->values[i] - ((a.values[i*rows + j] / a.values[i+i*rows]) * x.values[j]));                
                x.values[i] = y->values[i] + (1-omega)*x.values[i];                
            }
        }
        // Calculate 2-norm between current and previous iteration
        for (int i = 0; i < x.size_of_values; ++i) {
		    x_prev->values[i] = x.values[i] - x_prev->values[i]; //might replace this by overloaded operator -
            norm += pow(x_prev->values[i], 2)/x_prev->rows;
	    } 
		norm = sqrt(norm);
    }
}

template <class T>
void Matrix<T>::jacobi(Matrix<double>& a, Matrix<double>& b, Matrix<double>& x) {
//  Calculates one Jacobi iteration
//  Parameters:
//    Input:
//      <int> rows: Number of rows of A
//      <int> cols: Number of cols of A 
//      For the system Ax = b:
//      Matrix<double>* a, x_0, b 
//    Output:
//      Matrix<double> x_new: Update x vector after Jacobi iteration
//  Assumptions:
//    Ax = b has a unique solution
//    A has no zeros on its main diagonal
//    A is diagonally dominant 

    // norm and tolerance used to determine number of iterations
    double norm = 1.0;
    double tolerance = 10e-10;

    int rows = a.rows;

    //updated values after Jacobi iteration
    auto* x_new = new Matrix<double>(rows, 1, true);
    //difference between previous and current values
    auto* x_diff = new Matrix<double>(rows, 1, true);

    //loop until hitting tolerance
    while (norm>tolerance){
        norm = 0;
        //iterate through rows - Jacobi algorithm
        for (int i = 0; i < rows; i++ ){
            x_new->values[i] = b.values[i];
            for (int j = 0; j < rows; j++ ){
                if (j != i ){
                    x_new->values[i] = x_new->values[i] - a.values[i+j*rows] * x.values[j];
                }
            }
            x_new->values[i] = x_new->values[i] / a.values[i+i*rows];
        }

        // calculate difference between current and previous iterations 
        // update norm value
        for (int i = 0; i < rows; i++ ){
            x_diff->values[i] = x_new->values[i] - x.values[i];
            x.values[i] = x_new->values[i];
            norm += pow(x_diff->values[i], 2)/x_diff->rows;
        }
		norm = sqrt(norm);
    }
};

template <class T>
void Matrix<T>::conjgrad(T* x_guess, T* b, double tol, int maxiter, std::string keyword)
{

	/*
	 CONJUGATE GRADIENT SOLVER FOR A LINEAR SYSTEM Ax = b, where A is a vector, x and b are vectors.

	- INPUTS: (1) pointer: x initial guess vector; (2) pointer: known b vector, (3) double: solver tolerance, (4) maximum no.  of iterations
	- NOTE: the input matrix MUST be symmetric positive definite
	- OUTPUS: None (but passing in pointers so alters input x_guess outside the scope of the function)
	*/
	

	//// 1. INITIALISE SOME VARIABLES ////
	// (i) some arrays
	auto* r = new T[this->cols]; // residual array
	auto* r_old = new T[this->cols]; // array for storing previous iterations residuals
	auto* r_part = new T[this->cols]; // array for storing part of the residual calculation: A * x_guess
	auto* p = new T[this->cols]; // array for storing CG search direction 
	auto* store = new T[this->cols]; // array for storing A*p in alpha calculation 

	// (ii) some scalars
	double alpha, beta; // alpha and beta some scalars we will calculate
	int count = 0; // loop counter to give number of iterations taken 
	double beta_top = 0; // for beta dot products (numerator of beta)
	double beta_bot = 0; // for p * A * p (denominator of beta)
	double r_dot = 0; // for  r*r dot product, numerator for alpha coefficient 
	double bottom = 0; // denominator for alpha coefficient 
	double err_2norm = 0;  // calculating L2-norm to compare with user-defined tolerance


	//// 2. SOME INITIAL CALCULATIONS (PRIOR TO ITERATING) //// 
	this->matVecMult(x_guess, r_part); // A * x_guess = b_current, store in r_part array

	for (int i = 0; i < this->cols; i++)
	{
		r[i] = b[i] - r_part[i]; // residuals for first iteration, r = b - A*x_guess
		p[i] = r[i]; // search direction p is equal to the residual array r for the first iteration only 
	}

	this->matVecMult(p, store); // matrix-vector multiplicatation: A*p, store result in "store" array 

	for (int i = 0; i < this->cols; i++)
	{
		r_dot += r[i] * r[i]; // calculating vector-vector product r*r scalar
		bottom += p[i] * store[i]; //we need p* A* p = p[i] * store[i]
		err_2norm += r[i] * r[i]; // square of each residual 
	}

	alpha = r_dot / bottom; // alpha scalar coefficient 
	err_2norm = sqrt(err_2norm / this->cols); // RMS error - if less than the tolerance then we will not iterate


	//// 3. NOW ITERATE & SOLVE FOR X //// 
	while (count < maxiter && tol < err_2norm) // stops when either (a) we reach the max. number of iterations or (b) we reach our specified tolerance (RMS error)
	{
		// (i) calculate a new guess for our x array, x_guess
		for (int i = 0; i < this->cols; i++) // as cols of matrix A must be equal to rows of vector x (for dimensions to match)
		{
			x_guess[i] = x_guess[i] + alpha * p[i]; // we update our x_guess. x_guess = x_guess(old) + alpha * p
		}


		// (ii) copy old residual & calculate A*p mat-vec product
		for (int i = 0; i < this->cols; i++) // loop over each element of the array
		{
			r_old[i] = r[i]; // store the old residuals array, as we need it for another calculation later
		}
		this->matVecMult(p, store); // matrix-vector multiplicatation: A*p, store result in "store" array 


		// (iii) calculate our new residual array, r
		for (int i = 0; i < this->cols; i++)
		{
			r[i] = r[i] - alpha * store[i]; // now we update our residual vector
		}

		 // (iv) calculate our beta scalar coefficient
		beta_top = 0;
		beta_bot = 0;
		for (int i = 0; i < this->cols; i++)
		{
			beta_top += r[i] * r[i]; // dot product of r * r
			beta_bot += r_old[i] * r_old[i]; // dot product of r_old * r_old (residuals from previous iteration)
		}
		beta = beta_top / beta_bot; // calculating beta coefficient


		// (v) calculate our new search direction array, p
		for (int i = 0; i < this->cols; i++)
		{
			p[i] = r[i] + beta * p[i]; // p = r + beta *p(old)

		}

		// (vi) calc new alpha scalar coefficient 
		this->matVecMult(p, store); // matrix-vector multiplicatation: A*p, store result in "store" array 
		bottom = 0; // set sums of numerator and denominator of alpha coefficient to zero
		r_dot = 0; // 
		err_2norm = 0; // must set all back to zero for each loop as summing up
		for (int i = 0; i < this->cols; i++)
		{
			//r_dot += r[i] * r[i]; // calculating vector-vector product r*r scalar BOTH SAME AS BETA TOP = beta_top
			bottom += p[i] * store[i]; //we need p* A* p = p[i] * store[i]
			//err_2norm += r[i] * r[i]; // square of each residual  = beta_top
		}
		alpha = beta_top / bottom; // alpha scalar coefficient 

	   // (vii) calculate new RMS error of residual vector
		//err_2norm = sqrt(err_2norm / this->cols); // RMS error 
		err_2norm = sqrt(beta_top / this->cols); // RMS error 


		// (viii) Add one to iteration count & if maximum number of iterations reached then give error and exit
		count++;
		/// 
		if (count == maxiter)
		{
			cerr << "Convergence not reached after " << maxiter << " iterations. Exiting. " << endl;
			exit(1);
		}


	}


	//// 6. DELETE OUR POINTER ARRAYS FROM HEAP //// 
	delete[] r;
	delete[] r_old;
	delete[] r_part;
	delete[] p;
	delete[] store;

}
