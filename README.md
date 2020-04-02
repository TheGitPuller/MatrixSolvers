# acse-5-assignment-git_pullers 

# Implementing Linear Solvers

This advanced programming project aims to implement, test and optimize linear solvers for dense and sparse matrices with a view towards high performance computing and effective library design. This program is written in C++ and contains two libraries that pertain to dense and sparse matrices (each containing their respective `.cpp` and `.h` files), with the program running on the `main.cpp` which seeks to highlight this program's functionality.

## Getting Started

This program assumes a working understanding of matrices, matrix types and linear solvers for Symmetric Positive Definite (SPD) systems.

## Prerequisites

### For Linux/Mac users:
* Your compiler of C++ choice (it is recommended to use g++ (*GNU Project*) for its optimization capabilities).

### For Windors users:
* Microsoft Visual Studio Community IDE, or;
* Any choice of text Editor and Compiler (g++ recommended)

## Installing

* Clone the GitHub repository to your local machine.

## Running the program

### Running `main.cpp`:

* For Windows users:
	* For Microsoft Visual Studio Community IDE users:
		-- Open the Visual Studio in the cloned repository in the IDE;
		-- Run using the 'Release' version.

	* Run on Windows Subsystem for Linux (WSL) following Linux instructions below.

* For Linux/Mac Users
	* Enter the `Matrix/` directory in the cloned repository and compile the code using your compiler of choice in BASH terminal. We recommend using g++ with optimization flags turned on (namely `-O3` to optimize for execution time):

	```
	$> g++ -O3 main.cpp -o [my_file_name]
	```

	And run using:

	```
	$> ./[my_file_name]
	```

This runs the GUI to allow the user to easily run various linear solvers for the matrix system and configuration of choice, and to demonstrate the other functionality of the Matrix classes.
Running `main.cpp` also gives the user the opportunity to see everything else that the linear solvers library can do, namely:
* Time tests;
* Accuracy tests;
* Exhibitionary functions that demonstrate the functionality of the software and how to run each method.

Alongside the header files for Matrix and CSRMatrix classes, the latter function acts as an excellent examples list for how to run the solver suite!

## Matrix Types:

### Declaring Matrices:

#### Dense Matrices:
These matrices have every value at every index stored, regardless of whether the value is zero.

Example declaration:
```
int rows = 4;
int cols = 4;

double A_data[] = { 1., 0., 3., 7., 2., 1., 0., 4., 5., 4., 1., -2., 4., 1., 6., 2. };

double* A_vals = A_data;

auto* dense_mat = new Matrix<double>(rows, rows, A_vals);

[...]

delete dense_mat;
```

This declares a dense matrix, constructing its matrix from the values pointer (`A_vals`), which have been entered in **row-major order**.

The preallocation constructor may also be used to reserve memory space of a dense matrix in memory space to be filled in later by the user or a function when passed in as an argument.

For example:
```
auto* output_matrix = new Matrix<double>(rows, cols, true);

[...]

delete output_matrix;
```

#### Sparse Matrices:

This class inherits from dense class, such that its inherited classes have been polymorphed to suit CSR Matrices.

Sparse matrices are ones in which the majority of elements are zero. To save memory, we allocate its values using *Compressed Sparse Row* format.

Assuming knowledge of CSR Formatting, an example implementation is:

```
int rows = 4;		// number of rows
int cols = 4;		// number of columns
int nnzs = 6;		// number of non-zeros
double val_data[] = { 1., 2., 5., 3., 1., 3. };
double* val_input = val_data;			// Pointer to input values data
int col_data[] = { 2, 3, 1, 2, 0, 1 };
int* col_input = col_data;				// Pointer to input column indices data
int row_data[] = { 0, 2, 3, 4, 6 };
int* row_input = row_data;				// Pointer to input row pointer data

auto* sparse_mat = new CSRMatrix<double>(rows, cols, nnzs, val_input, row_input, col_input);

[...]

delete sparse_mat;
```

#### SPD Declarations:
As this project and its solvers assume that the matrix systems are SPD, there exists a function that automatically generates an SPD of a user-defined range of numbers for which to run the linear solvers. It is called as follows:

```
auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "fixed");		// creates a matrix of fixed near-diagonal values each time it's called
```
Or
```
auto* A_dense = new Matrix<double>(rows, 500, 495, 2, 1, "random");		// creates a matrix of random values each time it's called, with the near-diagonal being larger than the off-diagonal
```
This uses the SPD constructor to declare the new matrix of dimension 'rows', with maximum and minimum diagonal elements of 500 and 495, and maximum and minimum off-diagonal elements of 2 and 1.

All of the solvers (both dense and sparse) obey the following conventions:
```
 void some_solver(some_matrix_class<type>& A, some_matrix_class<type>& b, some_matrix_class<type>& x);

```

The libraries associated with this release are templated classes Matrix<T> and CSRMatrix<T>. Therefore, an example call to a linear solver may be:
```
A->gauss_seidel(*A, *b, *x);
```

Where the inputs (`*A`, `*b`, `*x`) are dereferenced sparse or dense matrix pointers, where x is a generic matrix class pointer whose dereferenced values get overwritten. The syntax has deliberately been written to remind the user that A (left) and b (right) are used to solve for x.

To conserve computational cost and memory storage associated with matrices over simple pointers to arrays, an exception to this convention is made for Conjugate Gradient, which accepts (dereferenced pointers to) arrays for `x` and `b`, which the  Here the syntax is:

```
A->conjgrad(*x, *b, tol, maxiter);
```

Where x and b are pointers to arrays, and maxiter (int) and print (string) arguments specify the maximum number of iterations and whether the matrices are printed, respectively.

**It is always recommended to refer to the header file to check syntax.**

Examples are given in `main.cpp` to demonstrate functionality.

## Authors

* **[Alex Campbell](https://github.com/acse-ac6915)**
* **[Jorge Garcias](https://github.com/acse-jg719)**
* **[Rory Johnston](https://github.com/acse-rej19)**

## Acknowledgements

* Dr. Steven Dargaville - Lecturer in Advanced Programming Methods. Special thanks also to his team of Teaching Assistants for their support and advice.
* Professor Matthew Piggott - ACSE-3 Linear Solvers Course Material from which LU Factorisation and Gaussian Elimination Solvers were adapted.
* Stack Overflow for its eternal fount of knowledge.
