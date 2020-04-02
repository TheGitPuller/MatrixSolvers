#ifndef CSRM_H
#define CSRM_H

#include "Matrix.h"

using namespace std;

template <class T>
class CSRMatrix : public Matrix<T>
{
public:
    // Default Constructor for zero input arguments- using an initialisation list here
    CSRMatrix();
    // constructor where we want to preallocate ourselves
    CSRMatrix(int rows, int cols, int nnzs, bool preallocate);
    // constructor where we already have allocated memory outside
    CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index);
    // Copy constructor
    CSRMatrix(const CSRMatrix& old_obj);
    // Construct a sparse matrix class out of a given dense matrix class
    CSRMatrix(Matrix<T>& output);

    // destructor
    ~CSRMatrix();

    // Print out the values in our matrix
    virtual void printValues();

    // Print out the matrix info in CSR format
    virtual void printMatrix();

    // Print out the matrix in dense format
    void printAsDense();

    // Perform some operations with our matrix
    void matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output);
    // Perform some operations with our matrix
    void matVecMult(T* input, T* output);

    // Convert sparse matrix to dense matrix
    void sparse2dense(Matrix<T>& output);

    // Gauss-Seidel solver 
    virtual void gauss_seidel(CSRMatrix<double>& a, Matrix<double>& b, Matrix<double>& x);

    // This creates a new slab of memory on the heap which the user is now responsible for.
    virtual CSRMatrix<T>* operator* (CSRMatrix<T>& obj);


    // Explicitly using the C++11 nullptr here
    int* row_position = nullptr;
    int* col_index = nullptr;

    // How many non-zero entries we have in the matrix
    int nnzs = -1;
};

#endif // !CSRM_H