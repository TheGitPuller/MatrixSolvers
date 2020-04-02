#include "CSRMatrix.h"

using namespace std;

// Default Constructor for zero input arguments- using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix() : Matrix<T>(0, 0, false), nnzs(0)
{
    // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
    // rows * cols in our base matrix class
    // So then we need to set it to the real value we had passed in
    this->preallocated = false;

    // If we want to handle memory ourselves
    if (this->preallocated)
    {
        // Must remember to delete this in the destructor
        this->values = nullptr;
        this->row_position = nullptr;
        this->col_index = nullptr;
    }
}

// Constructor - using an initialisation list here
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, bool preallocate) : Matrix<T>(rows, cols, false), nnzs(nnzs)
{
    // If we don't pass false in the initialisation list base constructor, it would allocate values to be of size
    // rows * cols in our base matrix class
    // So then we need to set it to the real value we had passed in
    this->preallocated = preallocate;

    // If we want to handle memory ourselves
    if (this->preallocated)
    {
        // Must remember to delete this in the destructor
        this->values = new T[this->nnzs];
        this->row_position = new int[this->rows + 1];
        this->col_index = new int[this->nnzs];
    }
}

// Constructor - now just setting the value of our T pointer
template <class T>
CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index) : Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index)
{}

// Sparsifier - converts dense matrix into sparse
template <class T>
CSRMatrix<T>::CSRMatrix(Matrix<T>& input) : Matrix<T>(input.rows, input.cols, false)
{
    // Calculate how many nnzs by iterating over all the values in the dense input and incrementing
    // the zeros counter by 1 each time we encounter a non-zero number
    this->nnzs = 0;
    for (int i = 0; i < (input.rows * input.cols); i++) {
        if (input.values[i] != 0) {
            this->nnzs++;
        }
    }
    // Declare array pointers on heap based on what we've just learnt about how many nnzs there are.
    this->values = new T[this->nnzs];
    this->row_position = new int[this->rows + 1];
    this->col_index = new int[this->nnzs];

    // Fill in row_position and col_index arrays values based on matrix
    int count = 0;
    this->row_position[0] = 0;
    for (int i = 0; i < this->rows; i++) {
        int nnzs_in_row = 0;    // Reset the row nnzs count to zero at the start of each new row
        for (int j = 0; j < this->cols; j++) {
            // Only do stuff IFF we encounter a non-zero
            if (input.values[i * input.cols + j] != 0) {
                // Fill the values array in order of how they occur in the input matrix.
                this->values[count] = input.values[i * input.cols + j];
                // Fill in the column index at the same rate as the values array gets filled
                // with the column index at this point.
                this->col_index[count] = j;
                // Increment the nnzs in row count by one
                nnzs_in_row++;
                count++;
            }
        }
        // row position is effectively a running count of nnzs we encounter in each row.
        this->row_position[i+1] = this->row_position[i] + nnzs_in_row;    // row fill is always one index ahead
    }
}

// Copy constructor
template<class T>
CSRMatrix<T>::CSRMatrix(const CSRMatrix& old_obj) {
    *this = old_obj;
    // This `new` here is now detaching our `.values` attribute from the `old_obj`'s
    // to a new slab of memory and then writing in our old_obj values to this into
    // this new slab of memory. So the `delete` called during the constructor deletes
    // this slab of memory but keeps the (dangling) pointer.
    this->values = new T[old_obj.nnzs];
    this->col_index = new T[old_obj.nnzs];
    this->row_position = new T[old_obj.rows + 1];
    for (int i = 0; i < old_obj.nnzs; ++i) {
        this->values[i] = old_obj.values[i];
        this->col_index[i] = old_obj.col_index[i];
    }
    for (int j = 0; j < old_obj.rows; j++) {
        this->row_position[j] = old_obj.row_position[j];
    }
}

// destructor
template <class T>
CSRMatrix<T>::~CSRMatrix()
{
    // Delete the values array
    if (this->preallocated) {
        delete[] this->row_position;
        delete[] this->col_index;
    }
    // The super destructor (one that contains pointers to col, row and values), 
    // is called after we finish here. This will delete this->values if preallocated is true
}

// Explicitly print out values
template <class T>
void CSRMatrix<T>::printValues()
{
    for (int j = 0; j < this->nnzs; j++)
    {
        std::cout << this->values[j] << " ";
    }

    std::cout << std::endl;

}

// Explicitly print out matrix in CSR format
template <class T>
void CSRMatrix<T>::printMatrix()
{
    std::cout << "Values: ";
    for (int j = 0; j < this->nnzs; j++)
    {
        std::cout << this->values[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "row_position: ";
    for (int j = 0; j < this->rows + 1; j++)
    {
        std::cout << this->row_position[j] << " ";
    }
    std::cout << std::endl;
    std::cout << "col_index: ";
    for (int j = 0; j < this->nnzs; j++)
    {
        std::cout << this->col_index[j] << " ";
    }
    std::cout << std::endl;
}

// Explicitly print out matrix in dense format
template <class T>
void CSRMatrix<T>::printAsDense()
{
    auto* A = new T[this->rows * this->cols];
    // Initiate depo array to zeros
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            A[i * this->cols + j] = 0.;
        }
    }

    // Fill in A matrix with values from sparse array
    for (int i = 0; i < this->rows; i++) {
        int nnzs_row = this->row_position[i];
        int nnzs_next = this->row_position[i+1];
        for (int j = nnzs_row; j < nnzs_next; j++) {
            A[i * this->cols + this->col_index[j]] = this->values[j];
        }
    }

    // Print Matrix
    for (int i = 0; i < this->rows; i++) {
        for (int j = 0; j < this->cols; j++) {
            cout << A[i * this->cols + j] << " ";
        }
        cout << endl;
    }

    // Delete new pointer to A
    delete[] A;
}

// Do a matrix-vector product
// output = this * input
template<class T>
void CSRMatrix<T>::matVecMult(T* input, T* output)
{
    if (input == nullptr || output == nullptr)
    {
        std::cerr << "Input or output haven't been created" << std::endl;
        return;
    }

    // Set the output to zero
    for (int i = 0; i < this->rows; i++)
    {
        output[i] = 0.0;
    }

    int val_counter = 0;
    // Loop over each row
    for (int i = 0; i < this->rows; i++)
    {
        // Loop over all the entries in this col
        for (int val_index = this->row_position[i]; val_index < this->row_position[i + 1]; val_index++)
        {
            // This is an example of indirect addressing
            // Can make it harder for the compiler to vectorise!
            output[i] += this->values[val_index] * input[this->col_index[val_index]];

        }
    }
}

// Do matrix matrix multiplication
// output = this * mat_right
template <class T>
void CSRMatrix<T>::matMatMult(CSRMatrix<T>& mat_right, CSRMatrix<T>& output)
{
    /* Source code to do matrix-matrix multiplication on two matrices that are stored in CSR (Yale) format.
    This takes is called as follows:
    ` A_mat->matMatMult(*B_mat, *out_mat); `
    OR:
    ` out_mat = *A_mat * *B_mat; ` using overloaded operator *

    Inputs:
    - A_mat, B_mat and out_mat: (dereferenced pointers to) CSR matrix classes.
    (- out_mat can be a a null class (instantiated with a default constructor, as its attributes get overwritten accordingly.)

    Outputs:
    - Void (out_mat overwritten*/

    // Check our dimensions match
    if (this->cols != mat_right.rows)
    {
        std::cerr << "Input dimensions for matrices don't match" << std::endl;
        return;
    }

    // Check if our output matrix has had space allocated to it
    if (output.values != nullptr)
    {
        // Check our dimensions match
        if (this->rows != output.rows || this->cols != output.cols)
        {
            std::cerr << "Input dimensions for matrices don't match" << std::endl;
            return;
        }
    }
    // The output hasn't been preallocated, so we are going to do that
    else
    {
        output.values = new T[this->nnzs];
        // Set to zeros
        for (int i = 0; i < nnzs; i++) {
            output.values[i] = 0.;
        }
    }

    vector<T> row_index_vec;
    vector <tuple<int, int, T>> xyz_v;

    // Loop over all rows of matrix
    for (int i = 0; i < this->rows; i++) {
        // Find for how many non-zeros we'll have to search for per row
        int start = this->row_position[i];
        int end = this->row_position[i + 1];
        // Loop for how many non-zeros there are, from first position to last position per row
        for (int j = start; j < end; j++) {
            // Set the starting point of our column index search
            int col_index_start = this->col_index[j];
            // Iterate start and end positions per row of the right-matrix
            for (int k = mat_right.row_position[col_index_start]; k < mat_right.row_position[col_index_start + 1]; k++) {
                row_index_vec.push_back(i);     // Append the row index of the non-zero produces
                xyz_v.push_back(tuple<int, int, T>(i, mat_right.col_index[k], this->values[j] * mat_right.values[k]));  // Create a tuple of row_pos, col_pos, value of each
                                                                                                                        // newly produced non-zero value (in dense-matrix format)
                
            }
        }
    }

    // Rows are sorted as a result of the algorithm structure.
    // Remove redundant (repeated) entries.
    row_index_vec.erase(unique(row_index_vec.begin(), row_index_vec.end()), row_index_vec.end());

    int beginning_of_this = 0;      // Value that bookmarks the starting index of each row
    int beginning_of_next = 0;      // Value that bookmarks the starting index of the next row
    vector<T> cum_sum;              // Vector to store non-zero values
    vector<T> col_index_vec;        // Vector to store column index values
    vector<T> row_position_count;   // Vector to how many non-zeros pertain to each row.
    // Keep track of how many non-zeros we encounter as we move through entire matrix.
    // This is updated using the nnzs_per_row counter set later.
    int running_total = 0;
    // Set a start flag to indicate whether this is the first loop.
    // This is used at the end of each loop to determine whether the first row_positions need to be padded
    // (in the case of 0 non-zeros on the first n rows)

    bool start = true;

    // Find braces of interest (the indeces of the column that corresponds to a single row)
    // Iterate only through rows that contain values, then fill in the gaps along the way (done at the end of each loop)
    for (auto row_ind : row_index_vec) {
        xyz_v.push_back(tuple<int, int, T>(0, 0, 0.));   // Now add in a pad at the end so we don't reach overflow
        // Find for how many indices in vector our current row spans
        for (int i = beginning_of_this; i < xyz_v.size(); i++) {
            if (get<0>(xyz_v[i]) != get<0>(xyz_v[i + 1])) {     // Find where the row changes
                beginning_of_next = i + 1;                      // Set a bookmark on the next row to indicate where to stop searching in the vectors
                break;                                          // break immediately once we find the next row
            }
        }
        
        // If we reach the end, don't continue anymore
        if (beginning_of_next == beginning_of_this) { break; }
        // Declare a non-zeros counter and set to zero at the beginning of each loop to keep track of how many zeros we encounter
        int nnzs_per_row = 0;
        // Sort all values by column index that are in the same row (defined by "beginning_of..." scopes
        sort(xyz_v.begin() + beginning_of_this, xyz_v.begin() + beginning_of_next, [](const tuple<int, int, T> a, const tuple<int, int, T> b) { return (get<1>(a) < get<1>(b)); });
        // This considers everything in the array between the beginning of this row and the beginning of the next row
        // then shuffles everything in place until it's ordered (increasing positive along the vector) based on the 2nd column (column position)

        // Declare a skip_count to tell the third inner loop where to begin searching from each time so we don't look back over numbers that we've already stored.
        int skip_count = 0;         // start at zero at the beginning of each row
        for (int j_ind = beginning_of_this; j_ind < beginning_of_next; j_ind++) {       // Search from the index that marks the beginning of the row under consideration
                                                                                        // to the index that marks the beginning of the next row
            // We have "pinned" the column vector index j_ind (corresponding to value column index xyz_v[j_ind])
            // and we must loop through all subsequent value column indices in this row, searching for repeated value column indices.
            // Search for duplicate (i, j) pairs and append add all duplicates together
            bool repeat_flag = true;        // Set a flag that tells the conditional gates whether the value we've encountered in the vector is a duplicate or the
                                            // first encounter of it.

            // Go from the index corresponding to the beginning of this line, plus the skip_count, which tells us how far we've already searched
            for (int val_ind = beginning_of_this + skip_count; val_ind < beginning_of_next; val_ind++) {
                if (get<1>(xyz_v[val_ind]) == get<1>(xyz_v[j_ind])) {           // if column index of current item is the same as the pinned value...
                    if (repeat_flag == true) {                                  // > if this is our first encounter of this (i, j)
                        cum_sum.push_back(get<2>(xyz_v[j_ind]));                // >> Append the corresponding value to our cumulative array and...
                        col_index_vec.push_back(get<1>(xyz_v[j_ind]));          // >> ...focus on this (i,j) until there are no more repeats
                        nnzs_per_row++;                                         // >> Increment the non-zeros counter when we encounter the first of a value
                    }
                    else {                                                      // > Otherwise...
                        cum_sum[cum_sum.size() - 1] += get<2>(xyz_v[val_ind]);  // >> Add to the value of interest until there are no more repeats
                    }
                    skip_count++;                                               // > Increment the skip_count by one to tell the outer loop where to start its search from next time
                                                                                // > this reduces how many loops and how many conditional gates we have
                    repeat_flag = false;                                        // > Set repeat flag to false to tell the recurring conditional gates that the values we're searching for
                                                                                //   are duplicates of before
                }                                                                                                                                       
            }                                                                   // Otherwise don't do anything and move on to the next index search.  
        }

        running_total += nnzs_per_row;                                          // Increment the running_total by one to tell the row_position_count vector how many non-zeros we've
                                                                                // encountered so far.

        // Pad with zeros at the start if row index doesn't start at zero
        if (start == true && get<0>(xyz_v[beginning_of_this]) > 0) {            // if this is the start of the matrix and the row_index isn't zero
            for (int i = 0; i < get<0>(xyz_v[beginning_of_this]); i++) {        
                row_position_count.push_back(0);                                // pad the start with zeros for how many rows we're skipping over
            }

        }
        // Then begin filling in the vector with running totals
        if (get<0>(xyz_v[beginning_of_next]) - get<0>(xyz_v[beginning_of_this]) > 0) {      // if we're in the body of the matrix (as given by the jump from this to next lines being > 0               
            for (int i = 0; i < get<0>(xyz_v[beginning_of_next]) - get<0>(xyz_v[beginning_of_this]); i++) {
                row_position_count.push_back(running_total);                                // fill in the values for how many non-zeros we'll encounter between now and the next line that
                                                                                            // contains non-zeros
            }
        } else {                                                                            // otherwise, we must be at the end of the rows that contain non-zeros
            for (int i = 0; i < this->rows - get<0>(xyz_v[beginning_of_this]); i++) {           
                // Find where the row changes
                row_position_count.push_back(running_total);                                // so append the running total to every remaining row that we won't pass through
            }
        }

        beginning_of_this = beginning_of_next;                                              // Reset the start of the next search to be the next line
        start = false;                                                                      // and signal that we are now searching in the body of the matrix
    }

    // Now fill in the values of our output CSR Matrix based on what we've learnt
    // CSRMatrix<T>::CSRMatrix(int rows, int cols, int nnzs, T* values_ptr, int* row_position, int* col_index)
    // : Matrix<T>(rows, cols, values_ptr), nnzs(nnzs), row_position(row_position), col_index(col_index);
    output.nnzs = cum_sum.size();
    output.rows = this->rows;
    output.cols = mat_right.cols;
    output.values = new T[output.nnzs];
    output.col_index = new int[output.nnzs];
    output.row_position = new int[this->rows + 1];      // pad the row_position array with an extra zero
    // these new integers will get deleted when the destructor for this class is called by the compiler.

    // Fill in values and column index arrays
    for (int i = 0; i < cum_sum.size(); i++) {
        output.values[i] = cum_sum[i];
        output.col_index[i] = col_index_vec[i];
    }
    output.row_position[0] = 0; 
    for (int i = 0; i < row_position_count.size(); i++) {
        output.row_position[i + 1] = row_position_count[i];     // miss the first zero of deposit array and fill in values from the row_position_count vector
    }
}

template <class T>
void CSRMatrix<T>::sparse2dense(Matrix<T>& output)
{
    /* Converts current sparse matrix class into a dense matrix, and overwrites these values into the
    given output dense matrix.
    Output is an overwritten matrix class, thus the matrix class can be a null class, instantiated using the default constructor:
    - ` Matrix<T> output()  `. */
    output.rows = this->rows;
    output.cols = this->cols;
    output.values = new T[this->rows * this->cols];


    // Initiate depo array to zeros
    for (int k = 0; k < this->rows; k++) {
        for (int l = 0; l < this->cols; l++) {
            output.values[k * this->cols + l] = 0.;
        }
    }

    // Fill in A matrix with values from sparse array
    for (int i = 0; i < this->rows; i++) {

        int nnzs_row = this->row_position[i];
        int nnzs_next = this->row_position[i + 1];
        for (int j = nnzs_row; j < nnzs_next; j++) {
            output.values[i * this->cols + this->col_index[j]] = this->values[j];
        }
    }
}

template<class T>
CSRMatrix<T>* CSRMatrix<T>::operator* (CSRMatrix<T>& obj) // matrix passed by reference
{
    // Declare a null class to be used to store the outputs of the algorithm and return it to the user.

    CSRMatrix<T>* output_prealloc = new CSRMatrix<T>();
    this->matMatMult(obj, *output_prealloc);

    // This returns a shared pointer to the new matrix
    return output_prealloc;
}

template<class T>
void CSRMatrix<T>::gauss_seidel(CSRMatrix<double>& a, Matrix<double>& b, Matrix<double>& x) {
    //  Calculates 1 iteration of Gauss Seidel on a sparsely stored matrix
    //  Parameters:
    //    Input (for the system Ax = b):
    //      CSRMatrix<double>* a 
    //      Matrix<double>* x_0, b 
    //    Output:
    //      Vector x is updated after a sparseGS iteration
    //  Assumptions:
    //    Ax = b has a unique solution
    //    A has no zeros on its main diagonal
    //    A is diagonally dominant 
    // perform a single iteration of sparse Gauss Seidel 
    int rows = a.rows;
    double tolerance = 10e-10;
    //temporarily store the trace value here
    double trace_val = 0.0;
    int val_counter = 0;

    //store values x_next
    auto* y = new Matrix<double>(a.rows, 1, true); //next iteration
    auto* x_diff = new Matrix<T>(rows, 1, true); // used to compute tolerance
    T norm = 1.0;
    double omega = 1;
    auto* x_prev = new Matrix<T>(rows, 1, true);
    // Iterate until convergence 
    while (norm > tolerance) {
        for (int i = 0; i < x.rows; i++) {
            x_prev->values[i] = x.values[i];
        }
        norm = 0;
        // Loop over rows
        for (int i = 0; i < rows; i++) {
            // Find trace element a[i][i]:
            for (int val_index = a.row_position[i]; val_index < a.row_position[i + 1]; val_index++) {
                if (a.col_index[val_index] == i) {
                    trace_val = a.values[val_index];
                }
            }
            // Compute n[i] = b[i]/a[i][i]
            y->values[i] = b.values[i] / trace_val;
            // Loop over non-zero columns of A
            for (int val_index = a.row_position[i]; val_index < a.row_position[i + 1]; val_index++) {
                // Compute n[i] = n[i] - a[i][j]/a[i][i]*x[j]
                if (a.col_index[val_index] == i) continue;
                y->values[i] = y->values[i] - a.values[val_index] * x.values[a.col_index[val_index]] / trace_val;
                x.values[i] = y->values[i] + (1 - omega) * x.values[i];
            }
        }

        // Calculate 2-norm of residual
        for (int i = 0; i < y->size_of_values; ++i) {
            x_prev->values[i] = x.values[i] - x_prev->values[i];
            norm += pow(x_prev->values[i], 2);
        }
        norm = sqrt(norm / y->rows);
    }
    delete y;
    delete x_diff;
    delete x_prev;
}
