//#include "Matrix.h"
//#include "Matrix.cpp"	// mandatory if using TEMPLATING
//#include "CSRMatrix.h"
//#include "CSRMatrix.cpp"
//
//using namespace std;
//int main()
//{
//	//int rows = 4;
//	//int cols = 4;
//
//	//double A_data[] = { 1., 0., 3., 7., 2., 1., 0., 4., 5., 4., 1., -2., 4., 1., 6., 2. };
//	//double* A_vals = A_data;
//	//auto* dense_mat = new Matrix<double>(rows, rows, A_vals);
//
//	//double b_data[] = { 1., 2., -3., 2. };
//	//double* b_vals = b_data;
//	//auto* b = new Matrix<double>(rows, 1, b_vals);
//
//	//auto* output = new Matrix<double>(rows, 1, true);
//
//	//std::cout << "Doing matrix multiplication between: \n";
//	//std::cout << "A: \n";
//	//dense_mat->printMatrix();
//	//std::cout << "and B: \n";
//	//b->printMatrix();
//	//// EXAMPLE FOR HOW TO RUN CODE NORMALLY (WITHOUT OVERLOADING)
//	//dense_mat->matMatMult(*b, *output);		// <==> *(dense_mat).matMatmult(...)
//	//std::cout << "Result: \n";
//	//output->printMatrix();
//
//	//std::cout << "\nOverwriting the * operator to do matrix multiplication between: \n";
//	//std::cout << "A: \n";
//	//dense_mat->printMatrix();
//	//std::cout << "and B: \n";
//	//b->printMatrix();
//	//// EXAMPLE FOR HOW TO RUN CODE USING OVERLOADING
//	//// You are now responsible for this new chunk of data - DOCUMENT THIS!!!
//	//output = *dense_mat * *b;
//	//std::cout << "Result: \n";
//	//output->printMatrix();
//
//	//std::cout << "Original Dense Matrix" << endl;
//	//dense_mat->printMatrix();
//	//output = dense_mat->transpose();
//	//std::cout << "Transposed Dense Matrix" << endl;
//	//output->printMatrix();
//	//delete output;
//	//output = nullptr;
//
//	//std::cout << "Print identity: \n";
//	//auto* identity = new Matrix<double>(5, "eye");
//	//identity->printMatrix();
//	//delete identity;
//
//	//std::cout << "Print zeros: \n";
//	//auto* zeros = new Matrix<double>(5, "zeros");
//	//zeros->printMatrix();
//	//delete zeros;
//
//	//std::cout << "Print diagonal: \n";
//	//auto* diag = new Matrix<double>(b_vals, 4, "diag");
//	//diag->printMatrix();
//	//delete diag;
//
//	//output = new Matrix<double>(rows, 1, true);
//
//	//std::cout << "\nSolving Ax = b system for x: \n";
//
//	//std::cout << "GAUSSIAN ELIMINATION: \n";
//	////DOCUMENT THAT THE OUTPUT NEEDS TO BE DELETED
//	//dense_mat->gaussian_elimination_solve(*dense_mat, *b, *output);
//	//std::cout << "GE solution for x: \n";
//	//output->printMatrix();
//
//	//output = new Matrix<double>(rows, 1, true);
//
//	//std::cout << "LU DECOMPOSITION: \n";
//	////DOCUMENT THAT THE OUTPUT NEEDS TO BE DELETED
//	//dense_mat->LU_solve(*dense_mat, *b, *output);
//	//std::cout << "LU solution for x: \n";
//	//output->printMatrix();
//
//	///////////////
//	//cout << endl;
//	//cout << endl;
//	//cout << "Dense Matrix to be turned sparse: \n";
//	//dense_mat->printMatrix();
//	//cout << endl;
//	//cout << endl;
//	//auto* sparsified = new CSRMatrix<double>(*dense_mat);
//	//cout << "New matrix: \n";
//	//sparsified->printMatrix();
//	//cout << endl;
//	//cout << endl;
//	//cout << "Sparse matrix, printed as dense: \n";
//	//sparsified->printAsDense();
//	//cout << endl;
//	//cout << endl;
//	////cout << "Testing Alex's code" << endl << endl << endl << endl << endl;
//	//auto* A = new Matrix<double>(rows, 500, 490, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified
//	//A->test_accuracy("dense");
//	///////////////
//
//	//// Delete created classes
//	//delete dense_mat;
//	//delete b;
//	//delete output;
//	//delete sparsified;
//	//delete A;
//
//	//// TESTING SPARSE MATRICES
//
//	// Define some sparse class matrices
//	int rows_sparse = 4;
//	int cols_sparse = 4;
//	int nnzs = 6;
//	
//	double val_data[] = { 1., 2., 5., 3., 1., 3. };
//	double* val_input = val_data;
//	int col_data[] = { 2, 3, 1, 2, 0, 1 };
//	int* col_input = col_data;
//	int row_data[] = { 0, 2, 3, 4, 6 };
//	int* row_input = row_data;
//
//	// Define sparse class
//	auto* sparse_mat = new CSRMatrix<double>(rows_sparse, cols_sparse, nnzs, val_input, row_input, col_input);
//	std::cout << "Print Sparse Values: \n";
//	sparse_mat->printValues();
//	std::cout << endl;
//
//	std::cout << "Print Sparse Matrix:\n";
//	sparse_mat->printMatrix();
//	std::cout << endl;
//
//	std::cout << "Print Sparse Matrix as dense:\n";
//	sparse_mat->printAsDense();
//	std::cout << endl;
//
//	auto* output_to_dense = new Matrix<double>();
//	sparse_mat->sparse2dense(*output_to_dense);
//	std::cout << "Print converted Matrix:\n";
//	output_to_dense->printMatrix();
//	std::cout << endl;
//	delete output_to_dense;
//
//	double* output_matVec = new double[rows_sparse];
//	double input_vals_matVec[] = { 1., 3., 7., 2., };
//	double* input_matVec = input_vals_matVec;
//	std::cout << "Doing matVecMult\n";
//	sparse_mat->matVecMult(input_matVec, output_matVec);
//	std::cout << "Printing sparse_mat_poly: \n";
//	for (int i = 0; i < rows_sparse; i++) {
//		std::cout << output_matVec[i] << '\n';
//	}
//	std::cout << endl;
//	delete[] output_matVec;
//	
//	//double A_val_data[] = { 1., 2., 3., 4., 5. };
//	//double* A_val_input = A_val_data;
//	//int A_col_data[] = { 0, 2, 3, 2, 1 };
//	//int* A_col_input = A_col_data;
//	//int A_row_data[] = { 0, 3, 3, 4, 5};
//	//int* A_row_input = A_row_data;
//
//	double A_val_data[] = { 5., 8., 3., 6. };
//	double* A_val_input = A_val_data;
//	int A_col_data[] = { 0, 1, 2, 1 };
//	int* A_col_input = A_col_data;
//	int A_row_data[] = { 0, 0, 2, 3, 4};
//	int* A_row_input = A_row_data;
//
//	// Define sparse class
//	CSRMatrix<double>* A_mat = new CSRMatrix<double>(4, 5, 5, A_val_input, A_row_input, A_col_input);
//	A_mat->printAsDense();
//
//	double B_val_data[] = { 10., 20., 30., 40., 50., 60., 70., 80.};
//	double* B_val_input = B_val_data;
//	int B_col_data[] = { 0, 1, 1, 3, 2, 3, 4, 5 };
//	int* B_col_input = B_col_data;
//	int B_row_data[] = { 0, 2, 4, 7, 8, 8};
//	int* B_row_input = B_row_data;
//
//
//	// Define sparse class
//	CSRMatrix<double>* B_mat = new CSRMatrix<double>(5, 6, 8, B_val_input, B_row_input, B_col_input);
//	//B_mat->printAsDense();
//	CSRMatrix<double>* out_mat = new CSRMatrix<double>();
//	cout << "Doing mat-mat mult\n";
//	cout << "Between A = \n";
//	A_mat->printAsDense();
//	cout << "and b = \n";
//	B_mat->printAsDense();
//	//A_mat->matMatMult(*B_mat, *out_mat);
//	cout << "...using the overwritten * operator" << endl;
//	out_mat = *A_mat * *B_mat;
//
//	cout << "Printing output matrix: \n";
//	out_mat->printMatrix();
//	cout << endl;
//
//	cout << "Printing output matrix as dense: \n";
//	out_mat->printAsDense();
//	cout << endl;
//
//	delete sparse_mat;
//	delete A_mat;
//	delete B_mat;
//	delete out_mat;
//
//	return(0);
//
//}
//
