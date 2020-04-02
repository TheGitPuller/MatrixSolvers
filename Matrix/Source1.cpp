//#include <iostream>
//#include <math.h>
//#include <ctime>
//#include <string>
//#include "Matrix.h"
//#include "Matrix.cpp"	// mandatory if using TEMPLATING
//#include "CSRMatrix.h"
//#include "CSRMatrix.cpp"
//#include <stdio.h>
//using namespace std;
//
//
//
//int main()
//{
//
//	const int rows = 2000; // altering rows and cols with each loop 
//	const int cols = 2000;
//
//	auto* A = new Matrix<double>(rows, rows, rows, 2, 1, "fixed"); // CALLING SPD MATRIX CONSTRUCTOR max/min diag and off-diag specified. Set a fixed matrix
//	auto* B = new Matrix<double>(rows, rows, rows, 2, 1, "fixed");
//	auto* output = new Matrix<double>(rows, rows, true);
//
//	A->printMatrix();
//	B->printMatrix();
//	A->matMatMult(*B, *output); // outputting our exac5 b vector from a known x and A
//
//
//	output->printMatrix();
//	delete A;
//	delete B;
//	delete output;
//	return 0;
//}
