//Functions for generating/updating coefficient matrices
//file: matrixgen.hpp
//author: Michael Stumpf

#ifndef _MATRIXGEN_ 
#define _MATRIXGEN_
	//General settings for handability
	#include <Eigen/SparseCore>
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;
	
	//Function prototypes
	SpMat createPredCoeffMat();
#endif
