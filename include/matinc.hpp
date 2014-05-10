//Necessary includes for Eigen classes and typedefs
//file: matinc.hpp
//author: Michael Stumpf

#ifndef _MATINC_ 
#define _MATINC_
	#include <Eigen/SparseCore>
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;	//Sparse Matrix
	typedef Eigen::Matrix< double,Eigen::Dynamic,1> VectorXd;	//Vector
#endif
