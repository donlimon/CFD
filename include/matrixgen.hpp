//Functions for generating/updating coefficient matrices
//file: matrixgen.hpp
//author: Michael Stumpf

#ifndef _MATRIXGEN_ 
#define _MATRIXGEN_
	//General settings for handability
	#include <Eigen/SparseCore>
	#include <Eigen/Dense>
	typedef Eigen::SparseMatrix<double,Eigen::RowMajor> SpMat;	//Sparse Matrix
	typedef Eigen::Matrix< double,Eigen::Dynamic,1> VectorXd;
	
	//Function prototypes
	void createCoeffMat(SpMat &mat,double a, double b);
	void updateRHS(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp);
	void updateRHSopt(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp);
	inline double interpolate(VectorXd &field, int i, int j, char comp);
#endif
