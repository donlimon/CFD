//Functions for generating/updating coefficient matrices
//file: matrixgen.hpp
//author: Michael Stumpf

#ifndef _MATRIXGEN_ 
#define _MATRIXGEN_
	//General settings for handability
	#include <matinc.hpp>
	
	//Function prototypes
	void createCoeffMat(SpMat &mat,double a, double b);
	void updateRHSmom(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp);
	void updateRHSpoi(VectorXd &rhs, VectorXd &u, VectorXd &v);
	inline double interpolate(VectorXd &field, int i, int j, char comp);
#endif
