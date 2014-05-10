//Testing/benchmarking functions
//file: testing.hpp
//author: Michael Stumpf

#ifndef _TESTING_ 
#define _TESTING_
	#include <matinc.hpp>
	double get_wall_time();
	double get_cpu_time();
	void calcError(SpMat &A, VectorXd &x, VectorXd &b);
	void compareVec(VectorXd &vec1, VectorXd &vec2);
#endif
