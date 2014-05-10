//Testing/benchmarking functions
//file: testing.hpp
//author: Michael Stumpf

#ifndef _TESTING_ 
#define _TESTING_
	#include <matrixgen.hpp>
	void fillVector(VectorXd &vector);
	void testVector(VectorXd &vec);
	void testArray(double* vec);
	double get_wall_time();
	double get_cpu_time();
	void compareVec(VectorXd &vec1, VectorXd &vec2);
#endif
