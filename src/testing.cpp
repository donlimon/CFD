/*
 * testing.cpp
 * This file is part of fracstep
 *
 * Copyright (C) 2014 - Michael Stumpf
 *
 * fracstep is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * fracstep is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with fracstep. If not, see <http://www.gnu.org/licenses/>.
 */

//Functions used for testing
//file: testing.cpp
//author: Michael Stumpf

#include <iostream>
#include <cstdlib>
#include <sys/time.h>
#include <testing.hpp>
#include <settings.hpp>

using namespace std;
using namespace Eigen;

double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

void calcErrorSolve(SpMat &A, VectorXd &x, VectorXd &b){
	double relative_error = (A*x - b).norm() / b.norm(); // norm() is L2 norm
	cout << "The relative error is:\n" << relative_error << endl;
}

void compareVec(VectorXd &vec1, VectorXd &vec2){
	double eps=1e-10;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			//cout << vec1(i+j*NGP) << endl;
			if(abs(vec1(i+j*NGP)-vec2(i+j*NGP)) > eps){
				cout << "Error: i=" << i << " j=" << j << endl;
			}
		}
	}
	cout << "Done comparing!" << endl;
}
