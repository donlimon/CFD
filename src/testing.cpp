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

void fillVector(VectorXd &vec){
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			vec(i+(j*NGP)) = 10;
		}
	}
}

void testVector(VectorXd &vec){
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			vec(i+(j*NGP)) = rand();
		}
	}
	//time measurement
	double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    for(int j=1;j<NGP-1;j++){
		for(int i=1;i<NGP-1;i++){
			vec(i+(j*NGP)) = vec(i-1+(j*NGP)) + vec(i+1+(j*NGP))
							+vec(i+(j+1)*NGP) + vec(i+(j-1)*NGP);
		}
	}
	double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
	cout << "Eigen" << endl;
    cout << "Wall Time = " << wall1 - wall0 << endl;
    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
}

void testArray(double* vec){
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			vec[i+(j*NGP)] = rand();
		}
	}
	//time measurement
	double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
    for(int j=1;j<NGP-1;j++){
		for(int i=1;i<NGP-1;i++){
			vec[i+(j*NGP)] = vec[i-1+(j*NGP)] + vec[i+1+(j*NGP)]
							+vec[i+(j+1)*NGP] + vec[i+(j-1)*NGP];
		}
	}
	double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
	cout << "Array" << endl;
    cout << "Wall Time = " << wall1 - wall0 << endl;
    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
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
