//Functions used for testing
//file: testing.cpp
//author: Michael Stumpf

#include <testing.hpp>
#include <settings.hpp>

using namespace std;
using namespace Eigen;

void fillVector(VectorXd &vec){
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			vec(i+(j*NGP)) = 10;
		}
	}
}
