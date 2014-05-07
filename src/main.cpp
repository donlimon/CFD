//CFD project	
//Fractional step method for the incompressible	Navier-Stokes equations
//file: main.cpp
//author: Michael Stumpf

#include <iostream>
#include <matrixgen.hpp>
#include <settings.hpp>

using namespace std;

int main(){
	SpMat mat = createPredCoeffMat();
	cout << mat;
	return 0;
}
