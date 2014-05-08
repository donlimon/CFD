//Functions for input/output
//file: io.cpp
//author: Michael Stumpf

#include <iostream>
#include <settings.hpp>

using namespace std;

void printSettings(){
	cout << "> Settings" << endl
		 << "---------------------" << endl
		 << "Grid Points (in 1D): " << NGP << endl
		 << "Grid Points (total): " << NGP*NGP << endl
		 << "Grid width:          " << DX << endl
		 << "Time step:           " << DT << endl
		 << "Reynolds number:     " << RE << endl;
}
