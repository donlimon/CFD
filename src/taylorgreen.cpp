//Functions used for Taylor-Green calculations
//file: taylorgreen.cpp
//author: Michael Stumpf

#include <cmath>
#include <taylorgreen.hpp>
#include <settings.hpp>

using namespace std;
using namespace Eigen;

void calcTaylorGreen(VectorXd &velocity,char comp,double t){
	double nu = (LREF*UREF)/RE;
	//u-component
	if(comp=='u'){
		for(int j=0;j<NGP;j++){
			double y=((double)j)/(NGP-1); //y-coordinate
			for(int i=0;i<NGP;i++){
				double x=((double)i)/(NGP-1); //x-coordinate
				velocity(i+j*NGP) =
					sin(2*M_PI*x)*cos(2*M_PI*y)*exp(-8*M_PI*M_PI*nu*t);
			}
		}
	}
	else if(comp=='v'){
		for(int j=0;j<NGP;j++){
			double y=((double)j)/(NGP-1); //y-coordinate
			for(int i=0;i<NGP;i++){
				double x=((double)i)/(NGP-1); //x-coordinate
				velocity(i+j*NGP) =
					-sin(2*M_PI*y)*cos(2*M_PI*x)*exp(-8*M_PI*M_PI*nu*t);
			}
		}
	}
}
