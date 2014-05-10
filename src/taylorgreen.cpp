//Functions used for Taylor-Green calculations
//file: taylorgreen.cpp
//author: Michael Stumpf

#include <iostream>
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

void calcTaylorError(VectorXd &u, VectorXd &v){
	VectorXd *u_tg = new VectorXd(NGP*NGP),
			*v_tg = new VectorXd(NGP*NGP);
	calcTaylorGreen(*u_tg,'u',TSMAX*DT);
	calcTaylorGreen(*v_tg,'v',TSMAX*DT);
	double relative_error_u = (u - *u_tg).norm() / v_tg->norm();
	double relative_error_v = (v - *v_tg).norm() / v_tg->norm();
	cout << "Relative error at time step #" << TSMAX << endl
		 << "velocity u: " << relative_error_u << endl
		 << "velocity v: " << relative_error_v << endl;
	delete u_tg;
	delete v_tg;
}
