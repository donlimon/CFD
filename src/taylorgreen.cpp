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
	double nu = (LREF*UREF)/RE,
		   x,y;
	//u-component
	if(comp=='u'){
		for(int j=0;j<NGP;j++){
			y=j*DX; //y-coordinate
			for(int i=0;i<NGP;i++){
				x=i*DX + DX/2; //x-coordinate (shifted grid)
				velocity(i+j*NGP) =
					sin(2*M_PI*x)*cos(2*M_PI*y)*exp(-8*M_PI*M_PI*nu*t);
			}
		}
	}
	else if(comp=='v'){
		for(int j=0;j<NGP;j++){
			y=j*DX + DX/2; //y-coordinate (shifted grid)
			for(int i=0;i<NGP;i++){
				x=i*DX; //x-coordinate
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
	double u_min = abs((*u_tg).minCoeff()),
		   u_max = abs((*u_tg).maxCoeff()),
		   v_min = abs((*v_tg).minCoeff()),
		   v_max = abs((*v_tg).maxCoeff()),
		   u_norm = max(u_min,u_max),
		   v_norm = max(v_min,v_max);
	double u_err=0, v_err=0;
	for(int i=0; i<NGP-1; i++){
		for(int j=0; j<NGP-1; j++){
			u_err += pow(u(i+j*NGP) - (*u_tg)(i+j*NGP),2);
			v_err += pow(v(i+j*NGP) - (*v_tg)(i+j*NGP),2);
		}
	}
	u_err = sqrt(u_err)/(NGP*NGP)/u_norm;
	v_err = sqrt(v_err)/(NGP*NGP)/v_norm;
	cout << "Normalized relative error at time step #" << TSMAX << endl
		 << "velocity u: " << u_err << endl
		 << "velocity v: " << v_err << endl;
	cout << "Log(DX)\tLog(Err)" << endl
		 << log(DX) << "\t" << log(u_err) << endl;
	delete u_tg;
	delete v_tg;
}
