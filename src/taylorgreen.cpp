/*
 * taylorgreen.cpp
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

//Functions used for Taylor-Green calculations
//file: taylorgreen.cpp
//author: Michael Stumpf

#include <iostream>
#include <cmath>
#include <taylorgreen.hpp>
#include <settings.hpp>
#include <io.hpp>

using namespace std;
using namespace Eigen;

/* FUNCTIONS PROVIDED FOR OTHER CPP */

void calcTaylorGreen(VectorXd &velocity,char comp,double t){
	double nu = (LREF*UREF)/RE;
	//u-component
	if(comp=='u'){
		for(int j=0;j<NGP;j++){
			double y=((double)j)/NGP; //y-coordinate
			for(int i=0;i<NGP;i++){
				double x=((double)i)/NGP + DX/2; //x-coordinate (shifted grid)
				velocity(i+j*NGP) =
					sin(2*M_PI*x)*cos(2*M_PI*y)*exp(-8*M_PI*M_PI*nu*t);
			}
		}
	}
	else if(comp=='v'){
		for(int j=0;j<NGP;j++){
			double y=((double)j)/NGP + DX/2; //y-coordinate (shifted grid)
			for(int i=0;i<NGP;i++){
				double x=((double)i)/NGP; //x-coordinate
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
	u_err = sqrt(u_err/(NGP*NGP))/u_norm;
	v_err = sqrt(v_err/(NGP*NGP))/v_norm;
	cout << "Normalized relative error at time step #" << TSMAX << endl
		 << "velocity u: " << u_err << endl
		 << "velocity v: " << v_err << endl;
	cout << "DX\tErr" << endl
		 << DX << "\t" << u_err << endl;
	cout << "DT\tErr" << endl
		 << DT << "\t" << u_err << endl;
		 
	writeVelocityToBinary(u,v,1);
	writeVelocityToBinary(*u_tg,*v_tg,2);
	delete u_tg;
	delete v_tg;
}
