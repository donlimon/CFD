//CFD project	
//Fractional step method for the incompressible	Navier-Stokes equations
//file: main.cpp
//author: Michael Stumpf

#include <iostream>
#include <settings.hpp>
#include <matrixgen.hpp>
#include <io.hpp>
#include <taylorgreen.hpp>

#include <testing.hpp>

using namespace std;

int main(){
	printSettings();
	
	//create coefficient matrix for u,v
	SpMat coeffMatMomentum = 
		createCoeffMat(
			1/DT+2/(RE*DX*DX),
			-1/(2*RE*DX*DX)
		);
	//create coefficient matrix for phi
	SpMat coeffMatPoisson = 
		createCoeffMat(
			1,
			-0.25
		);
		
	//Initialize velocity field
	VectorXd *u_current = new VectorXd(NGP*NGP),
			*u_previous = new VectorXd(NGP*NGP),
			*v_current = new VectorXd(NGP*NGP),
			*v_previous = new VectorXd(NGP*NGP),
			*phi_current = new VectorXd(NGP*NGP),
			*u_rhs = new VectorXd(NGP*NGP),
			*v_rhs = new VectorXd(NGP*NGP),
			*phi_rhs = new VectorXd(NGP*NGP);
	//Initialize with Taylor-Green at t=0	
	calcTaylorGreen(*u_current,'u',0);
	calcTaylorGreen(*u_previous,'u',0);
	calcTaylorGreen(*v_current,'v',0);
	calcTaylorGreen(*v_previous,'v',0);
	
	//Update righthand-side
	updateRHS(*u_rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'u');
	updateRHS(*v_rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'v');	
	cout << "Done!";
	//char stop; cin >> stop;
	//Free heap
	delete u_current;
	delete u_previous;
	delete v_current;
	delete v_previous;
	delete phi_current;
	delete u_rhs;
	delete v_rhs;
	delete phi_rhs;
	
	return 0;
}
