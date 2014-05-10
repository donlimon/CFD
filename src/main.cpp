//CFD project	
//Fractional step method for the incompressible	Navier-Stokes equations
//file: main.cpp
//author: Michael Stumpf

#include <iostream>
#include <Eigen/SparseCholesky> //solver
#include <settings.hpp>
#include <matrixgen.hpp>
#include <io.hpp>
#include <taylorgreen.hpp>

#include <testing.hpp>

using namespace std;
using namespace Eigen;

int main(){
	printSettings();
	
	//coefficient matrix for momentum/poisson eq.
	SpMat* coMatMom = new SpMat(NGP*NGP,NGP*NGP);
	SpMat* coMatPoi = new SpMat(NGP*NGP,NGP*NGP);
	createCoeffMat(*coMatMom,1/DT+2/(RE*DX*DX),-1/(2*RE*DX*DX));
	createCoeffMat(*coMatPoi,1,-0.25);
		
	//Allocate memory
	VectorXd *u_current = new VectorXd(NGP*NGP),
			*u_previous = new VectorXd(NGP*NGP),
			*u_prelim = new VectorXd(NGP*NGP),
			*v_current = new VectorXd(NGP*NGP),
			*v_previous = new VectorXd(NGP*NGP),
			*v_prelim = new VectorXd(NGP*NGP),
			*phi = new VectorXd(NGP*NGP),
			*rhs = new VectorXd(NGP*NGP);
			
	//Initialize with Taylor-Green at t=0	
	calcTaylorGreen(*u_current,'u',0);
	calcTaylorGreen(*u_previous,'u',0);
	calcTaylorGreen(*v_current,'v',0);
	calcTaylorGreen(*v_previous,'v',0);
		
	//Update righthand-side

	//Solve equations
	SimplicialLDLT<SparseMatrix<double> > solver;
	solver.analyzePattern(*coMatMom);	//Pattern is equal for poisson and momentum eq.
	for(int ts=0;ts<TSMAX;ts++){
		//Solve momentum eq.
		solver.factorize(*coMatMom);
		//Solve for u*
		updateRHSmom(*rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'u');
		*u_prelim = solver.solve(*rhs);
		//calcError(*coMatMom,*u_prelim,*rhs);
		//Solve for v*
		updateRHSmom(*rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'v');
		*v_prelim = solver.solve(*rhs);
		//calcError(*coMatMom,*v_prelim,*rhs);
		//Solve poisson eq.
		solver.factorize(*coMatPoi);
		updateRHSpoi(*rhs, *u_prelim, *v_prelim);
		*phi = solver.solve(*rhs);
		//calcError(*coMatPoi,*phi,*rhs);
	}

	cout << "Done!";
	//char stop; cin >> stop;
	
	//Free heap
	delete coMatMom;
	delete coMatPoi;
	delete u_current;
	delete u_previous;
	delete u_prelim;
	delete v_current;
	delete v_previous;
	delete v_prelim;
	delete phi;
	delete rhs;
	
	return 0;
}
