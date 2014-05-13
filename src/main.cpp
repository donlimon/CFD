//CFD project	
//Fractional step method for the incompressible	Navier-Stokes equations
//file: main.cpp
//author: Michael Stumpf

#include <iostream>
#include <Eigen/SparseCholesky> //solver
#include <omp.h>

#include <settings.hpp>
#include <matrixgen.hpp>
#include <io.hpp>
#include <taylorgreen.hpp>

#include <testing.hpp>

using namespace std;
using namespace Eigen;

int main(){
	initParallel();	//Eigen parallel intialization
	printSettings();
	double cpuTime0 = get_cpu_time();
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
			*rhs1 = new VectorXd(NGP*NGP),
			*rhs2 = new VectorXd(NGP*NGP);
	
	cout << "Initializing velocity field..." << endl;
	//Initialize with Taylor-Green at t=0	
	calcTaylorGreen(*u_current,'u',0);
	calcTaylorGreen(*u_previous,'u',0);
	calcTaylorGreen(*v_current,'v',0);
	calcTaylorGreen(*v_previous,'v',0);

<<<<<<< HEAD
	cout << "Initializing solver..." << endl;
	//Initialize solver
	SimplicialLDLT<SparseMatrix<double> > solverMom;
	SimplicialLDLT<SparseMatrix<double> > solverPoi;
	#pragma omp parallel sections
	{
		{
			solverMom.analyzePattern(*coMatMom);
			solverMom.factorize(*coMatMom);
		}
		#pragma omp section
		{
			solverPoi.analyzePattern(*coMatMom);
			solverPoi.factorize(*coMatPoi);
		}
	}
	
	cout << "Starting calculations..." << endl;
	for(int ts=0;ts<TSMAX;ts++){
		printProgress(ts);
		//Solve momentum eq.
		#pragma omp parallel sections
		{
			{
				//Solve for u*
				updateRHSmom(*rhs1, *u_current, *u_previous, 
								*v_current, *v_previous, 'u');
				*u_prelim = solverMom.solve(*rhs1);
			}
			#pragma omp section
			{
				//Solve for v*
				updateRHSmom(*rhs2, *u_current, *u_previous, 
								*v_current, *v_previous, 'v');
				*v_prelim = solverMom.solve(*rhs2);
			}
		}
		//Solve poisson eq.
		updateRHSpoi(*rhs1, *u_prelim, *v_prelim);
		*phi = solverPoi.solve(*rhs1);
=======
	//Solve equations
	SimplicialLDLT<SparseMatrix<double> > solverMom;
	solverMom.analyzePattern(*coMatMom);	
	solverMom.factorize(*coMatMom);
	SimplicialLDLT<SparseMatrix<double> > solverPoi;
	solverPoi.analyzePattern(*coMatPoi);
	solverPoi.factorize(*coMatPoi);
	for(int ts=0;ts<TSMAX;ts++){
		printProgress(ts);
		//Solve momentum eq.
		//Solve for u*
		updateRHSmom(*rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'u');
		*u_prelim = solverMom.solve(*rhs);
		//Solve for v*
		updateRHSmom(*rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'v');
		*v_prelim = solverMom.solve(*rhs);
		//Solve poisson eq.
		updateRHSpoi(*rhs, *u_prelim, *v_prelim);
		*phi = solverPoi.solve(*rhs);
		
		//shift time levels
		delete u_previous;
		delete v_previous;
		u_previous = u_current;
		v_previous = v_current;
		u_current = new VectorXd(NGP*NGP);
		v_current = new VectorXd(NGP*NGP);
>>>>>>> master
		
		#pragma omp parallel sections
		{
			{
				//shift time levels
				delete u_previous;
				u_previous = u_current;
				u_current = new VectorXd(NGP*NGP);
				//solve corrector
				solveCorrector(*u_current,*u_prelim,*phi,'u');
			}
			#pragma omp section
			{
				//shift time levels
				delete v_previous;
				v_previous = v_current;
				v_current = new VectorXd(NGP*NGP);
				//solve corrector
				solveCorrector(*v_current,*v_prelim,*phi,'v');
			}
		}
	}
	calcTaylorError(*u_current,*v_current);
	
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
<<<<<<< HEAD
	delete rhs1;
	delete rhs2;

	cout << "Done!" << endl;
	//char stop; cin >> stop;
=======
	delete rhs;
	cout << "CPU time: " << get_cpu_time()-cpuTime0 << endl;
>>>>>>> master
	return 0;
}
