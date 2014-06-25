/*
 * main.cpp
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

//CFD project	
//Fractional step method for the incompressible	Navier-Stokes equations
//file: main.cpp
//author: Michael Stumpf

#include <iostream>
#include <Eigen/SparseCholesky>
#include <omp.h>

#include <settings.hpp>
#include <matrixgen.hpp>
#include <io.hpp>
#include <taylorgreen.hpp>
#include <vortex.hpp>
#include <structs.hpp>

using namespace std;
using namespace Eigen;

/* LOCAL PROTOTYPES */
inline void checkDivergence(VectorXd &u, VectorXd &v);

/* MAIN FUNCTION */
int main(){
	initParallel();	//Eigen parallel intialization
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
			*rhs = new VectorXd(NGP*NGP),
			*omega = new VectorXd(NGP*NGP);
	location *vorPos = new location[NUMVOR];
	
	cout << "Initializing velocity field..." << endl;
	if(TYPE=='t'){
		//Initialize with Taylor-Green at t=0	
		calcTaylorGreen(*u_current,'u',0);
		calcTaylorGreen(*u_previous,'u',0);
		calcTaylorGreen(*v_current,'v',0);
		calcTaylorGreen(*v_previous,'v',0);
	}
	else if(TYPE=='v'){
		//Initialize with vortices
		initializeVortex(*u_current,*v_current);
		initializeVortex(*u_previous,*v_previous);
	}
	
	//Check CFL condition at t=0
	checkStabilityCriteria(*u_current,*v_current);
	
	//Save settings, grid data and initial velocity field
	writeInfoToBinary(0);
	writeGridToBinary();
	if(TYPE=='v')
		writeVelocityToBinary(*u_current,*v_current,0);

	//Initialize solver
	cout << "Initializing solver..." << endl;
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
	
	//Calculate solutions
	cout << "Starting calculations..." << endl;
	for(int ts=0;ts<TSMAX;ts++){
		printProgress(ts);
		//Solve momentum eq.
		updateRHSmom(*rhs, *u_current, *u_previous, 
						*v_current, *v_previous, 'u');
		*u_prelim = solverMom.solve(*rhs);
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
		
		//solve corrector
		solveCorrector(*u_current,*u_prelim,*phi,'u');
		solveCorrector(*v_current,*v_prelim,*phi,'v');
		
		//Check if solution diverges
		checkDivergence(*u_current,*v_current);
		
		//Calculate vorticity
		if(TYPE=='v')
			calcVorticity(*omega,*u_current,*v_current);
		//Locate vortex centers and save to binary file
		if(TYPE=='v'){
			locateVortex2(vorPos,*omega,ts+1);
			writeVortexLocationToBinary(vorPos,ts+1);	
		}
		//save field data for visual post processing
		if(TYPE=='v' && ((ts+1)%SAVEINT==0 || ts==TSMAX-1))
			saveData(*u_current,*v_current,*phi,*omega,ts+1);
	}
	//Check error if Taylor-Green
	if(TYPE=='t')
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
	delete rhs;
	delete omega;
	delete vorPos;

	cout << "Done!" << endl;

	return 0;
}

/* LOCAL FUNCTIONS (ONLY ACCESSIBLE IN THIS CPP) */
inline void checkDivergence(VectorXd &u, VectorXd &v){
	if(u.hasNaN()){
		cout << "ERROR: Divergence detected!" << endl;
		exit(1);
	}
	else if(v.hasNaN()){
		cout << "ERROR: Divergence detected!" << endl;
		exit(1);
	}
}
