/*
 * matrixgen.hpp
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


//Functions for generating/updating coefficient matrices
//file: matrixgen.hpp
//author: Michael Stumpf

#ifndef _MATRIXGEN_ 
#define _MATRIXGEN_
	//General settings for handability
	#include <matinc.hpp>
	
	//Function prototypes
	void createCoeffMat(SpMat &mat,double a, double b);
	void updateRHSmom(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp);
	void updateRHSpoi(VectorXd &rhs, VectorXd &u, VectorXd &v);
	void solveCorrector(VectorXd &velcur, VectorXd &velpre, 
						VectorXd &phi, char comp);
	inline double interpolate(VectorXd &field, int i, int j, char comp);
#endif
