/*
 * testing.hpp
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


//Testing/benchmarking functions
//file: testing.hpp
//author: Michael Stumpf

#ifndef _TESTING_ 
#define _TESTING_
	#include <matinc.hpp>
	double get_wall_time();
	double get_cpu_time();
	void calcErrorSolve(SpMat &A, VectorXd &x, VectorXd &b);
	void compareVec(VectorXd &vec1, VectorXd &vec2);
#endif
