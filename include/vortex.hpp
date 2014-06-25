/*
 * vortex.hpp
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


//Implementation of two counterrotating vortices
//file: vortex.hpp
//author: Michael Stumpf

#ifndef _VORTEX_
#define _VORTEX_

#include <matinc.hpp>
#include <structs.hpp>

void initializeVortex(VectorXd &u, VectorXd &v);
void calcVorticity(VectorXd &omega, VectorXd &u, VectorXd &v);
void locateVortex(location *newPos, VectorXd &phi, int ts);	//uses candidates and threshold
void locateVortex2(location *newPos, VectorXd &omega, int ts); //uses global min/max of vorticity

#endif
