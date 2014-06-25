/*
 * settings.hpp
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


//General settings
//file: settings.hpp
//author: Michael Stumpf

#ifndef _SETTINGS_ 
#define _SETTINGS_

//General settings
const char TYPE='v';		//t: Taylor-Green, v: vortices
const int NGP=500;			//number of grid points in one dimension 
const int TSMAX=80000;		//number of time steps
const double DT=0.0000125; 	//time step size
const double LREF=1;		//characteristic length (=lenght of domain)
const double DX=LREF/NGP;	//grid width (grid is uniform)
const double RE=50;			//Reynolds number
//Taylor Green
const double UREF=1;		//characteristic speed
//Vortex
const double A=0.01;		//characteristic radius
const double LRATIO=10;		//Length ratio b/a
const double NU=1;			//Viscosity
const int NUMVOR=2;			//number of vortices
const int THRESHOLD=3;		//threshold for vortex matching (gridpoints per timestep)
//IO settings
const char OUTDIR[]="data/";
const int SAVEINT=250;		//saving interval (field data)
const int PROGSIZE=10;		//print progress # times
	
#endif
