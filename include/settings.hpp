//General settings
//file: settings.hpp
//author: Michael Stumpf

#ifndef _SETTINGS_ 
#define _SETTINGS_

//General settings
const char TYPE='t';		//t: Taylor-Green, v: vortices
const int NGP=400;			//number of grid points in one dimension 
const int TSMAX=1000;		//number of time steps
const double DT=0.0001; 	//time step size
const double LREF=1;		//characteristic length (=lenght of domain)
const double DX=LREF/(NGP-1);	//grid width (grid is uniform)
const double RE=10;			//Reynolds number
//Taylor Green
const double UREF=1;		//characteristic speed
//Vortex
const double A=0.01;		//characteristic radius
const double LRATIO=10;		//Length ratio b/a
const double NU=1;			//Viscosity
const int NUMVOR=2;			//number of vortices
const int THRESHOLD=3;		//threshold for vortex matching (gridpoints per timestep)
//IO settings
const char OUTDIR[]="/home/michael/Dokumente/KIT/Numerische Strömungsmechanik II/project/data/";
const int SAVEINT=0;		//saving interval (field data)
const int PROGSIZE=10;		//print progress # times
	
#endif
