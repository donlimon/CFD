//General settings
//file: settings.hpp
//author: Michael Stumpf

#ifndef _SETTINGS_ 
#define _SETTINGS_

//General settings
const char TYPE='v';		//t: Taylor-Green, v: vortices
const int NGP=500;			//number of grid points in one dimension 
const int TSMAX=5000;		//number of time steps
const double DT=0.000002; 	//time step size
const double LREF=1;		//characteristic length (=lenght of domain)
const double DX=LREF/(NGP-1);	//grid width (grid is uniform)
//Taylor Green
const double RE=10;			//Reynolds number
const double UREF=1;		//characteristic speed
//Vortex
const double A=0.01;		//characteristic radius
const double LRATIO=10;		//Length ratio b/a
const double REV=10;		//Reynoldsnumber of vortex
const double NU=1;			//Viscosity
const int NUMVOR=2;			//number of vortices
const int THRESHOLD=3;		//threshold for vortex matching (gridpoints per timestep)
//IO settings
const char OUTDIR[]="/home/michael/Dokumente/KIT/Numerische Strömungsmechanik II/project/data/";
const int SAVEINT=50;		//saving interval (field data)
const int PROGSIZE=100;		//print progress # times
	
#endif
