//General settings
//file: settings.hpp
//author: Michael Stumpf

#ifndef _SETTINGS_ 
#define _SETTINGS_
	const int NGP=300;			//number of grid points in one dimension 
	const double RE=10;			//Reynolds number
	const double LREF=1;		//characteristic length (=lenght of domain)
	const double UREF=1;		//characteristic speed
	const int TSMAX=100;		//number of time steps
	const double DT=0.0000001; 	//time step size
	const double DX=LREF/(NGP-1);	//grid width (grid is uniform)
	//Vortex
	const double A=0.01;			//characteristic radius
	const double LRATIO=30;		//Length ratio b/a
	const double REV=1;		//Reynoldsnumber of vortex
	const double NU=1;			//Viscosity
	//IO settings
	const char OUTDIR[]="/home/michael/Dokumente/KIT/Numerische Strömungsmechanik II/project/data/";
#endif
