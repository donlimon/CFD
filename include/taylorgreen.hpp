//Functions used for Taylor-Green calculations
//file: taylorgreen.hpp
//author: Michael Stumpf

#ifndef _TAYLORGREEN_ 
#define _TAYLORGREEN_
	#include <matinc.hpp>
	void calcTaylorGreen(VectorXd &velocity,char comp,double t);
	void calcTaylorError(VectorXd &u, VectorXd &v);
#endif
