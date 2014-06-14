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
