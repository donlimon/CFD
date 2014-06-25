/*
 * io.hpp
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

//IO functions
//file: io.hpp
//author: Michael Stumpf

#ifndef _IOFUNC_ 
#define _IOFUNC_

#include <structs.hpp>

void printSettings();	//Print settings from settings.hpp on screen
void printProgress(int ts); //Prints progress on screen
void checkStabilityCriteria(VectorXd &u, VectorXd &v);
//ASCII file functions (for debugging)
void writeGridToFile();	//writes gridpoints for u,v to file
void writeVelocityToFile(VectorXd &u,VectorXd &v); //writes velocity to ascii file
void writePressureToFile(VectorXd &phi);	//writes pressure field to ascii file
void writePhiToFile(VectorXd &phi);			//write phi to ascii file
void writeVorticityToFile(VectorXd &omega);		//write vorticity to ascii file
//binary file functions (for data export to matlab)
void writeGridToBinary();
void writeVelocityToBinary(VectorXd &u,VectorXd &v,int ts);
void writePressureToBinary(VectorXd &phi, int ts);
void writeVorticityToBinary(VectorXd &omega, int ts);
void writeInfoToBinary(int ts);
void writeVortexLocationToBinary(location *locData, int ts);
void saveData(VectorXd &u, VectorXd &v, VectorXd &phi, VectorXd &omega, int ts);
	
#endif
