//IO functions
//file: io.hpp
//author: Michael Stumpf

#ifndef _IOFUNC_ 
#define _IOFUNC_
	void printSettings();	//Print settings from settings.hpp on screen
	void printProgress(int ts); //Prints progress on screen
	void writeGridToFile();	//writes gridpoints for u,v to file
	void writeVelocityToFile(VectorXd &u,VectorXd &v,char comp); //writes velocity to ascii file
	void writePressureToFile(VectorXd &phi);	//writes pressure field to ascii file
#endif
