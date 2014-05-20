//IO functions
//file: io.hpp
//author: Michael Stumpf

#ifndef _IOFUNC_ 
#define _IOFUNC_
	void printSettings();	//Print settings from settings.hpp on screen
	void printProgress(int ts); //Prints progress on screen
	void checkCFL(VectorXd &u, VectorXd &v, bool recommendTS);
	void writeGridToFile();	//writes gridpoints for u,v to file
	void writeVelocityToFile(VectorXd &u,VectorXd &v); //writes velocity to ascii file
	void writePressureToFile(VectorXd &phi);	//writes pressure field to ascii file
	void writePhiToFile(VectorXd &phi);			//write phi to ascii file
	void writeGridToBinary();
	void writeVelocityToBinary(VectorXd &u,VectorXd &v,int ts);
	void writePressureToBinary(VectorXd &phi, int ts);
	void writeInfoToBinary(int ts);
	void saveData(VectorXd &u, VectorXd &v, VectorXd &phi, int ts);
#endif
