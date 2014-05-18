//Functions for input/output
//file: io.cpp
//author: Michael Stumpf

#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <matinc.hpp>
#include <io.hpp>
#include <settings.hpp>

using namespace std;
using namespace Eigen;

bool progressInit = false;

void printSettings(){
	cout << "Settings" << endl
		 << "--------------------- " << endl
		 << "Grid Points (in 1D):  " << NGP << endl
		 << "Grid Points (total):  " << NGP*NGP << endl
		 << "Reynolds number:      " << RE << endl
		 << "Domain length:        " << LREF << endl
		 << "Reference speed:      " << UREF << endl
		 << "Grid width:           " << DX << endl
		 << "Time step size:       " << DT << endl
		 << "Number of time steps: " << TSMAX << endl << endl;
}

void printProgress(int ts){
	static double cpuTime0;
	static time_t time0;
	if(progressInit==false){
		cpuTime0=(double)clock()/CLOCKS_PER_SEC;
		time0=time(NULL);
		progressInit=true;
	}
	int interval = TSMAX/10;	//show progress every ~10%
	if((ts+1)%interval==0){
		double cpuTime=(double)clock()/CLOCKS_PER_SEC,
				progress=((double)ts+1)/TSMAX,
				elapTime=cpuTime-cpuTime0;
		//calculate ETA
		time_t eta = time0 + (int)(elapTime/progress);
		struct tm *etas;
		etas = localtime(&eta);
		
		cout << "Time step #" << ts+1 << " of " << TSMAX << endl
			 << "Progress:     " << progress*100 << "%" << endl
			 << "Elapsed time: " << (int)elapTime << " sec" << endl
			 << "Estimated ending time: " << asctime(etas) << endl;		 
	}
}

void writeGridToFile(){
	ofstream xpoints, ypoints;
	string filename;
	
	filename = OUTDIR;
	filename += "grid-x.dat";
	xpoints.open(filename.c_str());
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			xpoints << i*DX << " ";
		}
		xpoints << endl;
	}
	xpoints.close();
	
	filename = OUTDIR;
	filename += "grid-y.dat";
	ypoints.open(filename.c_str());
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			ypoints << j*DX << " ";
		}
		ypoints << endl;
	}
	ypoints.close();
}

void writeVelocityToFile(VectorXd &u,VectorXd &v,char comp){
	ofstream outfile;
	string filename=OUTDIR;
	if(comp=='u'){
		filename += "data-u.dat";
		outfile.open(filename.c_str());
		for(int j=0;j<NGP;j++){
			outfile << (u(j*NGP)+u((NGP-1)+j*NGP))/2 << " ";
			for(int i=1;i<NGP;i++){
				outfile << (u(i+j*NGP)+u((i-1)+j*NGP))/2 << " ";
			}
			outfile << endl;
		}
		outfile.close();
	}
	else if(comp=='v'){
		filename += "data-v.dat";
		outfile.open(filename.c_str());
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				if(j!=0){
					outfile << (v(i+j*NGP)+v(i+(j-1)*NGP))/2 << " ";
				}
				else{
					outfile << (v(i)+v(i+(NGP-1)*NGP))/2 << " ";
				}
			}
			outfile << endl;
		}
		outfile.close();
	}
	else if(comp=='m'){		//magnitude
		filename += "data-velmag.dat";
		double u_inp, v_inp;
		outfile.open(filename.c_str());
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				if(i!=0){
					u_inp = (u(i+j*NGP)+u((i-1)+j*NGP))/2;
				}
				else{
					u_inp = (u(j*NGP)+u((NGP-1)+j*NGP))/2;
				}
				if(j!=0){
					v_inp = (v(i+j*NGP)+v(i+(j-1)*NGP))/2;
				}
				else{
					v_inp = (v(i)+v(i+(NGP-1)*NGP))/2;
				}
				outfile << sqrt(pow(u_inp,2)+pow(v_inp,2)) << " ";
			}
			outfile << endl;
		}
		outfile.close();
	}
}

void writePressureToFile(VectorXd &phi){
	double phi_nab;
	ofstream outfile;
	string filename=OUTDIR;
	filename += "data-p.dat";
	outfile.open(filename.c_str());
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i==0){
				phi_nab = phi((i+1)+j*NGP)+phi((NGP-1)+j*NGP)-4*phi(i+j*NGP);
			}
			else if(i==(NGP-1)){
				phi_nab = phi(0+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);
			}
			else{
				phi_nab = phi((i+1)+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);
			}
			if(j==0){
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(NGP-1)*NGP);
			}
			else if(j==(NGP-1)){
				phi_nab += phi(i+0*NGP)+phi(i+(j-1)*NGP);
			}
			else{
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(j-1)*NGP);
			}
			outfile << phi(i+j*NGP)-(DT/(RE*2))*phi_nab << " ";
		}
		outfile << endl;
	}
	outfile.close();
}
