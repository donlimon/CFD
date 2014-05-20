//Functions for input/output
//file: io.cpp
//author: Michael Stumpf

#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <sstream>
#include <matinc.hpp>
#include <io.hpp>
#include <settings.hpp>

using namespace std;
using namespace Eigen;

//global variables
bool progressInit = false;
//prototypes for local functions
string getFilenameDigits(int);

void printSettings(){
	cout << "Settings" << endl
		 << "---------------------- " << endl;
	cout << "Simulation type:       ";
	if(TYPE=='t')
		cout << "Taylor-Green" << endl;
	else if(TYPE=='v')
		cout << "Vortex" << endl;
	else
		cout << "ERROR" << endl;
		
	cout << "Grid Points (in 1D):   " << NGP << endl
		 << "Grid Points (total):   " << NGP*NGP << endl
		 << "Grid width:            " << DX << endl
		 << "Time step size:        " << DT << endl
		 << "Number of time steps:  " << TSMAX << endl
		 << "Domain length:         " << LREF << endl;
	if(TYPE=='t'){
		cout << "Reynolds number:       " << RE << endl
		 	 << "Reference speed:       " << UREF << endl << endl;
	}
	else if(TYPE=='v'){
		cout << "Reynolds number:       " << REV << endl
		 	 << "Characteristic radius: " << A << endl
		 	 << "Length ratio:          " << LRATIO << endl << endl;
	}
	else{
		cout << "Aborting..." << endl;
		exit(2);
	}
}

void printProgress(int ts){
	static time_t time0;
	if(progressInit==false){
		time0=time(NULL);
		progressInit=true;
	}
	int interval = TSMAX/10;	//show progress every ~10%
	if(interval==0){
		return;
	}
	else if((ts+1)%interval==0){
		time_t curTime;
		double progress, elapTime;
		
		curTime = time(NULL);
		progress = ((double)ts+1)/TSMAX,
		elapTime = difftime(curTime,time0);
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

void checkCFL(VectorXd &u, VectorXd &v,bool recommendTS){
	double vmag_max = 0, vmag_tmp, u_inp, v_inp, CFL;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i!=0)
				u_inp = (u(i+j*NGP)+u((i-1)+j*NGP))/2;
			else
				u_inp = (u(j*NGP)+u((NGP-1)+j*NGP))/2;
			if(j!=0)
				v_inp = (v(i+j*NGP)+v(i+(j-1)*NGP))/2;
			else
				v_inp = (v(i)+v(i+(NGP-1)*NGP))/2;
				
			vmag_tmp = sqrt(pow(u_inp,2)+pow(v_inp,2));
			if(vmag_tmp > vmag_max)
				vmag_max = vmag_tmp;
		}
	}
	CFL = DT*vmag_max/DX;
	cout << "Current CFL number: " << CFL << endl;
	if(CFL>1)
		cout << "WARNING: simulation might be unstable (CFL > 1)" << endl;
	if(recommendTS){
		cout << "Recommended time step size: "  << endl
		 	 << "CFL = 0.5: " << 0.5*DX/vmag_max << endl
		 	 << "CFL = 0.1: " << 0.1*DX/vmag_max << endl; 
	}
}

void writeGridToFile(){
	string path = OUTDIR;
	string xFile = path + "grid-x.dat";
	string yFile = path + "grid-y.dat";
	
	ofstream xData(xFile.c_str());
	ofstream yData(yFile.c_str());
	
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			xData << i*DX << " ";
		}
		xData << endl;
	}
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			yData << j*DX << " ";
		}
		yData << endl;
	}
	xData.close();
	yData.close();
}

void writeVelocityToFile(VectorXd &u,VectorXd &v){
	string path = OUTDIR;
	string uFile = path + "data-u.dat";
	string vFile = path + "data-v.dat";
	string mFile = path + "data-m.dat";
	
	ofstream uData(uFile.c_str());
	ofstream vData(vFile.c_str());
	ofstream mData(mFile.c_str());
	for(int j=0;j<NGP;j++){
		uData << (u(j*NGP)+u((NGP-1)+j*NGP))/2 << " ";
		for(int i=1;i<NGP;i++){
			uData << (u(i+j*NGP)+u((i-1)+j*NGP))/2 << " ";
		}
		uData << endl;
	}
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(j!=0)
				vData << (v(i+j*NGP)+v(i+(j-1)*NGP))/2 << " ";
			else
				vData << (v(i)+v(i+(NGP-1)*NGP))/2 << " ";
		}
		vData << endl;
	}
	double u_inp, v_inp;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i!=0)
				u_inp = (u(i+j*NGP)+u((i-1)+j*NGP))/2;
			else
				u_inp = (u(j*NGP)+u((NGP-1)+j*NGP))/2;
			if(j!=0)
				v_inp = (v(i+j*NGP)+v(i+(j-1)*NGP))/2;
			else
				v_inp = (v(i)+v(i+(NGP-1)*NGP))/2;
			mData << sqrt(pow(u_inp,2)+pow(v_inp,2)) << " ";
		}
		mData << endl;
	}
	uData.close();
	vData.close();
	mData.close();
}

void writePressureToFile(VectorXd &phi){
	string path = OUTDIR;
	string pFile = path + "data-p.dat";
	ofstream pData(pFile.c_str());
	
	double phi_nab;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i==0)
				phi_nab = phi((i+1)+j*NGP)+phi((NGP-1)+j*NGP)-4*phi(i+j*NGP);
			else if(i==(NGP-1))
				phi_nab = phi(0+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);
			else
				phi_nab = phi((i+1)+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);
			if(j==0)
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(NGP-1)*NGP);
			else if(j==(NGP-1))
				phi_nab += phi(i+0*NGP)+phi(i+(j-1)*NGP);
			else
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(j-1)*NGP);
			pData << phi(i+j*NGP)-(DT/(RE*2))*phi_nab << " ";
		}
		pData << endl;
	}
	pData.close();
}

void writePhiToFile(VectorXd &phi){
	string path = OUTDIR;
	string phiFile = path + "data-phi.dat";
	ofstream phiData(phiFile.c_str());
	
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			phiData << phi(i+j*NGP) << " ";
		}
		phiData << endl;
	}
	phiData.close();
}

void writeGridToBinary(){
	string path = OUTDIR;
	string xFile = path + "grid-x.bin";
	string yFile = path + "grid-y.bin";
	ofstream xData(xFile.c_str(), ios::out | ios::binary);
	ofstream yData(yFile.c_str(), ios::out | ios::binary);

	double buffer;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			buffer=i*DX;
			xData.write((char*)&buffer,sizeof(buffer));
		}
	}
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			buffer=j*DX;
			yData.write((char*)&buffer,sizeof(buffer));
		}
	}
	xData.close();
	yData.close();
}

void writeVelocityToBinary(VectorXd &u,VectorXd &v,int ts){
	string path = OUTDIR;
	string fileNum = getFilenameDigits(ts);
	string uFile = path + "data-u-" + fileNum + ".bin";
	string vFile = path + "data-v-" + fileNum + ".bin";
	string mFile = path + "data-vmag-" + fileNum + ".bin";
	ofstream uData(uFile.c_str(), ios::out | ios::binary);
	ofstream vData(vFile.c_str(), ios::out | ios::binary);
	ofstream mData(mFile.c_str(), ios::out | ios::binary);
	
	double buffer;
	for(int j=0;j<NGP;j++){
		buffer = (u(j*NGP)+u((NGP-1)+j*NGP))/2;
		uData.write((char*)&buffer,sizeof(buffer));
		for(int i=1;i<NGP;i++){
			buffer = (u(i+j*NGP)+u((i-1)+j*NGP))/2;
			uData.write((char*)&buffer,sizeof(buffer));
		}
	}
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(j!=0)
				buffer = (v(i+j*NGP)+v(i+(j-1)*NGP))/2;
			else
				buffer = (v(i)+v(i+(NGP-1)*NGP))/2;
			vData.write((char*)&buffer,sizeof(buffer));
		}
	}
	double u_inp, v_inp;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i!=0)
				u_inp = (u(i+j*NGP)+u((i-1)+j*NGP))/2;
			else
				u_inp = (u(j*NGP)+u((NGP-1)+j*NGP))/2;
			if(j!=0)
				v_inp = (v(i+j*NGP)+v(i+(j-1)*NGP))/2;
			else
				v_inp = (v(i)+v(i+(NGP-1)*NGP))/2;
			buffer = sqrt(pow(u_inp,2)+pow(v_inp,2));
			mData.write((char*)&buffer,sizeof(buffer));
		}
	}
	
	uData.close();
	vData.close();
	mData.close();
}

void writePressureToBinary(VectorXd &phi, int ts){
	string path = OUTDIR;
	string fileNum = getFilenameDigits(ts);
	string pFile = path + "data-p-" + fileNum + ".bin";
	ofstream pData(pFile.c_str(), ios::out | ios::binary);
	
	double phi_nab, buffer;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i==0)
				phi_nab = phi((i+1)+j*NGP)+phi((NGP-1)+j*NGP)-4*phi(i+j*NGP);
			else if(i==(NGP-1))
				phi_nab = phi(0+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);
			else
				phi_nab = phi((i+1)+j*NGP)+phi((i-1)+j*NGP)-4*phi(i+j*NGP);

			if(j==0)
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(NGP-1)*NGP);
			else if(j==(NGP-1))
				phi_nab += phi(i+0*NGP)+phi(i+(j-1)*NGP);
			else
				phi_nab += phi(i+(j+1)*NGP)+phi(i+(j-1)*NGP);
			buffer = phi(i+j*NGP)-(DT/(RE*2))*phi_nab;
			pData.write((char*)&buffer,sizeof(buffer));
		}
	}
	pData.close();
}

void writeInfoToBinary(int ts){
	string path = OUTDIR;
	if(ts==0){
		struct Buffer{
			char type;
			int ngp,tsmax;
			double lref,dt,dx,re,uref,a,lratio,rev,nu;	
		} buffer;
		buffer.type = TYPE;
		buffer.ngp = NGP;
		buffer.tsmax = TSMAX;
		buffer.lref = LREF;
		buffer.dt = DT;
		buffer.dx = DX;
		if(TYPE=='t'){
			buffer.re = RE;
			buffer.uref = UREF;
			buffer.a = 0;
			buffer.lratio = 0;
			buffer.rev = 0;
			buffer.nu = 0;
		}
		else if(TYPE=='v'){
			buffer.re = 0;
			buffer.uref = 0;
			buffer.a = A;
			buffer.lratio = LRATIO;
			buffer.rev = REV;
			buffer.nu = NU;
		}
		string fileNum = getFilenameDigits(ts);
		string infoFile = path + "settings.bin";
		ofstream infoData(infoFile.c_str(), ios::out | ios::binary);
		infoData.write((char*)&buffer,sizeof(buffer));
		infoData.close();
	}
	else{
		string fileNum = getFilenameDigits(ts);
		string infoFile = path + "info-" + fileNum + ".bin";
		double buffer = ts*DT;
		ofstream infoData(infoFile.c_str(), ios::out | ios::binary);
		infoData.write((char*)&buffer,sizeof(buffer));
		infoData.close();
	}
}

void saveData(VectorXd &u, VectorXd &v, VectorXd &phi, int ts){
	writeInfoToBinary(ts);
	writeVelocityToBinary(u,v,ts);
	writePressureToBinary(phi,ts);
}

string getFilenameDigits(int ts){
	//determine number of digits needed for filename
	int digits = 1;
	int max = TSMAX;
	while (max/=10)
	   digits++;
	stringstream fileNum;
	fileNum.fill('0');
	fileNum.width(digits);
	fileNum << ts;
	return fileNum.str();
}
