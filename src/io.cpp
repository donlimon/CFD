//Functions for input/output
//file: io.cpp
//author: Michael Stumpf

#include <iostream>
#include <time.h>
#include <settings.hpp>

using namespace std;

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
