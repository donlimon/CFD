//Implementation of two counterrotating vortices
//file: vortex.cpp
//author: Michael Stumpf

#include <iostream>
#include <cmath>
#include <vector>

#include <settings.hpp>
#include <vortex.hpp>

using namespace std;

/* LOCAL PROTOTYPES */
int determineDirection(double a,double b,double c,double d);
void sortCandidates(vector<int>& cand_i, vector<int>& cand_j, vector<double>& cand_val);
void checkVortexOrder(location prevPos[], location newPos[]);

/* FUNCTIONS PROVIDED FOR OTHER CPP */
void initializeVortex(VectorXd &u, VectorXd &v){
	//two vortices with distance A are being initialized
	//each velocity field is calculated independently
	//vortex center is located on u grid!
	double v_circ, 				//circumferential velocity
		   xCoordU, yCoordU, 	//coordinates of current grid points (u component)
		   xCoordV, yCoordV,	//coordinates of current grid points (v component)
		   vortexCenterX, vortexCenterY,	//coordinates of vortex center
		   delX, delY,			//Delta X and Delta Y
		   r;					//distance from vortex center
	double gamma=REV*NU, b=LRATIO*A;
	//Vortex 1
	vortexCenterX=((NGP-1)/2)*DX + DX/2 - b/(2*sqrt(2));	//lower left
	vortexCenterY=((NGP-1)/2)*DX - b/(2*sqrt(2));
	for(int i=0;i<NGP;i++){
		xCoordU=i*DX + DX/2;
		xCoordV=i*DX;
		for(int j=0;j<NGP;j++){
			yCoordU=j*DX;
			yCoordV=j*DX + DX/2;
			//Calculate u velocity
			delX=xCoordU-vortexCenterX;
			delY=yCoordU-vortexCenterY;
			if(delX > LREF/2)
				delX=-(LREF-delX);
			else if(delX < (-LREF/2))
				delX=(LREF-abs(delX));
				
			if(delY > LREF/2)
				delY=-(LREF-delY);
			else if(delY < (-LREF/2))
				delY=(LREF-abs(delY));
				
			r=sqrt(pow(delX,2)+pow(delY,2));
			if(r!=0){		
				v_circ=gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
				u(i+j*NGP)=v_circ*(delY/r);
			}
			else{
				v_circ=0;
				u(i+j*NGP)=0;
			}
			//Calculate v velocity
			delX=xCoordV-vortexCenterX;
			delY=yCoordV-vortexCenterY;
			if(delX > LREF/2)
				delX=-(LREF-delX);
			else if(delX < (-LREF/2))
				delX=(LREF-abs(delX));
				
			if(delY > LREF/2)
				delY=-(LREF-delY);
			else if(delY < (-LREF/2))
				delY=(LREF-abs(delY));
				
			r=sqrt(pow(delX,2)+pow(delY,2));	
			if(r!=0){
				v_circ=gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
				v(i+j*NGP)=v_circ*(-delX/r);
			}
			else{
				v_circ=0;
				v(i+j*NGP)=0;
			}
		}
	}
	
	//Vortex 2
	vortexCenterX+=b/sqrt(2);	//upper right
	vortexCenterY+=b/sqrt(2);
	for(int i=0;i<NGP;i++){
		xCoordU=i*DX + DX/2;
		xCoordV=i*DX;
		for(int j=0;j<NGP;j++){
			yCoordU=j*DX;
			yCoordV=j*DX + DX/2;
			//Calculate u velocity
			delX=xCoordU-vortexCenterX;
			delY=yCoordU-vortexCenterY;
			if(delX > LREF/2)
				delX=-(LREF-delX);
			else if(delX < (-LREF/2))
				delX=(LREF-abs(delX));
				
			if(delY > LREF/2)
				delY=-(LREF-delY);
			else if(delY < (-LREF/2))
				delY=(LREF-abs(delY));
			
			r=sqrt(pow(delX,2)+pow(delY,2));
			if(r!=0){		
				v_circ=-gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
				u(i+j*NGP)+=v_circ*(delY/r);
			}
			else{
				v_circ=0;
				u(i+j*NGP)+=0;
			}
			//Calculate v velocity
			delX=xCoordV-vortexCenterX;
			delY=yCoordV-vortexCenterY;
			if(delX > LREF/2)
				delX=-(LREF-delX);
			else if(delX < (-LREF/2))
				delX=(LREF-abs(delX));
				
			if(delY > LREF/2)
				delY=-(LREF-delY);
			else if(delY < (-LREF/2))
				delY=(LREF-abs(delY));
			
			r=sqrt(pow(delX,2)+pow(delY,2));	
			if(r!=0){
				v_circ=-gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
				v(i+j*NGP)+=v_circ*(-delX/r);
			}
			else{
				v_circ=0;
				v(i+j*NGP)+=0;
			}
		}
	}
}

void calcVorticity(VectorXd &omega, VectorXd &u, VectorXd &v){
	double u_up, u_low, v_up, v_low,
			i_up, i_low, j_up, j_low,
			dvdx, dudy;
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i!=0 && i!=NGP-1)
				i_up = i+1, i_low = i-1;
			else if(i==0)
				i_up = i+1, i_low = NGP-1;
			else
				i_up = 0, i_low = i-1;
				
			if(j!=0 && j!=NGP-1)
				j_up = j+1, j_low = j-1;
			else if(j==0)
				j_up = j+1, j_low = NGP-1;
			else
				j_up = 0, j_low = j-1;
			
			u_up = (u(i+j_up*NGP)+u(i_low+j_up*NGP))/2;
			u_low = (u(i+j_low*NGP)+u(i_low+j_low*NGP))/2;
			v_up = (v(i_up+j*NGP)+v(i_up+j_low*NGP))/2;
			v_low = (v(i_low+j*NGP)+v(i_low+j_low*NGP))/2;
			
			dvdx = (v_up - v_low)/(2*DX);
			dudy = (u_up - u_low)/(2*DX);
			
			omega(i+j*NGP) = dvdx - dudy;
		}
	}
}

void locateVortex(location *newPos, VectorXd &omega, int ts){
	int i_low, i_up, j_low, j_up;
	static bool initialized = false;
	bool alternative = false;
	static location prevPos[NUMVOR]; // coordinates of maxima at last ts
	//create vectors for possible candidates
	const int buffer = 1000; //number of candidates should not exceed buffer
	vector<int> candidates_i;
	vector<int> candidates_j;
	vector<double> candidates_val;
	candidates_i.reserve(buffer);
	candidates_j.reserve(buffer);
	candidates_val.reserve(buffer);
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			if(i!=0 && i!=(NGP-1))
				i_low = i-1, i_up = i+1;
			else if(i==0)
				i_low = NGP-1, i_up = i+1;
			else
				i_low = i-1, i_up = 0;
		
			if(j!=0 && j!=(NGP-1))
				j_low = j-1, j_up = j+1;
			else if(j==0)
				j_low = NGP-1, j_up = j+1;
			else
				j_low = j-1, j_up = 0;
		
			// Check if surrounding points contain higher values of omega
			if(abs(omega(i+j*NGP)) > abs(omega(i_low+j*NGP)) 
				&& abs(omega(i+j*NGP)) > abs(omega(i+j_low*NGP))
				&& abs(omega(i+j*NGP)) > abs(omega(i_up+j*NGP)) 
				&& abs(omega(i+j*NGP)) > abs(omega(i+j_up*NGP)))
			{
				candidates_i.push_back(i);
				candidates_j.push_back(j);
				candidates_val.push_back(abs(omega(i+j*NGP)));
			}
		}
	}
	if(candidates_i.size() > NUMVOR){
		//cout << "Found " << candidates_i.size() << " maxima @ time step #" << ts << endl;
		if(candidates_i.size() != candidates_j.size()){
			cerr << "ERROR: size of canditates vectors is not equal @ time step #" << ts << endl;
			exit(3);
		}
		
		//sort candidates in ascending order
		sortCandidates(candidates_i,candidates_j,candidates_val);
		//assign candidates with largest vorticity to new position
		for(int i=0;i<NUMVOR;i++){
			newPos[i].i = candidates_i[candidates_i.size()-1-i];
			newPos[i].j = candidates_j[candidates_j.size()-1-i];
		}
		//check if positions have not changed (else switch)
		if(initialized)
			checkVortexOrder(prevPos,newPos);
		//assign time step
		for(int i=0;i<NUMVOR;i++){
			newPos[i].ts = ts;
		}
		//copy new position to old position
		for(int i=0;i<NUMVOR;i++){
			prevPos[i].i = newPos[i].i;
			prevPos[i].j = newPos[i].j;
			prevPos[i].ts = newPos[i].ts;
		}
		//set flags
		initialized = true;
	}
	else if(candidates_i.size() < NUMVOR){
		cerr << "Found less than " << NUMVOR << " maxima @ time step #" << ts << endl;
	}
	else{
		initialized = false;
	}
//	cout << "Position of Minimum 1 @ time step #" << newPos[0].ts << " i=" << newPos[0].i << " j="
//			 << newPos[0].j << endl;
//	cout << "Position of Minimum 2 @ time step #" << newPos[0].ts << " i=" << newPos[1].i << " j="
//			 << newPos[1].j << endl;
}

/* LOCAL FUNCTIONS (ONLY ACCESSIBLE IN THIS CPP) */
int determineDirection(double a,double b,double c,double d){
	double tmp[4] = {a,b,c,d}, maxVal;
	int index;
	maxVal=tmp[0];
	index=0;
	for(int i=1;i<4;i++){
		if(tmp[i] > maxVal){
			maxVal = tmp[i];
			index = i;
		}
	}
	return ++index;
}

void sortCandidates(vector<int>& cand_i, vector<int>& cand_j, vector<double>& cand_val){
	//vectors usually contain less than 10 entries, so bubblesort is fine
	int tmp;
	for(int i=0;i<cand_val.size()-1;i++) 
	{
		for(int j=0;j<cand_val.size()-i-1;j++) 
		{
			if(cand_val[j]>cand_val[j+1]) 
		    {
		    	//swap values
		 		tmp = cand_val[j];
		 		cand_val[j] = cand_val[j+1];
		 		cand_val[j+1] = tmp;
		 		//swap i
		 		tmp = cand_i[j];
		 		cand_i[j] = cand_i[j+1];
		 		cand_i[j+1] = tmp;
		 		//swap j
		 		tmp = cand_j[j];
		 		cand_j[j] = cand_j[j+1];
		 		cand_j[j+1] = tmp;
	 	    }
		}
	}
}

void checkVortexOrder(location prevPos[], location newPos[]){
	double distance, distanceMin;
	int minPosition;
	location tmpPos[NUMVOR];
	bool posBool[NUMVOR]={false};
	//look for minimum distance from previous location and assign to temporary variable
	for(int i=0;i<NUMVOR;i++){
		for(int j=0;j<NUMVOR;j++){
			distance = sqrt(pow(newPos[j].i-prevPos[i].i,2)
								+pow(newPos[j].j-prevPos[i].j,2));
			if(j==0){
				distanceMin = distance;
				minPosition = 0;
			}
			else if(distance < distanceMin){
				distanceMin = distance;
				minPosition = j;
			}
		}
		tmpPos[i] = newPos[minPosition];
		if(posBool[minPosition]==true)
			cerr << "WARNING: conflict while assigning vortex order" << endl;
		else
			posBool[minPosition] = true;
	}
	//copy temporary variable to new location
	for(int i=0;i<NUMVOR;i++){
		newPos[i] = tmpPos[i];
	}
}
