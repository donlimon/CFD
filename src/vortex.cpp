//Implementation of two counterrotating vortices
//file: vortex.cpp
//author: Michael Stumpf

#include <iostream>
#include <cmath>
#include <settings.hpp>
#include <vortex.hpp>

using namespace std;

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
	vortexCenterX=((NGP-1)/2)*DX + DX/2 - sqrt(2)*b/2;	//lower left
	vortexCenterY=((NGP-1)/2)*DX - sqrt(2)*b/2;
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
				v(i+j*NGP)=v_circ*(delX/r);
			}
			else{
				v_circ=0;
				v(i+j*NGP)=0;
			}
		}
	}
	
	//Vortex 2
	vortexCenterX+=sqrt(2)*b;	//upper right
	vortexCenterY+=sqrt(2)*b;
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
				v(i+j*NGP)+=v_circ*(delX/r);
			}
			else{
				v_circ=0;
				v(i+j*NGP)+=0;
			}
		}
	}
}
