//Implementation of two counterrotating vortices
//file: vortex.cpp
//author: Michael Stumpf

void initializeVortex(VectorXd &u, VectorXd &v){
	//Two vortices are being initialized
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
	vortexCenterX=((NGP-1)/2)*DX + DX/2;
	vortexCenterY=((NGP-1)/2)*DX;
	for(int i=0;i<NGP;i++){
		xCoordU=i*DX + DX/2;
		xCoordV=i*DX;
		for(int j=0;j<NGP;j++){
			yCoordU=j*DX;
			yCoordV=j*DX + DX/2;
			//Calculate u velocity
			delX=xCoordU-vortexCenterX;
			delY=yCoordU-vortexCenterY;
			r=sqrt(pow(delX,2)+pow(delY,2));			
			v_circ=gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
			u(i+j*NGP)=v_circ*(delY/r);
			//Calculate v velocity
			delX=xCoordV-vortexCenterX;
			delY=yCoordV-vortexCenterY;
			r=sqrt(pow(delX,2)+pow(delY,2));			
			v_circ=gamma/(2*M_PI*r)*(1-exp(-(r*r)/(A*A)));
			v(i+j*NGP)=v_circ*(delX/r);
		}
	}
}
