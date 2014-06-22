//Functions for generating/updating coefficient matrices
//file: matrixgen.cpp
//author: Michael Stumpf

#include <iostream>
#include <settings.hpp>
#include <matrixgen.hpp>
#include <testing.hpp>

using namespace std;
using namespace Eigen;

/* LOCAL PROTOTYPES */
inline void assignSurroundingPoints(int i, int j,int &im1,int &ip1, int &jm1, int &jp1);

//Create coefficient matrix A (for momentum/poisson eq.)
//a -> parameter in main diagonal
//b -> parameter in side diagonals
void createCoeffMat(SpMat &mat,double a, double b){
	mat.reserve(VectorXi::Constant(NGP*NGP,5));	//5 entries per row
	for(int i=0;i<(NGP*NGP);i++){
		//main diagonal with irregularities
		mat.insert(i,i)=a;
		if((i%NGP)==0){
			//lines with i=1
			mat.insert(i,(i+1))=b;
			mat.insert(i,(i+NGP-1))=b;
		}
		else if(((i+1)%NGP)==0){
			//lines with i=nx
			mat.insert(i,(i-NGP+1))=b;
			mat.insert(i,(i-1))=b;
		}
		else{
			mat.insert(i,(i+1))=b;
			mat.insert(i,(i-1))=b;
		}
		
		//side diagonals
		if(i<(NGP*(NGP-1))){
			//upper diagonal
			mat.insert(i,(i+NGP))=b;
			if(i<NGP){
				//upper right part
				mat.insert(i,(NGP*(NGP-1)+i))=b;
			}
		}
		if(i>=NGP){
			//lower diagonal
			mat.insert(i,(i-NGP))=b;
			if(i>=(NGP*(NGP-1))){
				//lower left part
				mat.insert(i,(i-NGP*(NGP-1)))=b;
			}
		}
	}
	mat.makeCompressed();
}

void updateRHSpoi(VectorXd &rhs, VectorXd &u, VectorXd &v){
	// poisson eq.
	double  a = -DX/(4*DT);
	int im1,ip1,jm1,jp1;	//i-1,i+1,j-1,j+1
	for(int j=0;j<NGP;j++){
		for(int i=0;i<NGP;i++){
			assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
			rhs(i+j*NGP) = a*(u(i+j*NGP)-u(im1+j*NGP)+v(i+j*NGP)-v(i+jm1*NGP));
		}
	}
}

void updateRHSmom(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp){
	double  a = -1/(4*DX),
			b = 1/(2*RE*DX*DX);
	double Hn, Hp, Gn, In;
	int im1,ip1,jm1,jp1;	//i-1,i+1,j-1,j+1
	if(comp=='u'){
		// momentum eq. for u
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
					
				Hn = ucur(i+j*NGP)*(ucur(ip1+j*NGP)-ucur(im1+j*NGP))
					+interpolate(vcur,i,j,'u')*(ucur(i+jp1*NGP)-ucur(i+jm1*NGP));
				Hp = uprev(i+j*NGP)*(uprev(ip1+j*NGP)-uprev(im1+j*NGP))
					+interpolate(vprev,i,j,'u')*(uprev(i+jp1*NGP)-uprev(i+jm1*NGP));
				Gn = ucur(ip1+j*NGP)+ucur(im1+j*NGP)-4*ucur(i+j*NGP)
					+ucur(i+jp1*NGP)+ucur(i+jm1*NGP);
				
				//calculate I and update RHS
				In = 1/DT*ucur(i+j*NGP);
				rhs(i+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
			}
		} 
	}
	else if(comp=='v'){
		// momentum eq. for v
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
				
				Hn = interpolate(ucur,i,j,'v')*(vcur(ip1+j*NGP)-vcur(im1+j*NGP))
					+vcur(i+j*NGP)*(vcur(i+jp1*NGP)-vcur(i+jm1*NGP));
				Hp = interpolate(uprev,i,j,'v')*(vprev(ip1+j*NGP)-vprev(im1+j*NGP))
					+vprev(i+j*NGP)*(vprev(i+jp1*NGP)-vprev(i+jm1*NGP));
				Gn = vcur(ip1+j*NGP)+vcur(im1+j*NGP)-4*vcur(i+j*NGP)
					+vcur(i+jp1*NGP)+vcur(i+jm1*NGP);

				//calculate I and update RHS
				In = 1/DT*vcur(i+j*NGP);
				rhs(i+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
			}
		} 
	}
}


inline double interpolate(VectorXd &field, int i, int j, char comp){
	int im1,ip1,jm1,jp1;
	assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
	if(comp=='u'){
		//value of v on grid of u
		//v @ u(i,j) = 0.25*[v(i,j)+v(i+1,j)+v(i,j-1)+v(i+1,j-1)]
		return 0.25*(field(i+j*NGP)+field(ip1+j*NGP)
				+field(i+jm1*NGP)+field(ip1+jm1*NGP));
	}
	else if(comp=='v'){
		//value of u on grid of v
		//u @ v(i,j) = 0.25*[u(i-1,j+1)+u(i,j+1)+u(i-1,j)+u(i,j)]
		return 0.25*(field(i+j*NGP)+field(im1+j*NGP)
				+field(i+jp1*NGP)+field(im1+jp1*NGP));
	}
	else{
		cerr << "ERROR: interpolation error" << endl;
		return 0;
	}
}

void solveCorrector(VectorXd &velcur, VectorXd &velpre, VectorXd &phi, char comp){
	//v(n+1) = v* - a*[phi(i+1)-phi(i)]
	double a = DT/DX;
	int im1,ip1,jm1,jp1;
	if(comp=='u'){
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
				velcur(i+j*NGP) = velpre(i+j*NGP)
									-a*(phi(ip1+j*NGP)-phi(i+j*NGP));
			}
		}
	}
	else if(comp=='v'){
		for(int i=0;i<NGP;i++){
			for(int j=0;j<NGP;j++){
				assignSurroundingPoints(i,j,im1,ip1,jm1,jp1);
				velcur(i+j*NGP) = velpre(i+j*NGP)
									-a*(phi(i+jp1*NGP)-phi(i+j*NGP));
			}
		}
	}
}

/* LOCAL FUNCTIONS */
inline void assignSurroundingPoints(int i, int j,int &im1,int &ip1, int &jm1, int &jp1){
	if(i==0)
		im1=(NGP-1);
	else 
		im1=i-1;
		
	if(i==(NGP-1))
		ip1=0;
	else
		ip1=i+1;
	
	if(j==0)
		jm1=(NGP-1);
	else
		jm1=j-1;
	
	if(j==(NGP-1))
		jp1=0;
	else
		jp1=j+1;
}
