//Functions for generating/updating coefficient matrices
//file: matrixgen.cpp
//author: Michael Stumpf

#include <iostream>
#include <settings.hpp>
#include <matrixgen.hpp>

using namespace std;
using namespace Eigen;

//Create coefficient matrix A (for momentum/poisson eq.)
//a -> parameter in main diagonal
//b -> parameter in side diagonals
SpMat createCoeffMat(double a, double b){
	SpMat mat(NGP*NGP,NGP*NGP);
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
	return mat;
}

void updateRHS(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp){
	if(comp=='u'){
		// momentum eq. for u
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		double  a = -1/(4*DX),
				b = 1/(2*RE*DX*DX);
		double Hn, Hp, Gn, In;
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				//first part of H and G (i+-1)
				if(i==0){
					Hn = ucur(i+j*NGP)*(ucur((i+1)+j*NGP)-ucur((NGP-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev((i+1)+j*NGP)-uprev((NGP-1)+j*NGP));
					Gn = ucur((i+1)+j*NGP)+ucur((NGP-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				else if(i==(NGP-1)){
					Hn = ucur(i+j*NGP)*(ucur(0+j*NGP)-ucur((i-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev(0+j*NGP)-uprev((i-1)+j*NGP));
					Gn = ucur(0+j*NGP)+ucur((i-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				else{
					Hn = ucur(i+j*NGP)*(ucur((i+1)+j*NGP)-ucur((i-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev((i+1)+j*NGP)-uprev((i-1)+j*NGP));
					Gn = ucur((i+1)+j*NGP)+ucur((i-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				//second part of H and G (j+-1)
				if(j==0){
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+(j+1)*NGP)-ucur(i+(NGP-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+(j+1)*NGP)-uprev(i+(NGP-1)*NGP));
					Gn += ucur(i+(j+1)*NGP)+ucur(i+(NGP-1)*NGP);
				}
				else if(j==(NGP-1)){
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+0*NGP)-ucur(i+(j-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+0*NGP)-uprev(i+(j-1)*NGP));
					Gn += ucur(i+0*NGP)+ucur(i+(j-1)*NGP);
				}
				else{
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+(j+1)*NGP)-ucur(i+(j-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+(j+1)*NGP)-uprev(i+(j-1)*NGP));
					Gn += ucur(i+(j+1)*NGP)+ucur(i+(j-1)*NGP);
				}
				//calculate I and update RHS
				In = 1/DT*ucur(i+j*NGP);
				rhs(i+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
			}
		} 
	}
	else if(comp=='v'){
		// momentum eq. for v
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		double  a = -1/(4*DX),
				b = 1/(2*RE*DX*DX);
		double Hn, Hp, Gn, In;
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				//first part of H and G (i+-1)
				if(i==0){
					Hn = interpolate(ucur,i,j,'v')*(vcur((i+1)+j*NGP)-vcur((NGP-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev((i+1)+j*NGP)-vprev((NGP-1)+j*NGP));
					Gn = vcur((i+1)+j*NGP)+vcur((NGP-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				else if(i==(NGP-1)){
					Hn = interpolate(ucur,i,j,'v')*(vcur(0+j*NGP)-vcur((i-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev(0+j*NGP)-vprev((i-1)+j*NGP));
					Gn = vcur(0+j*NGP)+vcur((i-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				else{
					Hn = interpolate(ucur,i,j,'v')*(vcur((i+1)+j*NGP)-vcur((i-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev((i+1)+j*NGP)-vprev((i-1)+j*NGP));
					Gn = vcur((i+1)+j*NGP)+vcur((i-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				//second part of H and G (j+-1)
				if(j==0){
					Hn += vcur(i+j*NGP)*(vcur(i+(j+1)*NGP)-vcur(i+(NGP-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+(j+1)*NGP)-vprev(i+(NGP-1)*NGP));
					Gn += vcur(i+(j+1)*NGP)+vcur(i+(NGP-1)*NGP);
				}
				else if(j==(NGP-1)){
					Hn += vcur(i+j*NGP)*(vcur(i+0*NGP)-vcur(i+(j-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+0*NGP)-vprev(i+(j-1)*NGP));
					Gn += vcur(i+0*NGP)+vcur(i+(j-1)*NGP);
				}
				else{
					Hn += vcur(i+j*NGP)*(vcur(i+(j+1)*NGP)-vcur(i+(j-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+(j+1)*NGP)-vprev(i+(j-1)*NGP));
					Gn += vcur(i+(j+1)*NGP)+vcur(i+(j-1)*NGP);
				}
				//calculate I and update RHS
				In = 1/DT*vcur(i+j*NGP);
				rhs(i+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
			}
		} 
	}
	else if(comp=='p'){
		// poisson eq.
		// a*(I+J)
		double  a = -DX/(4*DT);
		double I,J;
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				//I
				if(i==0){
					I=ucur(i+j*NGP)-ucur((NGP-1)+j*NGP);
				}
				else{
					I=ucur(i+j*NGP)-ucur((i-1)+j*NGP);
				}
				//J
				if(j==0){
					J=vcur(i+j*NGP)-vcur(i+(NGP-1)*NGP);
				}
				else{
					J=vcur(i+j*NGP)-vcur(i+(j-1)*NGP);
				}
				//update RHS
				rhs(i+j*NGP) = a*(I+J);
			}
		} 
	}
}

double interpolate(VectorXd field, int i, int j, char comp){
	if(comp=='u'){
		//value of v on grid of u
		//v @ u(i,j) = 0.25*[v(i,j)+v(i+1,j)+v(i,j-1)+v(i+1,j-1)]
		//w1 -> i+1, w2 -> j-1
		int w1,w2;
		if(i==NGP-1){
			w1=0;
		}
		else{
			w1=i+1;
		}
		if(j==0){
			w2=NGP-1;
		}
		else{
			w2=j-1;
		}
		return 0.25*(field(i+j*NGP)+field(w1+j*NGP)
				+field(i+w2*NGP)+field(w1+w2*NGP));
	}
	else if(comp=='v'){
		//value of u on grid of v
		//u @ v(i,j) = 0.25*[u(i-1,j+1)+u(i,j+1)+u(i-1,j)+u(i,j)]
		//w1 -> i-1, w2 -> j+1
		int w1,w2;
		if(i==0){
			w1=NGP-1;
		}
		else{
			w1=i-1;
		}
		if(j==NGP-1){
			w2=0;
		}
		else{
			w2=j+1;
		}
		return 0.25*(field(i+j*NGP)+field(w1+j*NGP)
				+field(i+w2*NGP)+field(w1+w2*NGP));
	}
	else{
		return 0;
	}
}
