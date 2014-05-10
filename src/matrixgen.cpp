//Functions for generating/updating coefficient matrices
//file: matrixgen.cpp
//author: Michael Stumpf

#include <iostream>
#include <settings.hpp>
#include <matrixgen.hpp>
#include <testing.hpp>

using namespace std;
using namespace Eigen;

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

void updateRHS(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp){
	double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
	if(comp=='u'){
		// momentum eq. for u
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		double  a = -1/(4*DX),
				b = 1/(2*RE*DX*DX);
		double Hn, Hp, Gn, In;
		for(int j=0;j<NGP;j++){
			for(int i=0;i<NGP;i++){
				//first part of H and G (i+-1)
				if(i>0 && i<(NGP-1)){
					Hn = ucur(i+j*NGP)*(ucur((i+1)+j*NGP)-ucur((i-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev((i+1)+j*NGP)-uprev((i-1)+j*NGP));
					Gn = ucur((i+1)+j*NGP)+ucur((i-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				else if(i==(NGP-1)){
					Hn = ucur(i+j*NGP)*(ucur(0+j*NGP)-ucur((i-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev(0+j*NGP)-uprev((i-1)+j*NGP));
					Gn = ucur(0+j*NGP)+ucur((i-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				else{
					Hn = ucur(i+j*NGP)*(ucur((i+1)+j*NGP)-ucur((NGP-1)+j*NGP));
					Hp = uprev(i+j*NGP)*(uprev((i+1)+j*NGP)-uprev((NGP-1)+j*NGP));
					Gn = ucur((i+1)+j*NGP)+ucur((NGP-1)+j*NGP)-4*ucur(i+j*NGP);
				}
				//second part of H and G (j+-1)
				if(j>0 && j<(NGP-1)){
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+(j+1)*NGP)-ucur(i+(j-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+(j+1)*NGP)-uprev(i+(j-1)*NGP));
					Gn += ucur(i+(j+1)*NGP)+ucur(i+(j-1)*NGP);
				}
				else if(j==(NGP-1)){
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+0*NGP)-ucur(i+(j-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+0*NGP)-uprev(i+(j-1)*NGP));
					Gn += ucur(i+0*NGP)+ucur(i+(j-1)*NGP);
				}
				else{
					Hn += interpolate(vcur,i,j,'u')*(ucur(i+(j+1)*NGP)-ucur(i+(NGP-1)*NGP));
					Hp += interpolate(vprev,i,j,'u')*(uprev(i+(j+1)*NGP)-uprev(i+(NGP-1)*NGP));
					Gn += ucur(i+(j+1)*NGP)+ucur(i+(NGP-1)*NGP);
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
				if(i>0 && i<NGP-1){
					//inner region
					Hn = interpolate(ucur,i,j,'v')*(vcur((i+1)+j*NGP)-vcur((i-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev((i+1)+j*NGP)-vprev((i-1)+j*NGP));
					Gn = vcur((i+1)+j*NGP)+vcur((i-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				else if(i==(NGP-1)){
					//right boundary
					Hn = interpolate(ucur,i,j,'v')*(vcur(0+j*NGP)-vcur((i-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev(0+j*NGP)-vprev((i-1)+j*NGP));
					Gn = vcur(0+j*NGP)+vcur((i-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				else{
					//left boundary
					Hn = interpolate(ucur,i,j,'v')*(vcur((i+1)+j*NGP)-vcur((NGP-1)+j*NGP));
					Hp = interpolate(uprev,i,j,'v')*(vprev((i+1)+j*NGP)-vprev((NGP-1)+j*NGP));
					Gn = vcur((i+1)+j*NGP)+vcur((NGP-1)+j*NGP)-4*vcur(i+j*NGP);
				}
				//second part of H and G (j+-1)
				if(j>0 && j<(NGP-1)){
					//inner region
					Hn += vcur(i+j*NGP)*(vcur(i+(j+1)*NGP)-vcur(i+(j-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+(j+1)*NGP)-vprev(i+(j-1)*NGP));
					Gn += vcur(i+(j+1)*NGP)+vcur(i+(j-1)*NGP);
				}
				else if(j==(NGP-1)){
					//upper boundary
					Hn += vcur(i+j*NGP)*(vcur(i+0*NGP)-vcur(i+(j-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+0*NGP)-vprev(i+(j-1)*NGP));
					Gn += vcur(i+0*NGP)+vcur(i+(j-1)*NGP);
				}
				else{
					//lower boundary
					Hn += vcur(i+j*NGP)*(vcur(i+(j+1)*NGP)-vcur(i+(NGP-1)*NGP));
					Hp += vprev(i+j*NGP)*(vprev(i+(j+1)*NGP)-vprev(i+(NGP-1)*NGP));
					Gn += vcur(i+(j+1)*NGP)+vcur(i+(NGP-1)*NGP);
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
	double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
	cout << "Update normal" << endl;
    cout << "Wall Time = " << wall1 - wall0 << endl;
    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
}

void updateRHSopt(VectorXd &rhs, VectorXd &ucur, VectorXd &uprev, 
						VectorXd &vcur, VectorXd &vprev, char comp){
	double wall0 = get_wall_time();
    double cpu0  = get_cpu_time();
	if(comp=='u'){
		// momentum eq. for u
		// a*[3*H(n)-H(n-1)]+b*G(n)+I(n)
		double  a = -1/(4*DX),
				b = 1/(2*RE*DX*DX);
		double Hn, Hp, Gn, In;
		
		double i_up=NGP-1;
		double j_up=NGP-1;
		//lower boundary (j=0)  [no corner points]
		for(int i=1;i<(NGP-1);i++){
			Hn = ucur(i)*(ucur(i+1)-ucur(i-1))
				+interpolate(vcur,i,0,'u')*(ucur(i+NGP)-ucur(i+(NGP-1)*NGP));
			Hp = uprev(i)*(uprev(i+1)-uprev(i-1))
				+interpolate(vprev,i,0,'u')*(uprev(i+NGP)-uprev(i+(NGP-1)*NGP));
			Gn = ucur(i+1)+ucur(i-1)-4*ucur(i)
				+ucur(i+NGP)+ucur(i+(NGP-1)*NGP);
			In = 1/DT*ucur(i);
			rhs(i) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//upper boundary (j=NGP-1) [no corner points]
		for(int i=1;i<(NGP-1);i++){
			Hn = ucur(i+j_up*NGP)*(ucur(i+1+j_up*NGP)-ucur(i-1+j_up*NGP))
				+interpolate(vcur,i,j_up,'u')*(ucur(i)-ucur(i+(j_up-1)*NGP));
			Hp = uprev(i+j_up*NGP)*(uprev(i+1+j_up*NGP)-uprev(i-1+j_up*NGP))
				+interpolate(vprev,i,j_up,'u')*(uprev(i)-uprev(i+(j_up-1)*NGP));
			Gn = ucur(i+1+j_up*NGP)+ucur(i-1+j_up*NGP)-4*ucur(i+j_up*NGP)
				+ucur(i)+ucur(i+(j_up-1)*NGP);
			In = 1/DT*ucur(i+j_up*NGP);
			rhs(i+j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//left boundary (i=0) [no corner points]
		for(int j=1;j<(NGP-1);j++){
			Hn = ucur(j*NGP)*(ucur(1+j*NGP)-ucur((NGP-1)+j*NGP))
				+interpolate(vcur,0,j,'u')*(ucur((j+1)*NGP)-ucur((j-1)*NGP));
			Hp = uprev(j*NGP)*(uprev(1+j*NGP)-uprev((NGP-1)+j*NGP))
				+interpolate(vprev,0,j,'u')*(uprev((j+1)*NGP)-uprev((j-1)*NGP));
			Gn = ucur(1+j*NGP)+ucur((NGP-1)+j*NGP)-4*ucur(j*NGP)
				+ucur((j+1)*NGP)+ucur((j-1)*NGP);
			In = 1/DT*ucur(j*NGP);
			rhs(j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//right boundary (i=NGP-1) [no corner points]
		for(int j=1;j<(NGP-1);j++){
			Hn = ucur(i_up+j*NGP)*(ucur(j*NGP)-ucur((i_up-1)+j*NGP))
				+interpolate(vcur,i_up,j,'u')*(ucur(i_up+(j+1)*NGP)-ucur(i_up+(j-1)*NGP));
			Hp = uprev(i_up+j*NGP)*(uprev(j*NGP)-uprev((i_up-1)+j*NGP))
				+interpolate(vprev,i_up,j,'u')*(uprev(i_up+(j+1)*NGP)-uprev(i_up+(j-1)*NGP));
			Gn = ucur(j*NGP)+ucur((i_up-1)+j*NGP)-4*ucur(i_up+j*NGP)
				+ucur(i_up+(j+1)*NGP)+ucur(i_up+(j-1)*NGP);
			In = 1/DT*ucur(i_up+j*NGP);
			rhs(i_up+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//corner points
		//i=0,j=0
		Hn = ucur(0)*(ucur(1)-ucur(NGP-1))
			+interpolate(vcur,0,0,'u')*(ucur(NGP)-ucur((NGP-1)*NGP));
		Hp = uprev(0)*(uprev(1)-uprev(NGP-1))
			+interpolate(vprev,0,0,'u')*(uprev(NGP)-uprev((NGP-1)*NGP));
		Gn = ucur(1)+ucur(NGP-1)-4*ucur(0)
			+ucur(NGP)+ucur((NGP-1)*NGP);
		In = 1/DT*ucur(0);
		rhs(0) = a*(3*Hn-Hp)+b*Gn+In;
		//i=NGP-1,j=0
		Hn = ucur(i_up)*(ucur(0)-ucur(i_up-1))
			+interpolate(vcur,i_up,0,'u')*(ucur(i_up+NGP)-ucur(i_up+(NGP-1)*NGP));
		Hp = uprev(i_up)*(uprev(0)-uprev(i_up-1))
			+interpolate(vprev,i_up,0,'u')*(uprev(i_up+NGP)-uprev(i_up+(NGP-1)*NGP));
		Gn = ucur(0)+ucur(i_up-1)-4*ucur(i_up)
			+ucur(i_up+NGP)+ucur(i_up+(NGP-1)*NGP);
		In = 1/DT*ucur(i_up);
		rhs(i_up) = a*(3*Hn-Hp)+b*Gn+In;
		//i=0,j=NGP-1
		Hn = ucur(j_up*NGP)*(ucur(1+j_up*NGP)-ucur((NGP-1)+j_up*NGP))
			+interpolate(vcur,0,j_up,'u')*(ucur(0)-ucur((j_up-1)*NGP));
		Hp = uprev(j_up*NGP)*(uprev(1+j_up*NGP)-uprev((NGP-1)+j_up*NGP))
			+interpolate(vprev,0,j_up,'u')*(uprev(0)-uprev((j_up-1)*NGP));
		Gn = ucur(1+j_up*NGP)+ucur((NGP-1)+j_up*NGP)-4*ucur(j_up*NGP)
			+ucur(0)+ucur((j_up-1)*NGP);
		In = 1/DT*ucur(j_up*NGP);
		rhs(j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		//i=NGP-1,j=NGP-1
		Hn = ucur(i_up+j_up*NGP)*(ucur(j_up*NGP)-ucur((i_up-1)+j_up*NGP))
			+interpolate(vcur,i_up,j_up,'u')*(ucur(i_up)-ucur(i_up+(j_up-1)*NGP));
		Hp = uprev(i_up+j_up*NGP)*(uprev(j_up*NGP)-uprev((i_up-1)+j_up*NGP))
			+interpolate(vprev,i_up,j_up,'u')*(uprev(i_up)-uprev(i_up+(j_up-1)*NGP));
		Gn = ucur(j_up*NGP)+ucur((i_up-1)+j_up*NGP)-4*ucur(i_up+j_up*NGP)
			+ucur(i_up)+ucur(i_up+(j_up-1)*NGP);
		In = 1/DT*ucur(i_up+j_up*NGP);
		rhs(i_up+j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		//inner region
		for(int j=1;j<NGP-1;j++){
			for(int i=1;i<NGP-1;i++){
				Hn = ucur(i+j*NGP)*(ucur((i+1)+j*NGP)-ucur((i-1)+j*NGP))
					+interpolate(vcur,i,j,'u')*(ucur(i+(j+1)*NGP)-ucur(i+(j-1)*NGP));
				Hp = uprev(i+j*NGP)*(uprev((i+1)+j*NGP)-uprev((i-1)+j*NGP))
					+interpolate(vprev,i,j,'u')*(uprev(i+(j+1)*NGP)-uprev(i+(j-1)*NGP));
				Gn = ucur((i+1)+j*NGP)+ucur((i-1)+j*NGP)-4*ucur(i+j*NGP)
					+ucur(i+(j+1)*NGP)+ucur(i+(j-1)*NGP);
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
		
		double i_up=NGP-1;
		double j_up=NGP-1;
		//lower boundary (j=0)  [no corner points]
		for(int i=1;i<(NGP-1);i++){
			Hn = interpolate(ucur,i,0,'v')*(vcur(i+1)-vcur(i-1))
				+vcur(i)*(vcur(i+NGP)-vcur(i+(NGP-1)*NGP));
			Hp = interpolate(uprev,i,0,'v')*(vprev(i+1)-vprev(i-1))
				+vprev(i)*(vprev(i+NGP)-vprev(i+(NGP-1)*NGP));
			Gn = vcur(i+1)+vcur(i-1)-4*vcur(i)
				+vcur(i+NGP)+vcur(i+(NGP-1)*NGP);
			In = 1/DT*vcur(i);
			rhs(i) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//upper boundary (j=NGP-1) [no corner points]
		for(int i=1;i<(NGP-1);i++){
			Hn = interpolate(ucur,i,j_up,'v')*(vcur(i+1+j_up*NGP)-vcur(i-1+j_up*NGP))
				+vcur(i+j_up*NGP)*(vcur(i)-vcur(i+(j_up-1)*NGP));
			Hp = interpolate(uprev,i,j_up,'v')*(vprev(i+1+j_up*NGP)-vprev(i-1+j_up*NGP))
				+vprev(i+j_up*NGP)*(vprev(i)-vprev(i+(j_up-1)*NGP));
			Gn = vcur(i+1+j_up*NGP)+vcur(i-1+j_up*NGP)-4*vcur(i+j_up*NGP)
				+vcur(i)+vcur(i+(j_up-1)*NGP);
			In = 1/DT*vcur(i+j_up*NGP);
			rhs(i+j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//left boundary (i=0) [no corner points]
		for(int j=1;j<(NGP-1);j++){
			Hn = interpolate(ucur,0,j,'v')*(vcur(1+j*NGP)-vcur((NGP-1)+j*NGP))
				+vcur(j*NGP)*(vcur((j+1)*NGP)-vcur((j-1)*NGP));
			Hp = interpolate(uprev,0,j,'v')*(vprev(1+j*NGP)-vprev((NGP-1)+j*NGP))
				+vprev(j*NGP)*(vprev((j+1)*NGP)-vprev((j-1)*NGP));
			Gn = vcur(1+j*NGP)+vcur((NGP-1)+j*NGP)-4*vcur(j*NGP)
				+vcur((j+1)*NGP)+vcur((j-1)*NGP);
			In = 1/DT*vcur(j*NGP);
			rhs(j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//right boundary (i=NGP-1) [no corner points]
		for(int j=1;j<(NGP-1);j++){
			Hn = interpolate(ucur,i_up,j,'v')*(vcur(j*NGP)-vcur((i_up-1)+j*NGP))
				+vcur(i_up+j*NGP)*(vcur(i_up+(j+1)*NGP)-vcur(i_up+(j-1)*NGP));
			Hp = interpolate(uprev,i_up,j,'v')*(vprev(j*NGP)-vprev((i_up-1)+j*NGP))
				+vprev(i_up+j*NGP)*(vprev(i_up+(j+1)*NGP)-vprev(i_up+(j-1)*NGP));
			Gn = vcur(j*NGP)+vcur((i_up-1)+j*NGP)-4*vcur(i_up+j*NGP)
				+vcur(i_up+(j+1)*NGP)+vcur(i_up+(j-1)*NGP);
			In = 1/DT*vcur(i_up+j*NGP);
			rhs(i_up+j*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		}
		//corner points
		//i=0,j=0
		Hn = interpolate(ucur,0,0,'v')*(vcur(1)-vcur(NGP-1))
			+vcur(0)*(vcur(NGP)-vcur((NGP-1)*NGP));
		Hp = interpolate(uprev,0,0,'v')*(vprev(1)-vprev(NGP-1))
			+vprev(0)*(vprev(NGP)-vprev((NGP-1)*NGP));
		Gn = vcur(1)+vcur(NGP-1)-4*vcur(0)
			+vcur(NGP)+vcur((NGP-1)*NGP);
		In = 1/DT*vcur(0);
		rhs(0) = a*(3*Hn-Hp)+b*Gn+In;
		//i=NGP-1,j=0
		Hn = interpolate(ucur,i_up,0,'v')*(vcur(0)-vcur(i_up-1))
			+vcur(i_up)*(vcur(i_up+NGP)-vcur(i_up+(NGP-1)*NGP));
		Hp = interpolate(uprev,i_up,0,'v')*(vprev(0)-vprev(i_up-1))
			+vprev(i_up)*(vprev(i_up+NGP)-vprev(i_up+(NGP-1)*NGP));
		Gn = vcur(0)+vcur(i_up-1)-4*vcur(i_up)
			+vcur(i_up+NGP)+vcur(i_up+(NGP-1)*NGP);
		In = 1/DT*vcur(i_up);
		rhs(i_up) = a*(3*Hn-Hp)+b*Gn+In;
		//i=0,j=NGP-1
		Hn = interpolate(ucur,0,j_up,'v')*(vcur(1+j_up*NGP)-vcur((NGP-1)+j_up*NGP))
			+vcur(j_up*NGP)*(vcur(0)-vcur((j_up-1)*NGP));
		Hp = interpolate(uprev,0,j_up,'v')*(vprev(1+j_up*NGP)-vprev((NGP-1)+j_up*NGP))
			+vprev(j_up*NGP)*(vprev(0)-vprev((j_up-1)*NGP));
		Gn = vcur(1+j_up*NGP)+vcur((NGP-1)+j_up*NGP)-4*vcur(j_up*NGP)
			+vcur(0)+vcur((j_up-1)*NGP);
		In = 1/DT*vcur(j_up*NGP);
		rhs(j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;
		//i=NGP-1,j=NGP-1
		Hn = interpolate(ucur,i_up,j_up,'v')*(vcur(j_up*NGP)-vcur((i_up-1)+j_up*NGP))
			+vcur(i_up+j_up*NGP)*(vcur(i_up)-vcur(i_up+(j_up-1)*NGP));
		Hp = interpolate(uprev,i_up,j_up,'v')*(vprev(j_up*NGP)-vprev((i_up-1)+j_up*NGP))
			+vprev(i_up+j_up*NGP)*(vprev(i_up)-vprev(i_up+(j_up-1)*NGP));
		Gn = vcur(j_up*NGP)+vcur((i_up-1)+j_up*NGP)-4*vcur(i_up+j_up*NGP)
			+vcur(i_up)+vcur(i_up+(j_up-1)*NGP);
		In = 1/DT*vcur(i_up+j_up*NGP);
		rhs(i_up+j_up*NGP) = a*(3*Hn-Hp)+b*Gn+In;

		
		for(int j=1;j<NGP-1;j++){
			for(int i=1;i<NGP-1;i++){
				Hn = interpolate(ucur,i,j,'v')*(vcur((i+1)+j*NGP)-vcur((i-1)+j*NGP))
					+vcur(i+j*NGP)*(vcur(i+(j+1)*NGP)-vcur(i+(j-1)*NGP));
				Hp = interpolate(uprev,i,j,'v')*(vprev((i+1)+j*NGP)-vprev((i-1)+j*NGP))
					+vprev(i+j*NGP)*(vprev(i+(j+1)*NGP)-vprev(i+(j-1)*NGP));
				Gn = vcur((i+1)+j*NGP)+vcur((i-1)+j*NGP)-4*vcur(i+j*NGP)
					+vcur(i+(j+1)*NGP)+vcur(i+(j-1)*NGP);
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
				if(i!=0){
					I=ucur(i+j*NGP)-ucur((i-1)+j*NGP);
				}
				else{
					I=ucur(i+j*NGP)-ucur((NGP-1)+j*NGP);
				}
				//J
				if(j!=0){
					J=vcur(i+j*NGP)-vcur(i+(j-1)*NGP);
				}
				else{
					J=vcur(i+j*NGP)-vcur(i+(NGP-1)*NGP);
				}
				//update RHS
				rhs(i+j*NGP) = a*(I+J);
			}
		} 
	}
	double wall1 = get_wall_time();
    double cpu1  = get_cpu_time();
	cout << "Update opt" << endl;
    cout << "Wall Time = " << wall1 - wall0 << endl;
    cout << "CPU Time  = " << cpu1  - cpu0  << endl;
}

inline double interpolate(VectorXd &field, int i, int j, char comp){
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
