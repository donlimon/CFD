//Functions for generating/updating coefficient matrices
//file: matrixgen.cpp
//author: Michael Stumpf

#include <Eigen/SparseCore>
#include <iostream> //for debugging
#include <settings.hpp>

using namespace std;
using namespace Eigen;

SparseMatrix<double,RowMajor> createPredCoeffMat(){
	SparseMatrix<double,RowMajor> mat(NGP*NGP,NGP*NGP);
	double 	a=1/DT+2/(RE*DX*DX), 
			b=-1/(2*RE*DX*DX);
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
