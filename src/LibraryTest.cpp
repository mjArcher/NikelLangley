/* -I ~/Libraries/ */
//Matt Archer
//compile as follows: g++ LibraryTest.cpp -I ~/Libraries/ -o test

#include <sstream>
#include <fstream>
#include <iostream>
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Householder>
#include <Eigen/FFT>
#include "ElasticPrimState.h"

#include "SquareTensor3.h" //class for different size square tensors - Templates 

//blitz
/* #include "SquareMatrix.h" */
#include <blitz/tinyvec.h>
#include <blitz/tinyvec-et.h>

/* #include "clapack.h" */

using namespace std;

struct Domain {
	int Ni; //number of x cells
	int GNi; // number of x cells plus ghost cells
	int GC; // number of ghost cells
	int starti; //start of computational domain
	int endi; //end of computational domain
	double Lx; // length of domain
	double dx; // cell size

	Domain(const int a_Ni, const int a_GC, const double a_Lx){};
	Domain(){};
};

void tvmet_L(){
	tvmet::Vector<double, 3> v1(1,2,3);
	tvmet::Vector<double, 3> v2;        // construct the Vector<T,Sz> object at first
	v2 = v1;                     // ... and assign the contents of v1 to v2
	tvmet::Vector<double, 3> v3(v1);    // ... or simple use the copy constructor
	v3 += v1;
	std::cout <<" operator test" <<  v3 << std::endl;
}

void eigenComp()
{
	vector<Eigen::Matrix3d> test(3);
	Eigen::Matrix3d ident = Eigen::Matrix3d::Identity();
	test[0] = ident;
	test[1] = ident*2.;
	test[2] = ident*3.;

	SquareTensor3 tensor(test);
	cout << test << endl;

	tensor *= 3;
	cout << tensor << endl;

	cout << tensor[2](2,2) << endl;
	cout << tensor(1,2,2) << endl;

	cout << "Testing overloaded operators" << endl;

	tensor *= ident;	


	cout << tensor << endl;
	
	Eigen::Matrix3d mtest;
	mtest << 2, 2, 2, 2, 2, 2, 2, 2, 2;

	SquareTensor3 tensor1;
		
	tensor1 = mtest * tensor;

	cout << tensor1 << endl;

	cout << "cross check against simple matrix matrix multiplication" << endl;
	
	Eigen::Matrix3d te2;	
	te2(1,1) = 10;
	te2 = tensor[1];

	cout << te2 * mtest << endl;

	cout << "Matrix addition" << endl;
	Eigen::Matrix3d addTest;

	addTest = te2 + te2;

	cout << addTest << endl;

}

void eigen_L(){
	Eigen::Vector3d v1(1,2,3);
	cout << v1 << endl;
	cout << "get element " << v1[2] << endl;
	cout << "assignment " << endl;
	v1(1) = 4;
	cout << v1 << endl;
	Eigen::Matrix3d m1;
	m1   <<  v1(1), 2, 3,
			 4, 5, 6,
			 7, 8, 9;
	cout << m1 << endl;
	cout << "get element "  << endl;
	cout << m1(1,1) << endl;

	cout << "Alternative assignment " << endl;
	Eigen::Vector3d v2;	
	v2(0) = 1; 
	v2(1) = 2;
	v2(2) = 3;
	cout << v2 << endl;
	v2+=v1;
	cout << "Operator test " << v2 << endl;

	vector<Eigen::Matrix3d> vtest(3);
	Eigen::Matrix3d mtest;
	mtest << 1,2,3,4,5,6,7,8,9;
	vtest[0] = mtest;
	cout << "Indexing " << endl;
	cout << vtest[0](1,1) << endl;	
	cout << "identity " << endl;
	
	Eigen::Matrix3d identity = Eigen::Matrix3d::Identity();
	cout << identity << endl;	
	
	/* vector<Eigen::Matrix3d> itest(3); */
	/* itest[0] = v1; */
	
}

void blitz_L(){
	blitz::TinyVector<double, 3> v1(1,2,3);
	cout << v1 << endl;
	//compute the dot of a matrix times a vector
	blitz::TinyVector<ElasticPrimState, 14> primelastic;
  cout << primelastic << endl;
}

void Lapack_L()
{
  // worth testing this with euler material

}

void boost(){


}

int main(int argc, char* argv[])
{		
	cout << "Compute matrix and vector operations with tvmet" << endl;
	tvmet_L();
	/* cout << "Compute matrix and vector operations with eigen" << endl; */
	/* eigen_L(); */
	cout << "Compute matrix and vector operations with blitz" << endl;
	blitz_L();
	/* cout << "Compute matrix and vector operations with eigen" << endl; */
	/* eigen_L(); */

	cout << "Domain test" << endl;
	Domain dom;
	dom = Domain(100, 3, 0.1);

	Eigen::VectorXd test;
	test = Eigen::VectorXd(13);
	cout << test << endl;

	cout << "More advanced eigen tests" << endl;
	
	eigenComp();
}
