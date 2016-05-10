#ifndef __TENSOR3_H
#define __TENSOR3_H

/* #include "blitz/tinyvec.h" */
/* #include "blitz/tinyvec-et.h" */
#include "SymmetricMatrix.h"
#include "SquareMatrix.h"
using namespace blitz;

//! Third order tensor
class Tensor3
{

	double Data_[3][3][3];

	public:

	Tensor3() {}

	explicit Tensor3 ( const double scal )  {
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				for (int k=0; k<3; ++k) {
					Data_[i][j][k] = scal;
				}
			}
		}      
	}

	const double& operator()(int i, int j, int k) const
	{
		return Data_[i][j][k];
	}

	double& operator()(int i, int j, int k)
	{
		return Data_[i][j][k];
	}    

	SquareMatrix operator()(int i) const
	{
		SquareMatrix M;
		for(int j=0; j<3; ++j) 
		{
			for(int k=0; k<3; ++k) 
			{
				M(j,k)=Data_[i][j][k];
			}	
		}
		return M;
	}

	const Tensor3& operator-=(const Tensor3& T2) 
	{
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				for (int k=0; k<3; ++k) {
					Data_[i][j][k] -= T2(i,j,k);
				}
			}
		}            
		return *this;      
	}

	const Tensor3& operator*=(const double s) 
	{
		for (int i=0; i<3; ++i) {
			for (int j=0; j<3; ++j) {
				for (int k=0; k<3; ++k) {
					Data_[i][j][k] *= s;
				}
			}
		}            
		return *this;
	}

};

Tensor3 operator-(const Tensor3& T1, const Tensor3& T2); 

//! Third order tensor, symmetric in indices 23
class Symmetric23Tensor3
{

	double Data_[3][6];
	static const int Lookup_[3][3];

	public:

	Symmetric23Tensor3() {}

	explicit Symmetric23Tensor3 ( const double scal )  {
		for (int i=0; i<3; ++i) {
			for (int j=0; j<6; ++j) {
				Data_[i][j] = scal;	    
			}
		}      
	}

	inline double& operator() ( const int i, const int j, const int k )
	{
		return Data_[i][Lookup_[j][k]];
	}

	inline const double& operator() ( const int i, const int j, const int k ) const
	{
		return Data_[i][Lookup_[j][k]];
	}

	SymmetricMatrix operator()(int i) const
	{
		SymmetricMatrix M;
		for(int j=0; j<3; ++j) 
		{
			for(int k=j; k<3; ++k) 
			{
				M(j,k)=Data_[i][Lookup_[j][k]];
			}	
		}
		return M;
	}

	const Symmetric23Tensor3& operator*=(const double s) 
	{
		for (int i=0; i<3; ++i) {
			for (int j=0; j<6; ++j) {
				Data_[i][j] *= s;
			}
		}            
		return *this;
	}
};

ostream& operator<<(ostream& os, const Tensor3& T); 
ostream& operator<<(ostream& os, const Symmetric23Tensor3& T); 



#endif
