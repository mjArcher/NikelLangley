#ifndef SQUARETENSOR3_H
#define SQUARETENSOR3_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

class SquareTensor3
{
	private:

		std::vector<Eigen::Matrix3d> tensor;

	public:

		SquareTensor3(){}

		SquareTensor3(std::vector<Eigen::Matrix3d> tensorInit){
			tensor = tensorInit;
		}

		const double& operator()(int i, int j, int k) const
		{
			return tensor[i](j,k);
		}
		
		const Eigen::Matrix3d& operator[](int i) const
		{
			return tensor[i];
		}

		const SquareTensor3& operator*=(const double s) 
		{
			for (int i=0; i<3; ++i) {
				tensor[i] *= s;
			}            
			return *this;
		}
		
		const SquareTensor3& operator*=(const Eigen::Matrix3d mat2d)
		{
			for (int i=0; i<3; ++i) {
				tensor[i] *= mat2d; //overloaded in eigen
			}
			return *this;
		}

		const SquareTensor3& operator-=(const SquareTensor3& T2) 
		{
			for (int i=0; i<3; ++i) {
				tensor[i] -= T2[i];
			}            
			return *this;      
		}

};

std::ostream& operator<<(std::ostream& os, const SquareTensor3& T); 

SquareTensor3 operator*(const SquareTensor3& lhs, const Eigen::Matrix3d& s);

SquareTensor3 operator*(const Eigen::Matrix3d& s, const SquareTensor3& rhs);



#endif
