#ifndef ELASTICPRIMSTATE_H
#define ELASTICPRIMSTATE_H

#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip> //setw

#include "SquareTensor3.h"
 
class ElasticPrimState
{
	//velocity, deformation gradient, entropy
	public: 
		ElasticPrimState(Eigen::Vector3d u, Eigen::Matrix3d F, double S);

		ElasticPrimState();

		~ElasticPrimState();	

    //accessors
		Eigen::Vector3d u_() const;

		Eigen::Matrix3d F_() const;
		
		double S_() const;

    //setters
    void u(const Eigen::Vector3d&);

    void F(const Eigen::Matrix3d&);

    void S(const double&);

		Eigen::VectorXd getStateVector() const {return v;};
	
		//! Acoustic tensor derivatives

		SquareTensor3 dI_dG(const Eigen::Matrix3d&, const Eigen::Vector3d&) const;

		SquareTensor3 dI_dF(const Eigen::Matrix3d&, const Eigen::Vector3d&) const;

		Eigen::Matrix3d dsigma_dG(int, const Eigen::Vector3d&, const Eigen::Matrix3d&, double rho) const;

		Eigen::Matrix3d dG_dF(const Eigen::Matrix3d& G, const Eigen::Matrix3d& F, int j, int m) const;
		//! Overloaded operators
		
		double& operator[] (int x){return v[x];}

		ElasticPrimState& operator+=(const ElasticPrimState& consState);

		ElasticPrimState& operator-=(const ElasticPrimState& consState);

		ElasticPrimState& operator*=(const double& scalar);

		ElasticPrimState& operator/=(const double& scalar);
		
/* ElasticPrimState operator=(ElasticPrimState); */	

		friend std::ostream& operator<<(std::ostream&, const ElasticPrimState&);

		static const int e_size = 13;
	private:
		/* Initial states */
		
	//! Solution vector
		Eigen::VectorXd v; //Check tomorrow: Does this have equality built in?
};

		ElasticPrimState operator+(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

		ElasticPrimState operator-(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

		ElasticPrimState operator*(const ElasticPrimState& lhs, const double& s);

		ElasticPrimState operator*(const double& s, const ElasticPrimState& rhs);

		ElasticPrimState operator/(const ElasticPrimState& lhs, const double& s);

    bool operator==(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

    bool operator!=(const ElasticPrimState& lhs, const ElasticPrimState& rhs);

    


#endif
