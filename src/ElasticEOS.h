#ifndef ELASTICEOS_H
#define ELASTICEOS_H

#include <cstdlib>
#include <iostream>
#include <sstream>


#include <Eigen/Dense>
#include <Eigen/Core>
// this is the Romenski equation of state
//

class ElasticEOS {

	public:

		ElasticEOS(const std::string&);
		
		//! Default constuctor. Careful!
		ElasticEOS();

		double entropy(const Eigen::Vector3d invariants, double internalEnergy) const;	

		double internalEnergy(const Eigen::Vector3d invariants, double Entropy) const; //invariants, 

		void checkEosConstants() const;

		double soundSpeed() const;

		Eigen::Vector3d depsi_dI(const Eigen::Vector3d I, double S) const;

		Eigen::Matrix3d depsi_dI_dI(const Eigen::Vector3d& I, double S) const;

    Eigen::Vector3d depsi_dI_dS(const Eigen::Vector3d& I, double S) const;

		~ElasticEOS();

		double rho0; //this is public because it is required by  the density function in the system class 
	private:
		
		//! Equation of state parameters 	

		//! Density

		//! Heat capacity at constant volume
		double cv;

		//! Temperature
		double T0;

		//! squared bulk shear wave
		double B0;	

		//! squared bulk speed of sound
		double K0;

		//! Constants characterising the non-linear dependence of sound speeds and
		//! temperature on the mass density
		double alpha;

		double beta;

		double gamma;

};
#endif
