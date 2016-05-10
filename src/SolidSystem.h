#ifndef SOLIDSYSTEM_H
#define SOLIDSYSTEM_H

#include <Eigen/Dense>
#include <Eigen/Core>
#include <iostream>
#include <fstream>
#include <sstream>

#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "ElasticEOS.h"
#include "Utils.h"
#include "SquareTensor3.h"

//convert to class

typedef Eigen::Matrix<double, 13, 1> Vector13d;
//possible linker error?
struct eigen { 
  double lambda;
  ElasticPrimState eigenvecR;
  ElasticPrimState eigenvecL;
  
  //default constructor
  eigen(){}

  eigen(const double val, Vector13d vecL, Vector13d vecR){
    lambda = val;
    eigenvecL = convert(vecL);
    eigenvecR = convert(vecR); //normalised 
  }

  ElasticPrimState convert(Vector13d vecR)
  {
    ElasticPrimState temp;
    for(int i = 0; i < 13; i++)
    {
      temp[i] = vecR[i];
    }
    return temp;
  }

  //descending order 
  bool operator>(eigen const &other) const { 
    return lambda > other.lambda;
  }

  //ascending order 
  bool operator<(eigen const &other) const { 
    return lambda < other.lambda;
  }
};

class System {

//store invariants
// derivatives of invariants with respec to the components of the strain tensor -- here or EOS

	public:
		System(const ElasticEOS Eos);
		~System();

		ElasticEOS Eos;

		//! Convert to primitive state
		ElasticPrimState conservativeToPrimitive(const ElasticState& consState) const;

		//! Convert to conservative state
		ElasticState primitiveToConservative(const ElasticPrimState& primState) const;

    //! Calculate density from matrix
    double Density(const Eigen::Matrix3d& F) const;

		//! Calculate density from primitive state
		double Density(const ElasticPrimState& primState) const;
		
		//! Calculate density from primitive state
		double Density(const ElasticState& consState) const;

		//! Determine finger tensor
		Eigen::Matrix3d strainTensor(const Eigen::Matrix3d& F) const;

		//! Vector of invariants computed from finger tensor
    Eigen::Vector3d getInvariants(const Eigen::Matrix3d& F) const;

    ElasticState flux(const ElasticPrimState&) const;

		ElasticState flux(const ElasticState&) const;
		
    ElasticState flux(const ElasticState& consState, const ElasticPrimState& primState) const;

		//! Symmetric matrix of stresses.
		Eigen::Matrix3d stress(const ElasticPrimState&) const;

		//! Compute the eigen decomposition of the acoustic tensor - return a vector of eigenvalues
		//13 element vectors - a square matrix containing 13 x 13 elements
		Eigen::VectorXd stateEigenDecompose(const ElasticPrimState&, const int, std::vector<ElasticPrimState>&, std::vector<ElasticPrimState>&) const;

		//! Maximum wave speeds (correspond to maximum eigenvalue of acoustic tensor) required for the calculation of stable timestep
		double maxWaveSpeeds() const;

		//! Compute the acoustic tensor 
		Eigen::Matrix3d AcousticTensor(const ElasticPrimState& pri, const int dirn) const;

    ElasticState godunovFlux(const ElasticState& qL, const ElasticState& qR) const;
    
    Eigen::Vector3d B(const ElasticPrimState&) const;

    std::vector<eigen> stateEigenDecompose_A(const ElasticPrimState& pL, const ElasticPrimState& pR, const int dirn) const;

    Eigen::VectorXd stateEigenDecompose_A_old(const ElasticPrimState& pW, const int dirn, std::vector<ElasticPrimState>& Le, std::vector<ElasticPrimState>& Re) const;

    Eigen::VectorXd stateEigenDecompose_Omega(const ElasticPrimState& pW, const int dirn, std::vector<ElasticPrimState>& Le, std::vector<ElasticPrimState>& Re) const; 

		SquareTensor3 dstress_dF(const ElasticPrimState&, const Eigen::Matrix3d& G, const Eigen::Vector3d& inv, const int dirn) const;

		std::vector<Eigen::Matrix3d> dep_dF(const SquareTensor3, const Eigen::Matrix3d) const;

	  ElasticPrimState godunovState(const ElasticPrimState& pL, const ElasticPrimState& pR, const int dirn) const;

    double getMaxWaveSpeed(const ElasticPrimState& pL, const int dirn) const;

    double dotState(ElasticPrimState& state1, ElasticPrimState& state2) const;

  private:


    /* ElasticPrimState godunovState(const ElasticPrimState pL, const ElasticPrimState pR); */
    /* //prototype may need to be modified */ 


    /* ElasticState godunovFlux(const ElasticPrimState& pL, const ElasticPrimState& pR); */

		//! Derivatives required for the computation of the acoustic tensor
		//! Derivative of stress (sigma) w.r.t density
		/* Eigen::Tensor<double, 4> dstressrho(const ) const; */

		/* //! Derivative of density w.r.t deformation gradient F */
		/* Eigen::Tensor<double, 4> rhodF(const ) const; */

		/* //! Derivative of stress w.r.t the derivative of the internal energy density w.r.t the invariants I_p */
		/* Eigen::Tensor<double, 4> dsdepsI(const ) const; */

		/* //! Derivative of epsI (defined as the derivative of the internal energy density w.r.t each invariant) w.r.t deformation gradient F. */
		/* Eigen::Tensor<double, 4> depsIdF(const ) const; */

		
		//! Derivatives required for the computation of the acoustic tensor

		//! Derivative of stress (sigma) w.r.t density
		/* Tensor4 dstressrho(const ) const; */

		/* //! Derivative of density w.r.t deformation gradient F */
		/* rhodF(const ) const; *1/ */

		/* //! Derivative of stress w.r.t the derivative of the internal energy density w.r.t the invariants I_p */
		/* dsdepsI(const ) const; */

		/* //! Derivative of epsI (defined as the derivative of the internal energy density w.r.t each invariant) w.r.t deformation gradient F. */
		 /* depsIdF(const ) const; */

		
	
		
/* string getDirName(const string&, const Domain&); */

};
#endif
