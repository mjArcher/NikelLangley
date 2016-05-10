#include "SolidSystem.h"

//equation of state and system tests
//
//

// System Derivatives to test 
//dsigma/dG
//drho/dF
//dStress/dFe

//Equation of state derivatives to test
//depsi_dI
//d2epsi_dI_dI
//


int main(void) 
{
	//initialize example states
	//
	
  ElasticPrimState pri(SquareMatrix( 1.01, 0.01, 0, -0.1, 0.95, 0, 0, 0, 0.99),
		       TinyVector<double, 3>(0, 0, 0),
		       10, 0.6);
  SquareMatrix Fp(1.2,0.1,0.0,0.02,1.02,0,0,0,0.9);
  Fp/=cbrt(Fp.Determinant());

  ElasticPrimState priB(SquareMatrix( 1, 0.1, 0, -0.2, 1, 0, 0, 0, 1),
		       TinyVector<double, 3>(0, 0, 0),
		       0, 0);


	Matrix3d FL, FR;
	Vector3d uL, uR;
	double SL, SR;
	FL <<	 1.01, 0.01, 0, -0.1, 0.95, 0, 0, 0, 0.99;
	uL << 0, 0, 0;
	SL = 10.;

	FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
	uR << 0, 0, 0;
	SR = 0.;


}
