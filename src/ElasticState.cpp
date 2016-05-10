#include "ElasticState.h"

using namespace Eigen;
using namespace std;


//mom, rhoF, rhoE
//! Constructor
ElasticState::ElasticState(Vector3d mom, Matrix3d rhoF, double rhoE)
{
	v = VectorXd(ElasticState::e_size);
	v << mom(0), mom(1), mom(2),
		rhoF(0,0), rhoF(0,1), rhoF(0,2),
		rhoF(1,0), rhoF(1,1), rhoF(1,2),
		rhoF(2,0), rhoF(2,1), rhoF(2,2), rhoE;
}

	/* v << mom(0), mom(1), v[2] = mom(2), */
	/* v[3] = rhoF(0,0), v[4] = rhoF(0,1), v[5] = rhoF(0,2), */
	/* v[6] = rhoF(1,0), v[7] = rhoF(1,1), v[8] = rhoF(1,2), */
	/* v[9] = rhoF(2,0), v[10] = rhoF(2,1), v[11] = rhoF(2,2), */
	/* v[12] = rhoE; */

ElasticState::ElasticState(){
	v = VectorXd(ElasticState::e_size);
}
//destructor
ElasticState::~ElasticState(){}

//! Momentum vector
Vector3d ElasticState::mom() const
{
	Vector3d mom;
	mom << v[0], v[1], v[2];
	return mom;
}

//! Rho * deformation gradient matrix
Matrix3d ElasticState::rhoF() const
{
	Matrix3d rhoF;
	rhoF << v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11];
	return rhoF;
}

//! rhoE double
double ElasticState::rhoE() const
{
	return v(12);
}

//overload operators
std::ostream& operator<<(std::ostream& os, const ElasticState& param) 
{
	os.precision(3);
	Vector3d m = param.mom();
	Matrix3d rhoFe = param.rhoF();
	double rhoEn = param.rhoE();
	
	
	/* os << setw(3) << " |" << setw(sp) << m(0) << setw(3) << " |" << */ 
	
  os << setprecision(4) << "\t{ Mom   = ( " << setw(9) << m[0] << "  " << setw(9) << m[1] << "  " << setw(9) << m[2] << " ), rhoE  = " << setw(9) << rhoEn << "\n";
	os << "\t          [ " << setw(9) << rhoFe(0,0) << "  " << setw(9) << rhoFe(0,1) << "  " << setw(9) << rhoFe(0,2) << " ]\n";
	os << "\t  rhoFe = [ " << setw(9) << rhoFe(1,0) << "  " << setw(9) << rhoFe(1,1) << "  " << setw(9) << rhoFe(1,2) << " ]\n";
	os << "\t          [ " << setw(9) << rhoFe(2,0) << "  " << setw(9) << rhoFe(2,1) << "  " << setw(9) << rhoFe(2,2) << " ]\n";


	/* os.precision(3); */
	/* for (int i = 0; i < ElasticState::e_size; i++) */
	/* { */
	/* 	os << setw(3) << " |" << setw(7) << param[i]; */
	/* } */	
	/* os << setw(7) << " | "; */
	return os;
}

ElasticState& ElasticState::operator+=(const ElasticState& consState)
{
	// actual addition of rhs to *this
	v += consState.getStateVector();
	return *this;
}

ElasticState& ElasticState::operator-=(const ElasticState& consState)
{
	// actual addition of rhs to *this
	v -= consState.getStateVector();
	return *this;
}

ElasticState& ElasticState::operator*=(const double& scalar)
{
	// actual addition of rhs to *this
	v *= scalar;
	return *this;
}

ElasticState& ElasticState::operator/=(const double& scalar)
{
	// actual addition of rhs to *this
	v /= scalar;
	return *this;
}
//const safety 
ElasticState operator+(const ElasticState& lhs, const ElasticState& rhs)
{
	ElasticState C(lhs);	
  C += rhs;
  return C;
}
//const safety
ElasticState operator-(const ElasticState& lhs, const ElasticState& rhs)
{
	ElasticState C(lhs);	
  C -= rhs;
  return C;
}

ElasticState operator*(const ElasticState& lhs, const double& s)
{
  return s*lhs;
}

ElasticState operator*(const double& s, const ElasticState& rhs)
{
	ElasticState C(rhs);
  C *= s;
  return C;
}

ElasticState operator/(const ElasticState& lhs, const double& s)
{
  return lhs*(1.0/s);
}


/* class X { */
/*   X& operator+=(const X& rhs) */
/*   { */
/*     // actual addition of rhs to *this */
/*     return *this; */
/*   } */
/* }; */
/* inline X operator+(X lhs, const X& rhs) */
/* { */
/*   lhs += rhs; */
/*   return lhs; */
/* } */


//Will eventaully remove this function
//Requires the first few invariants and entropy

/* void ElasticState::prepareOutputs() */
/* { */
/* 	Matrix3d rhoF; */
/* 	ElasticState temp = *this; */	
/* 	rhoF <<  v[3], v[4], v[5], */
/* 			 v[6], v[7], v[8], */
/* 			 v[9], v[10], v[11]; */
/* 	/1* rho = pow(rhoF.determinant()/2., 2.); *1/ */
/* 	rho = sqrt(rhoF.determinant()/rho_0); */
/* 	//deformation gradient */
/* 	Matrix3d F = rhoF/rho; */
/* 	// Calculate finger tensor */
/* 	Matrix3d G = strainTensor(F); */

/* 	//calculate invariants */
/* 	/1* cout << B_0 << endl; *1/ */
/* 	invariants(G); */
/* 	dInvariants_G(); */

/* 	//Stress */
/* 	stress(G); */
/* } */

/* double ElasticState::rho_() */
/* { */
/* 	/1* rho = sqrt(rhoF.determinant()/rho_0); *1/ */
/* 	/1* cout << rho << endl; *1/ */
/* 	return rho; */
/* } */

// throw error if index is greater than certain value
