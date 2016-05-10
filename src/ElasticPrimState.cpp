#include "ElasticPrimState.h"

using namespace Eigen;
using namespace std;

//there is a redundancy here. Note the accessor and setter functions.
//
//1. Calc U,
ElasticPrimState::ElasticPrimState(Vector3d u, Matrix3d F, double S)
{
	v = VectorXd(ElasticPrimState::e_size);

	v << u(0), u(1), u(2),
		F(0,0), F(0,1), F(0,2),
		F(1,0), F(1,1), F(1,2),
		F(2,0), F(2,1), F(2,2),	S; /* cout << "calc Rho" << rhoNew << endl; */
}

// default constructor

ElasticPrimState::ElasticPrimState(){
	v = VectorXd::Zero(ElasticPrimState::e_size);
}

ElasticPrimState::~ElasticPrimState(){}

void ElasticPrimState::u(const Vector3d& vel)
{
  v(0) = vel(0), v(1) = vel(1), v(2) = vel(2); 
}

void ElasticPrimState::F(const Matrix3d& mat)
{
  v(3) = mat(0,0), v(4) = mat(0,1), v(5) = mat(0,2), 
  v(6) = mat(1,0), v(7) = mat(1,1), v(8) = mat(1,2),
  v(9) = mat(2,0), v(10) = mat(2,1), v(11) = mat(2,2);
}

void ElasticPrimState::S(const double& S)
{
  v(12) = S; 
}

//! velocity vector
Vector3d ElasticPrimState::u_() const
{
	Vector3d vel;
	vel << v(0), v(1), v(2);
	return vel;
}

//! Rho * deformation gradient matrix
Matrix3d ElasticPrimState::F_() const
{
	Matrix3d rhoF;
	rhoF << v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11);
	return rhoF;
}

//! rhoE double
double ElasticPrimState::S_() const
{
	return v(12);
}

std::ostream& operator<<(std::ostream& os, const ElasticPrimState& param)
{
  os.precision(7);
	Vector3d u = param.u_();
	Matrix3d Fe = param.F_();
	double S = param.S_();
  os << "{ Vel = "   << setw(14) << u << ", S  = "   << setw(14) << S << "\n";
  os << "        [ "  << setw(14) << Fe(0,0) << "\t"  << setw(14) << Fe(0,1) << "\t"  << setw(14) << Fe(0,2) << " ]\n";
  os << "   Fe = [ "  << setw(14) << Fe(1,0) << "\t"  << setw(14) << Fe(1,1) << "\t"  << setw(14) << Fe(1,2) << " ]\n";
  os << "        [ "  << setw(14) << Fe(2,0) << "\t"  << setw(14) << Fe(2,1) << "\t"  << setw(14) << Fe(2,2) << " ]}\n";
	/* os.precision(3); */
	/* for (int i = 0; i < ElasticPrimState::e_size; i++) */
	/* { */
	/* 	os << setw(3) << " |" << setw(7) << param[i]; */
	/* } */	
	/* os << setw(7) << " | "; */
	return os;
}

ElasticPrimState& ElasticPrimState::operator+=(const ElasticPrimState& primState)
{
	// actual addition of rhs to *this
	v += primState.getStateVector();
	return *this;
}

ElasticPrimState& ElasticPrimState::operator-=(const ElasticPrimState& primState)
{
	// actual addition of rhs to *this
	v -= primState.getStateVector();
	return *this;
}

ElasticPrimState& ElasticPrimState::operator*=(const double& scalar)
{
	// actual addition of rhs to *this
	v *= scalar;
	return *this;
}

ElasticPrimState& ElasticPrimState::operator/=(const double& scalar)
{
	// actual addition of rhs to *this
	v /= scalar;
	return *this;
}

// const safety
ElasticPrimState operator+(const ElasticPrimState& lhs, const ElasticPrimState& rhs)
{
	ElasticPrimState C(lhs);
  C += rhs;
  return C;
}

//const safety
ElasticPrimState operator-(const ElasticPrimState& lhs, const ElasticPrimState& rhs)
{
	ElasticPrimState C(lhs);
  C -= rhs;
  return C;
}

ElasticPrimState operator*(const ElasticPrimState& lhs, const double& s)
{
  return s*lhs;
}

ElasticPrimState operator*(const double& s, const ElasticPrimState& rhs)
{
	ElasticPrimState C(rhs);
	C *= s;
  return C;
}

ElasticPrimState operator/(const ElasticPrimState& lhs, const double& s)
{
  return lhs/s;
} 

//this needs some thought: 
bool operator==(const ElasticPrimState& p1, const ElasticPrimState& p2) 
{
  VectorXd v1 = p1.getStateVector();
  VectorXd v2 = p2.getStateVector();
  /* std::cout << v1.isApprox(v2) << std::endl; */
  /* std::cout.precision(15); */
  /* std::cout << v1 << " " << v2 << std::endl; */
  return v1.isApprox(v2, 1e-10);
  /* std::cout << p1.S_() << " " << p2.S_() << std::endl; */
  /* return p1.u_().isApprox(p2.u_())&&p1.F_().isApprox(p2.F_())&&p1.S_().isApprox(p2.S_()); */
}

bool operator!=(const ElasticPrimState& p1, const ElasticPrimState& p2)
{
  return !(p1==p2);
}

SquareTensor3 ElasticPrimState::dI_dG(const Matrix3d& G, const Vector3d& inv) const
{
	vector<Eigen::Matrix3d> dI_dG(3); 	
	Matrix3d ID = Matrix3d::Identity();	
	dI_dG[0] = ID; 
	dI_dG[1] = inv(0) * ID - G.transpose(); 
	dI_dG[2] = inv(2) * G.transpose().inverse();
	SquareTensor3 tensor(dI_dG);
	return tensor;
}

//final derivative
SquareTensor3 ElasticPrimState::dI_dF(const Matrix3d& G, const Vector3d& inv) const 
{
	vector<Matrix3d> dI_dF(3);
	Matrix3d Fe = F_();
	Matrix3d m2FinvT = Fe.inverse().transpose();
	Matrix3d ident = Matrix3d::Identity();

	dI_dF[0] = -2*G*m2FinvT;
	/* dI_dF[1] = -2*(G - G.trace()*ident)*m2FinvT; */
	dI_dF[1] = -2*G*(inv(0)*ident - G.transpose())*m2FinvT;
	dI_dF[2] = -2*inv(2)*m2FinvT;

	SquareTensor3 tensor(dI_dF);
	return tensor;
}



//assuming that i is always 1 for 1D 
//nneds to be changed for higher dimensions
Matrix3d ElasticPrimState::dsigma_dG(int k, const Vector3d& de_dI, const Matrix3d& G, double rho) const
{

	Matrix3d ds_dG;

	if(k == 0){
		ds_dG(0,0) = -(de_dI[0] + (G(1,1) + G(2,2))*de_dI[1] + (-G(1,2)*G(2,1) + G(1,1)*G(2,2))*de_dI[2]);

		ds_dG(0,1) = G(1,0)*de_dI[1]-(G(1,2)*G(2,0) - G(1,0)*G(2,2))*de_dI[2];

		ds_dG(0,2) = G(2,0)*de_dI[1]-(-G(1,1)*G(2,0) + G(1,0)*G(2,1))*de_dI[2];

		ds_dG(1,0) = G(0,1)*de_dI[1]-(G(0,2)*G(2,1) - G(0,1)*G(2,2))*de_dI[2];

		ds_dG(1,1) = -(G(0,0)*de_dI[1]+(-G(0,2)*G(2,0) + G(0,0)*G(2,2))*de_dI[2]);

		ds_dG(1,2) = -(G(0,1)*G(2,0) - G(0,0)*G(2,1))*de_dI[2];

		ds_dG(2,0) = G(0,2)*de_dI[1]-(-G(0,2)*G(1,1) + G(0,1)*G(1,2))*de_dI[2];

		ds_dG(2,1) = -(G(0,2)*G(1,0) - G(0,0)*G(1,2))*de_dI[2];

		ds_dG(2,2) = -(G(0,0)*de_dI[1]+(-G(0,1)*G(1,0) + G(0,0)*G(1,1))*de_dI[2]);
	}
	else if(k == 1){

		ds_dG(0,0) = 0;

		ds_dG(0,1) = -(de_dI[0] + G(2,2)*de_dI[1]);	

		ds_dG(0,2) = G(2,1)*de_dI[1];

		ds_dG(1,0) = 0;

		ds_dG(1,1) = 0;

		ds_dG(1,2) = 0;

		ds_dG(2,0) = 0;

		ds_dG(2,1) = G(0,2)*de_dI[1];

		ds_dG(2,2) = -G(0,1)*de_dI[1];
	}
	else if(k == 2){

		ds_dG(0,0) = 0;

		ds_dG(0,1) = G(1,2)*de_dI[1];

		ds_dG(0,2) = -(de_dI[0] + G(1,1) * de_dI[1]);

		ds_dG(1,0) = 0;

		ds_dG(1,1) = -G(0,2)*de_dI[1];

		ds_dG(1,2) = G(0,1)*de_dI[1];

		ds_dG(2,0) = 0;

		ds_dG(2,1) = 0;

		ds_dG(2,2) = 0;

	}
	ds_dG *= 2.*rho;
	return ds_dG;
  //more derivatives needed : check using 
}

//this is inefficient - may need to be changed
Matrix3d ElasticPrimState::dG_dF(const Matrix3d& G, const Matrix3d& F, int j, int m) const
{
	Matrix3d dGdF;	
	const Matrix3d Finv = F.inverse();
	
	for(int q = 0; q < 3; q++){
		for(int r = 0; r < 3; r++){
			dGdF(q,r) = -G(q,j) * Finv(m,r) - G(j,r)*Finv(m,q);
		}
	}
	return dGdF;
}

/* TinyVector<SymmetricMatrix, 3> ElasticPrimState::dI_dG(const SymmetricMatrix& G, const Invariants& inv) const { */ 
/*   const SymmetricMatrix ID = SymmetricMatrix::Identity; */
/*   return TinyVector<SymmetricMatrix,3>(ID, inv.I * ID - G.Transpose(), G.DetInverse().Transpose()); */
/* } */

//array containing 3 matrices 

//test functions - eigen library
/* Tensor<double, 4> dsigma_dG() */
/* { */
/* 	Tensor<double, 4>	t(10, 10, 10, 10); */	
	
/* 	return t; */
/* } */

/* Tensor<double, 4> dG_dF() */
/* { */
/* 	Tensor<double, 4>	t(10, 10, 10, 10); */	
/* 	return t; */
/* } */
