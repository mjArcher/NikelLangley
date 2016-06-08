#ifndef ELASTICSTATE_H
#define ELASTICSTATE_H

#include <math.h>
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
/* #include <Eigen/CXX11/Tensor> */
/* #include <Eigen/FFT> */
/* #include <Eigen/CXX11/Tensor> */
/* #include <unsupported/Eigen/CXX11/Tensor.h> */
/* #include <tvmet/Matrix.h> */
#include <fstream>
#include <vector>
#include <iomanip> 

class ElasticPrimState;

class ElasticState
{
	public:
		ElasticState(Eigen::Vector3d mom, Eigen::Matrix3d rhoFe, double rhoE);

		ElasticState();

		~ElasticState();	
		
		Eigen::Vector3d mom() const;

		Eigen::Matrix3d rhoF() const;
	
		double rhoE() const;

		Eigen::VectorXd getStateVector() const {return v;};

		//! Overloaded operators
		double& operator[] (int x) {return v[x];}

		ElasticState& operator+=(const ElasticState&);

		ElasticState& operator-=(const ElasticState&);

		ElasticState& operator*=(const double&);

		ElasticState& operator/=(const double&);

		/* ElasticState operator=(ElasticState); */	

		friend std::ostream& operator<<(std::ostream&, const ElasticState&);

		static const unsigned int e_size = 13;

	private:
		/* Initial states */

		Eigen::VectorXd v;
};

		ElasticState operator+(const ElasticState& lhs, const ElasticState& rhs);

		ElasticState operator-(const ElasticState& lhs, const ElasticState& rhs);

		ElasticState operator*(const ElasticState& lhs, const double& s);

		ElasticState operator*(const double& s, const ElasticState& rhs);

		ElasticState operator/(const ElasticState& lhs, const double& s);

/* inline ElasticState::ElasticState(Eigen::Vector3d u_i, Eigen::Matrix3d F_i, double rho_i, double S_0, double e_0):rho_0(rho_i), S(S_0) */
/* { */
/* 	u_0 = u_i; */
/* 	F_0 = F_i; */
/* 	Eigen::VectorXd temp(e_size); */
/* 	v = temp; */
/* 	rho = rho_0/(F_0.determinant()); */ 
/* 	/1* cout << "rho Calculated " << rho << " mod u " << mod_u() << endl; *1/ */
/* 	E_0 = e_0 + (u_0(0)*u_0(0)+u_0(1)*u_0(1)+u_0(2)*u_0(2))/2.; */
/* 	E = E_0; */
/* 	Eigen::Vector3d mom = rho*u_0; */
/*     Eigen::Matrix3d rhoF = rho*F_0; */
/* 	v[0] = mom(0), v[1] = mom(1), v[2] = mom(2), */
/* 	v[3] = rhoF(0,0), v[4] = rhoF(0,1), v[5] = rhoF(0,2), */
/* 	v[6] = rhoF(1,0), v[7] = rhoF(1,1), v[8] = rhoF(1,2), */
/* 	v[9] = rhoF(2,0), v[10] = rhoF(2,1), v[11] = rhoF(2,2), */
/* 	v[12] = rho*E_0; */
/* 	/1* cout << "calc Rho" << rhoNew << endl; *1/ */
/* } */

/* inline double ElasticState::mod_u() */
/* { */
/* 	return u_(0)*u_(0)+u_(1)*u_(1)+u_(2)*u_(2); */
/* } */

/* inline void ElasticState::initialStates(std::ostream& os) */
/* { */
/* 	os <<   "Initial States\n " << */ 
/* 			"Density Initial\t" << rho_0 <<"\n" */ 
/* 			"Energy Initial\t"	<< E_0 << "\n" */
/* 			"Velocity \t" 		<< u_0 << "\n"; */
/* } */
/* // calc /1* //default constructor *1/ */
/* inline ElasticState::ElasticState(){ */
/* 	Eigen::VectorXd temp(e_size); */
/* 	v = temp; */
/* } */
/* //destructor */
/* inline ElasticState::~ElasticState(){} */

/* inline ElasticState ElasticState::F() */
/* { */
/* 	//Calculate deformation gradient F(from conserved using density) */
/* 	//Deformation gradient: */
/* 	//first determine density */
/* 	Eigen::Matrix3d rhoF; */

/* 	Eigen::Matrix3d F = rhoF/rho; */
/* 	Eigen::Matrix3d G = strainTensor(F); */
/* 	//calculate invariants */
/* 	invariants(G); */
/* 	dInvariants_G(); */
/* 	ElasticState Fl; //flux */
/* 	//Stress */
/* 	stress(G); */
/* 	Eigen::Vector3d rhou; */
/* 	rhou << v[0], v[1], v[2]; */
/* 	Eigen::Vector3d u = rhou/rho; */
/* 	Fl[0]  = v[0]*u[0] - sigma(0,0); */
/* 	Fl[1]  = v[1]*u[0] - sigma(1,0); */
/* 	Fl[2]  = v[2]*u[0] - sigma(2,0); */
/* 	Fl[3]  = 0; */ 
/* 	Fl[4]  = 0; */ 
/* 	Fl[5]  = 0; */ 
/* 	Fl[6]  = rho*F(1,0)*u(0) - rho*F(0,0)*u(1); */
/* 	Fl[7]  = rho*F(1,1)*u(0) - rho*F(0,1)*u(1); */
/* 	Fl[8]  = rho*F(1,2)*u(0) - rho*F(0,2)*u(1); */
/* 	Fl[9]  = rho*F(2,0)*u(0) - rho*F(0,0)*u(2); */
/* 	Fl[10] = rho*F(2,1)*u(0) - rho*F(0,1)*u(2); */
/* 	Fl[11] = rho*F(2,2)*u(0) - rho*F(0,2)*u(2); */
/* 	Fl[12] = rho*u(0)*E_() - u(0)*sigma(0,0)-u(1)*sigma(0,1)-u(1)*sigma(0,2); */
/* 	return Fl; */
/* } */

/* Eigen::Matrix3d ElasticState::strainTensor(Matrix3d F) */
/* { */
/* 	Eigen::Matrix3d G = (F.inverse().transpose())*(F.inverse()); */	
/* 	return G; */
/* } */

/* void ElasticState::invariants(Eigen::Matrix3d G) */
/* { */	
/* 	I[0] = G.trace(); */
/* 	I[1] = 0.5*(pow(G.trace(), 2)-(G*G).trace()); */
/* 	I[2] = G.determinant(); */
/* 	dInvariants_G(); */
/* } */

/* void ElasticState::dInvariants_G() */
/* { */
/* 	/1* double alpha = 1.0, b_0 = 2.1, beta = 3.0; *1/ */
/* 	/1* double B_0 = b_0*b_0; *1/ */ 
/* 	/1* //dElasticState/DI_3, dW/dI_1, dW/dI_2, dW/dI_3 *1/ */

/* 	double I3a = pow(I[2], alpha/2.); */
/* 	double I3b = pow(I[2], beta/2.); */
/* 	double I3ad = pow(I[2], alpha/2. - 1.); */
/* 	double I3bd = pow(I[2], beta/2. - 1.); */
/* 	double I3gd = pow(I[2], gamma/2. - 1.); */
/* 	I[3] = (B_0/3.)*I3b*I[0]; */
/* 	I[4] = (-B_0/2.)*I3b; */
/* 	/1* I[5] = alpha*I3ad*(I3a-1.)+(alpha/2.)*c_v*T_0*I3gd*(exp(S_()/c_v)-1.); *1/ */
/* 	I[5] = K_0*I3ad*(I3a-1.)/(2.*alpha) + (gamma/2.)*c_v*T_0*I3gd*(exp(S_()/c_v)-1.)+(B_0*beta/4.)*I3bd*(pow(I[0], 2)/3. - I[1]); */

	
/* } */

/* inline double ElasticState::e_() */
/* { */
/* 	/1* cout << "internal from Energy "<< E_() - mod_u()/2. << endl; *1/ */
/* 	return E_() - mod_u()/2.; */
/* } */

/* /1* // derivative incorrect *1/ */ 
/* /1* inline void ElasticState::stress(Eigen::Matrix3d G) *1/ */
/* /1* { *1/ */
/* /1* 	Eigen::Matrix3d GI = G.inverse(); *1/ */
/* /1* 	double gProd1, gProd2; *1/ */
/* /1* 	for(int i = 0; i < 3; i++) *1/ */	
/* /1* 		for(int j = 0; j < 3; j++) *1/ */
/* /1* 		{ *1/ */
/* /1* 			for(int k = 0; k < 3; k++) *1/ */
/* /1* 			{ *1/ */
/* /1* 				gProd1 = G(i,k)*G(j,k); *1/ */
/* /1* 				gProd2 = G(i,k)*GI(j,k); *1/ */
/* /1* 			} *1/ */
/* /1* 			sigma(i,j) = -2*rho*(G(i,j)*I[3] + I[4]*(I[0]*G(i,j)- *1/ */
/* /1* 							gProd1) + I[5]*I[2]*gProd2); *1/ */
/* /1* 		} *1/ */
/* /1* 	stressi(G); *1/ */
/* /1* 	cout << "first\n" << sigma << endl; *1/ */
/* /1* } *1/ */

/* inline void ElasticState::stress(Eigen::Matrix3d G) */
/* { */
/* 	/1* I[3] = 0; *1/ */
/* 	Eigen::Matrix3d GI = G.inverse(); */
/* 	double sigma_t; */
/* 	double term1 = 0., term2 = 0., term3 = 0.; */
/* 	for(int i = 0; i < 3; i++) */
/* 		for(int j = 0; j < 3; j++) */
/* 		{ */	
/* 			sigma_t = 0; */
/* 			for(int k = 0; k < 3; k++) */
/* 			{ */
/* 				term1 = I[3]*kron(j,k); */
/* 				term2 =	I[4]*(I[0]*kron(j,k) - G(j,k)); */
/* 				term3 = I[5]*(I[2]*GI(j,k)); */
/* 				sigma_t += -2.*rho*G(i,k)*(term1 + term2 + term3); */
/* 			} */	
/* 			sigma(i,j) = sigma_t; */	
/* 		} */
/* 	/1* cout << "second\n" << sigma << endl; *1/ */

/* 	/1* cout << G << '\n'<< endl; *1/ */
	
/* 	/1* cout << GI << '\n' << endl; *1/ */

/* 	/1* cout << I[0] << " " << I[1] << " " << I[2] << " " << I[3] << " " << I[4] << " " << I[5] << "\n " << endl; *1/ */
/* } */

/* inline double ElasticState::kron(int j, int k) */
/* { */
/* 	if (j == k) */
/* 		return 1.; */
/* 	else */
/* 		return 0; */
/* } */
/* //Will eventaully remove this function */
/* //Requires the first few invariants and entropy */

/* inline void ElasticState::prepareOutputs() */
/* { */
/* 	Eigen::Matrix3d rhoF; */
/* 	ElasticState temp = *this; */	
/* 	rhoF <<  v[3], v[4], v[5], */
/* 		     v[6], v[7], v[8], */
/* 		     v[9], v[10], v[11]; */
/* 	/1* rho = pow(rhoF.determinant()/2., 2.); *1/ */
/* 	rho = sqrt(rhoF.determinant()/rho_0); */
/* 	//deformation gradient */
/* 	Eigen::Matrix3d F = rhoF/rho; */
/* 	// Calculate finger tensor */
/* 	Eigen::Matrix3d G = strainTensor(F); */

/* 	//calculate invariants */
/* 	/1* cout << B_0 << endl; *1/ */
/* 	invariants(G); */
/* 	dInvariants_G(); */
	
/* 	//Stress */
/* 	stress(G); */
/* } */

/* double ElasticState::soundSpeed() */
/* { */
/* 	return sqrt(K_0); */
/* } */

/* void ElasticState::checkEOSConstants() */
/* { */
/* 	cout << "alpha " << alpha << "\nbeta " << beta <<"\ngamma "<< gamma << "\nb_0 " << b_0 << "\nc_0 " << c_0 << "\nc_v " << c_v << "\nK_0 " << K_0 << endl; */
/* } */

/* inline double ElasticState::rho_() */
/* { */
/* 	/1* rho = sqrt(rhoF.determinant()/rho_0); *1/ */
/* 	/1* cout << rho << endl; *1/ */
/* 	return rho; */
/* } */

/* // throw error if index is greater than certain value */
/* inline double ElasticState::u_(int index) */
/* { */
/* 	return v[index]/rho_(); */
/* } */

/* inline double ElasticState::sigma_(int index_i, int index_j) */
/* { */
/* 	return sigma(index_i, index_j); */
/* } */

/* inline double ElasticState::E_() */
/* { */
/* 	// This must be calculated after the flux */
/* 	return v[12]/rho_(); */
/* } */
/* //explicit expression for entropy calculated for stress */
/* inline double ElasticState::S_() */
/* { */	
/* 	// This can be calculated by rearranging the */ 
/* 	// expression for the internal energy as a fn of S */
/* 	// and working out what e is from e =|u|^2/2 - E */
/* 	// (since we know what E is from the solved conserved */
/* 	double I_3a = pow(I[2], alpha/2.); */
/* 	double I_3b = pow(I[2], beta/2.); */
/* 	double I_3g = pow(I[2], gamma/2.); */
/* 	double term1 = (K_0/pow(alpha,2)/2.)*pow((I_3a-1),2); */
/* 	double term2 = (B_0/2.)*I_3b*(pow(I[0], 2)/3. - I[1]); */
/* 	double term3 = c_v*T_0*I_3g; */
/* 	return c_v*log(1. + (e_() - term1 - term2)/term3); */
/* } */

/* //overload operators */
/* inline ElasticState ElasticState::operator+(ElasticState param) */
/* { */
/* 	ElasticState temp = *this; */	
/* 	for (int i = 0; i < e_size; i++) */
/* 	{ */
/* 		temp[i] += param[i]; */
/* 	} */	
/* 	return temp; */ 
/* } */

/* inline ElasticState ElasticState::operator-(ElasticState param) */
/* { */       
/*     ElasticState temp = *this; */					
/*     for (int i = 0; i < e_size; i++) */
/*     { */
/* 	    temp[i] -= param[i]; */     
/*     } */       
/*     return temp; */
/* } */

/* inline ElasticState ElasticState::operator*(double scalar) */
/* { */
/*     ElasticState temp = *this; */
/*     for (int i = 0; i < e_size; i++) */
/*     { */
/*         temp[i] *= scalar; */        
/*     } */                                    
/*     return temp; */
/* } */

/* inline ElasticState operator*(double k, ElasticState v) */
/* { */
/* 	return v*k; */
/* } */

/* inline ElasticState ElasticState::operator/(ElasticState rhs) */
/* { */
/* 	ElasticState temp = *this; */
/* 	for(int i = 0; i < e_size; i++) */
/* 	{ */
/* 		temp[i] /= rhs[i]; */
/* 	} */
/* 	return temp; */
/* } */	
	
/* /1* inline ElasticState ElasticState::operator=(ElasticState rhs) *1/ */
/* /1* { *1/ */
/* /1* 	/2* for(int i = 0; i < e_size; i++) *2/ *1/ */
/* /1* 	/2* { *2/ *1/ */
/* /1* 	/2* 	(*this)[i] = rhs[i]; *2/ *1/ */
/* /1* 	/2* } *2/ *1/ */
/* /1* 	// i may decide to change v back to a vector *1/ */ 
/* /1* 	// *1/ */
/* /1* 	(*this).v = rhs.v; *1/ */
/* /1* } *1/ */

/* inline std::ostream& operator<<(std::ostream& os, ElasticState& param) */
/* { */
/* 	/1* param.initialStates(os); *1/ */
/* 	/1* os << "Current states \n" << *1/ */
/* 	/1* 	  "Density\t"    << param.rho_() << "\n" *1/ */
/* 	/1* 	  "Energy\t"     << param.E_() << "\n" *1/ */
/* 	/1* 	  "U or F vector\t " << "\n"; *1/ */
/* 	os.precision(3); */
/* 	for (int i = 0; i < ElasticState::e_size; i++) */
/* 	{ */
/* 		os << setw(3) << " |" << setw(7) << param[i]; */
/* 	} */	
/* 	os << setw(7) << " | "; */
/* 	return os; */
/* } */

/* inline double ElasticState::getInvariant(int index) */
/* { */
/* 	return I[index]; */
/* } */
#endif
