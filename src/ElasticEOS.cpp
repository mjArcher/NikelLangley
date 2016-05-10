#include "ElasticEOS.h"

using namespace std;
using namespace Eigen;

// internal energy, entropy and various derivatives - eventually extend to a inheriting class
//
//constructor
ElasticEOS::ElasticEOS(const string& mat)
{
	double c0, b0;

	if(mat == "copper")
	{
		rho0 = 8930.;
		c0 = 4600.;
		b0 = 2100.;
		cv = 390.;
		T0 = 300.;
		alpha = 1.001;
		beta = 3.001;
		gamma = 2.001;
	}
	else if(mat == "aluminium")
	{ //need to be changed
    rho0 = 2710; // kg m^-3
    c0   = 6220; // m s^-1
    cv   = 900; // J kg^-1 K^-1
    T0   = 300; // K
    b0   = 3160; // m s^-1
    alpha= 1.001;
    beta = 3.577;
    gamma = 2.088;
		//etc
	} //need to be changed
	else if(mat == "steel")
	{
    rho0 = 8030; // kg m^-3
    c0   = 5680; // m s^-1
    cv   = 500; // J kg^-1 K^-1
    T0   = 300; // K
    b0   = 3100; // m s^-1
    alpha= 0.596;
    beta = 2.437;
    gamma = 1.563;
	}
	else
  {
    cerr << "ERROR: Unknown material " << mat << " in Romenski elastic eos.\n";
    cerr << "       Options are: aluminium, steel, copper\n";
    exit(1);
  }

	B0 = b0*b0;
	K0 = c0*c0-(4./3.)*B0;

}
//! Destructor
ElasticEOS::~ElasticEOS(){}

void ElasticEOS::checkEosConstants() const
{
	cout << "alpha " << alpha << "\nbeta " << beta <<"\ngamma "<< gamma << "\nB0 " << B0 <<  "\ncv " << cv << "\nK0 " << K0 << endl;

}

//calculate invariantes in system class
//
//
double ElasticEOS::internalEnergy(Vector3d I, double S) const
{
	double U, W;
	/* cout << I[2] << " cv "; */
	/* cout << "Here 3 " << endl; */

	U = (K0/(2.*alpha*alpha))*pow((pow(I(2),alpha/2.))-1.,2.)+ cv*T0*pow(I(2),gamma/2.)*(exp(S/cv)-1.);
	W = (B0/2.)*pow(I(2),beta/2.)*(pow(I(0),2)/3. - I(1));

	/* cout << S << " " << alpha << " " << beta << " " << gamma << " " << cv << " " << B0 << " " << T0 << " " << K0 << endl; */
	/* cout << I << endl; */
	/* U = (K0/(2.*alpha*alpha))*pow((pow(I(2),alpha/2.))-1,2.)+ */
	/* 	cv*T0*pow(I(2),gamma/2.)*(exp(S/cv)-1.); */
	/* cout << " U " << U << endl; */
	/* W = (B0/2.)*pow(I(2),beta/2.)*(pow(I(0),2)/3. - I(1)); */
	/* cout << "W " << W << endl; */

	return U + W;	
}

double ElasticEOS::entropy(Vector3d I, double internalEnergy) const
{	
	// This can be calculated by rearranging the 
	// expression for the internal energy as a fn of S
	// and working out what e is from e =|u|^2/2 - E
	// (since we know what E is from the solved conserved
	double I_3a = pow(I[2], alpha/2.);
	double I_3b = pow(I[2], beta/2.);
	double I_3g = pow(I[2], gamma/2.);

	/* double term1 = 0.5*(K0/pow(alpha,2))*pow((I_3a-1.),2.); */
  double coeff = 0.5;
	double term1 = coeff*(K0/pow(alpha,2.))*(I_3a-1.)*(I_3a-1.);
	double term2 = (B0/2.)*I_3b*(pow(I[0], 2.)/3. - I[1]);
	double term3 = cv*T0*I_3g;

  /* std::cout << cv << " " << term1 << " " << term2 << " " << term3 << " " << internalEnergy << std::endl; */
	return cv*log(1. + (internalEnergy - term1 - term2)/term3);

}

//rethink where this function should go
//retun vector of derivatives
Vector3d ElasticEOS::depsi_dI(const Vector3d I, double S) const
{
	/* double alpha = 1.0, b_0 = 2.1, beta = 3.0; */
	/* double B_0 = b_0*b_0; */ 
	/* //dElasticState/DI_3, dW/dI_1, dW/dI_2, dW/dI_3 */
	Vector3d Id;

	double I3a = pow(I[2], alpha/2.);
	double I3b = pow(I[2], beta/2.);
	double I3ad = pow(I[2], alpha/2. - 1.);
	double I3bd = pow(I[2], beta/2. - 1.);
	double I3gd = pow(I[2], gamma/2. - 1.);

	Id[0] = (B0/3.)*I3b*I[0];
	Id[1] = (-B0/2.)*I3b;
	/* I[5] = alpha*I3ad*(I3a-1.)+(alpha/2.)*c_v*T_0*I3gd*(exp(S_()/c_v)-1.); */
	Id[2] = K0*I3ad*(I3a-1.)/(2.*alpha) + (gamma/2.)*cv*T0*I3gd*(exp(S/cv)-1.)+(B0*beta/4.)*I3bd*(pow(I[0], 2)/3. - I[1]);

	return Id;
}

//see notes on how this is arranged
Matrix3d ElasticEOS::depsi_dI_dI(const Vector3d& I, double S) const
{
	Matrix3d depsdI2;

	const double I3a = pow(I[2], alpha/2.);
	const double I3b = pow(I[2], beta/2.);
	const double I3ad = pow(I[2], alpha/2. - 1.);
	const double I3bd = pow(I[2], beta/2. - 1.);
	const double I3gd = pow(I[2], gamma/2. - 1.);

	const double I3aa2 = pow(I[2], alpha/2. - 2.);
	const double I3ab2 = pow(I[2], beta/2. - 2.);
	const double I3ag2 = pow(I[2], gamma/2. - 2.);
	//first row 
	depsdI2(0,0) = (1./3.)*B0*I3b; //dI1dI1
	depsdI2(0,1) = 0;  //dI1dI2
	depsdI2(0,2) = (beta/6.)*B0*I3bd*I[0]; //dI1dI3 
	//second row
	depsdI2(1,0) = 0; //dI2dI1
	depsdI2(1,1) = 0; //dI2dI2
	depsdI2(1,2) = -(1./4.)*B0*beta*I3bd; //dI2dI3 
	//third row
	double fact = B0*beta*I3bd;
	depsdI2(2,0) = (1./6.)*fact*I[0]; //dI3dI1
	depsdI2(2,1) = -(1./4.)*fact; //dI3dI2
  depsdI2(2,2) = 1./4.*B0*(beta/2.-1.)*beta*(I[0]*I[0]/3.-I[1])*I3ab2
                    + 1./4.*pow(I[2],alpha-2.)*K0+(alpha/2.-1.)*I3aa2*(I3a-1.)*K0/(2.*alpha)
                      + 1./2.*cv*T0*(exp(S/cv)-1.)*(gamma/2.-1.)*gamma*I3ag2;
	return depsdI2;
}

//only one non-zero derivative deps/dI3/dS
Vector3d ElasticEOS::depsi_dI_dS(const Vector3d& I, double S) const
{
  //differentiated with mathematica
  Vector3d depsi_dS;
  depsi_dS[0] = 0;
  depsi_dS[1] = 0;
  depsi_dS[2] =  0.5*T0*exp(S/cv)*gamma*pow(I[2],gamma/2. - 1);
  return depsi_dS;
}

double ElasticEOS::soundSpeed() const
{
	return sqrt(K0);
}


/* // derivative incorrect */ 
/*  void ConsState::stress(Matrix3d G) */
/* { */
/* 	Matrix3d GI = G.inverse(); */
/* 	double gProd1, gProd2; */
/* 	for(int i = 0; i < 3; i++) */	
/* 		for(int j = 0; j < 3; j++) */
/* 		{ */
/* 			for(int k = 0; k < 3; k++) */
/* 			{ */
/* 				gProd1 = G(i,k)*G(j,k); */
/* 				gProd2 = G(i,k)*GI(j,k); */
/* 			} */
/* 			sigma(i,j) = -2*rho*(G(i,j)*I[3] + I[4]*(I[0]*G(i,j)- */
/* 							gProd1) + I[5]*I[2]*gProd2); */
/* 		} */
/* 	stressi(G); */
/* 	cout << "first\n" << sigma << endl; */
/* } */
