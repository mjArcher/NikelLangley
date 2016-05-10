#include "SolidSystem.h"

// getInvariants
using namespace Eigen;
using namespace std; 

System::System(const ElasticEOS eos):Eos(eos){}

System::~System(){}

/* #define debug_ */

//solid problem worth addressing for efficiency purposes is that when convert from 

ElasticPrimState System::conservativeToPrimitive(const ElasticState& consState) const 
{
	//calculate rho from inside class or calculate here
	const double rho = Density(consState);
	const Vector3d u = consState.mom()/rho;
	const Matrix3d F = consState.rhoF()/rho;
	double E = consState.rhoE()/rho;
	Vector3d inv = getInvariants(F);
	//! calculate internl energy from consState
	double ie =	E - u.dot(u)/2.; 	
	/* cout << "Cons INTERNAL ENERGY " << ie << endl; */
	double S = Eos.entropy(inv, ie);
	ElasticPrimState primState(u, F, S);
	return primState;
}

ElasticState System::primitiveToConservative(const ElasticPrimState& primState) const
{
	//finger tensor, invariants, internal energy, Energy (using eigen dot), 
	const Matrix3d F = primState.F_();
	const Vector3d invariants = getInvariants(F);
	/* cout << "Primitive invariants " << invariants << " G " << G << endl; */
	const double ie =	Eos.internalEnergy(invariants, primState.S_());
	Vector3d u = primState.u_();
	double E = ie + u.dot(u)/2.;
	/* cout << " Prim INTERNAL ENERGY " << ie << endl; */
	double rho = Density(primState);
	ElasticState consState(rho*primState.u_(), rho*F, rho*E);
	return consState;
}

double System::Density(const Eigen::Matrix3d& F) const 
{
	double rho = Eos.rho0/(F.determinant());	
	return rho;
}

double System::Density(const ElasticPrimState& primState) const 
{
	Matrix3d F = primState.F_();
	double rho = Eos.rho0/(F.determinant());	
	return rho;
}

double System::Density(const ElasticState& consState) const
{
	Matrix3d rhoF = consState.rhoF();	
	double rho = sqrt(rhoF.determinant()/Eos.rho0);	
	return rho;
}

Matrix3d System::strainTensor(const Matrix3d& F) const
{
	Matrix3d G = (F.inverse().transpose())*(F.inverse());	
	/* Matrix3d G = F.inverse()*(F.inverse().transpose()); */	
	return G;
}

Vector3d System::getInvariants(const Matrix3d& F) const
{
	const Matrix3d G = strainTensor(F);
	// Calculate invariants
	Vector3d I;
	I[0] = G.trace();
	I[1] = 0.5*(G.trace()*G.trace()-(G*G).trace());
	I[2] = G.determinant();
	return I;
}

ElasticState System::flux(const ElasticPrimState& primState) const 
{
  return flux(primitiveToConservative(primState), primState);
}

ElasticState System::flux(const ElasticState& consState) const
{
	return flux(consState, conservativeToPrimitive(consState));
}

ElasticState System::flux(const ElasticState& consState, const ElasticPrimState& primState) const
{
	ElasticState Fl; //flux
	Matrix3d rhoF = consState.rhoF();	
	double rhoE = consState.rhoE();
	Vector3d m = consState.mom();
	Vector3d u = primState.u_();
	Matrix3d sigma = stress(primState);
	//but we can use both the primitive state and conservative state to compute the flux
	// write for loops here instead

	Fl[0]  = m(0)*u(0) - sigma(0,0);
	Fl[1]  = m(1)*u(0) - sigma(0,1);
	Fl[2]  = m(2)*u(0) - sigma(0,2);
	Fl[3]  = 0; 
	Fl[4]  = 0; 
	Fl[5]  = 0; 
	Fl[6]  = rhoF(1,0)*u(0) - rhoF(0,0)*u(1);
	Fl[7]  = rhoF(1,1)*u(0) - rhoF(0,1)*u(1);
	Fl[8]  = rhoF(1,2)*u(0) - rhoF(0,2)*u(1);
	Fl[9]  = rhoF(2,0)*u(0) - rhoF(0,0)*u(2);
	Fl[10] = rhoF(2,1)*u(0) - rhoF(0,1)*u(2);
	Fl[11] = rhoF(2,2)*u(0) - rhoF(0,2)*u(2);
	Fl[12] = rhoE*u(0) - u(0)*sigma(0,0)-u(1)*sigma(0,1)-u(2)*sigma(0,2);
	return Fl;

	/* Fl[0]  = v[0]*u[0] - sigma(0,0); */
	/* Fl[1]  = v[1]*u[0] - sigma(0,1); */
	/* Fl[2]  = v[2]*u[0] - sigma(0,2); */
	/* Fl[3]  = 0; */ 
	/* Fl[4]  = 0; */ 
	/* Fl[5]  = 0; */ 
	/* Fl[6]  = rho*F(1,0)*u(0) - rho*F(0,0)*u(1); */
	/* Fl[7]  = rho*F(1,1)*u(0) - rho*F(0,1)*u(1); */
	/* Fl[8]  = rho*F(1,2)*u(0) - rho*F(0,2)*u(1); */
	/* Fl[9]  = rho*F(2,0)*u(0) - rho*F(0,0)*u(2); */
	/* Fl[10] = rho*F(2,1)*u(0) - rho*F(0,1)*u(2); */
	/* Fl[11] = rho*F(2,2)*u(0) - rho*F(0,2)*u(2); */
	/* Fl[12] = rho*u(0)*E_() - u(0)*sigma(0,0)-u(1)*sigma(0,1)-u(1)*sigma(0,2); */
}


//kronecker in utility functions
Matrix3d System::stress(const ElasticPrimState& primState) const
{
	Matrix3d F = primState.F_();
	Vector3d I = getInvariants(F);
	double rho = Density(primState);
	// can we just get the entropy from e primitive state?
	// yes but we still need to calculate the stress
	
	/* double S = Eos.entropy(inv, EosVgc */
	double S = primState.S_();
	Vector3d Id = Eos.depsi_dI(I, S);
	Matrix3d G = strainTensor(F);
	Matrix3d GI = G.inverse();
	double sigma_t;
	double term1 = 0., term2 = 0., term3 = 0.;
	Matrix3d sigma;
	//this needs to be reimplemented to take advantage of library matrix multiplication
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
		{	
			sigma_t = 0;
			for(int k = 0; k < 3; k++)
			{
				term1 = Id[0]*kron(j,k);
				term2 =	Id[1]*(I[0]*kron(j,k) - G(j,k));
				term3 = Id[2]*(I[2]*GI(j,k));
				sigma_t += -2.*rho*G(i,k)*(term1 + term2 + term3);
			}	
			sigma(i,j) = sigma_t;	
		}
	return sigma;
}

Matrix3d System::AcousticTensor(const ElasticPrimState& primState, const int dirn) const
{
  double rho = Density(primState);
	const Matrix3d F = primState.F_();
	const Matrix3d G = strainTensor(F);
	const Vector3d I = getInvariants(F);
  Matrix3d omega;
	const SquareTensor3 dsdFe = dstress_dF(primState, G, I,dirn); //need to decide on a return type here
	//construct acoustic tensor	i, j, k, m : i does not change for 1D 
  for(int i=0; i<3; ++i){
    for(int j=0; j<3; ++j){
      omega(i,j) = 0.0;
      for(int k=0; k<3; ++k){	
				omega(i,j) += dsdFe(i,j,k) * F(dirn,k); //symmetric 
      }
      omega(i,j) /= rho;
    }
  }
  return omega;
}

// reconsider the implementation of this
vector<Matrix3d> System::dep_dF(const SquareTensor3 dI_dF, const Matrix3d depsi_dI_dI) const
{
	vector<Matrix3d> depsi_dF(3);	

	Matrix3d dEps_pdF;

	for(int p = 0; p < 3; ++p)//row in depsi_dI_dI
	{
		dEps_pdF =  depsi_dI_dI(p, 0)*dI_dF[0] + depsi_dI_dI(p, 1)*dI_dF[1] + depsi_dI_dI(p, 2)*dI_dF[2]; //
		depsi_dF[p] = dEps_pdF;
	}

	/* depsi_dF[0] = depsi_dI_dI[0] * dI_dF[0]; //row */ 
	/* depsi_dF[1] = depsi_dI_dI[1] * dI_dF[1]; */
	/* depsi_dF[2] = depsi_dI_dI[2] * dI_dF[2]; */
	return depsi_dF;
}

double System::getMaxWaveSpeed(const ElasticPrimState& prim, const int dirn) const
{
  //can do this a number of ways
  //we can take the eigenvalues from decomposition of A (this means performing the expensive decomposition again)
  //or we could decompose the acoustic tensor - this would be be more efficient
  //or we could just guess
  double a = 10000;
  return abs(prim.u_()[dirn]) + a;
}

// Aijk 
// this is incorrect?
/* SquareTensor3 System::dstress_dF(const ElasticPrimState& primState, const Matrix3d& G, const Vector3d& I) const */
/* { */
/* 	const double rho = Density(primState); */
/* 	const Matrix3d F = primState.F_(); */	
/* 	const double ie =	Eos.internalEnergy(I, primState.S_()); */
/* 	const Matrix3d dstress_drho = stress(primState)/rho; */
/* 	const Matrix3d drho_dF = - rho*(F.inverse()).transpose(); */
/* 	const Matrix3d m2rho = -2.*rho*G; */
/* 	const vector<Eigen::Matrix3d> depsdF = dep_dF(primState.dI_dF(G,I), Eos.depsi_dI_dI(I, primState.S_())); */
/* 	const Vector3d de_dI = Eos.depsi_dI(I, primState.S_()); */
/* 	const SquareTensor3 dsdeps = m2rho * primState.dI_dG(G,I); //cannot remember this derivation? */
/* 	double sigma_rho; */

/*   int dirn = 0;  //function argument eventually */

/* 	vector<Matrix3d> A(3); */
/* 	const Matrix3d sigma = stress(primState); */

/* 	for (int k = 0; k < 3; k++) */
/*   { */
/* 		sigma_rho = dstress_drho(0,k); */
/* 		Matrix3d ds_dG = primState.dsigma_dG(k, de_dI, G, rho); */ 
/*     for(int j = 0; j < 3; j++) */
/*     { */
/*       for(int m = 0; m < 3; m++){ */
/* 				/1* A(k, m, j) -- row, column, depth -- index using A[j](k,m) - this is a std array *1/ */
/* 				Matrix3d dGdF = primState.dG_dF(G,F,j,m); //returns 2d slice of tensor at j, m (denominator constant) */
/* 				A[k](j,m) = sigma_rho * drho_dF(j,m) */ 
/* 					 + depsdF[0](j,m) * dsdeps[0](dirn,k) */ 
/*            + depsdF[1](j,m) * dsdeps[1](dirn,k) */ 
/*            + depsdF[2](j,m) * dsdeps[2](dirn,k) */
/*         + (dGdF * ds_dG).trace(); */
/* 			} */
/* 		} */
/* 	} */ 
/* 	SquareTensor3 dstressdF(A); */
/* 	return dstressdF; */
/* } */

SquareTensor3 System::dstress_dF(const ElasticPrimState& primState, const Matrix3d& G, const Vector3d& I, const int dirn) const
{
	const double rho = Density(primState);
	const Matrix3d F = primState.F_();	
	const double ie =	Eos.internalEnergy(I, primState.S_());
	const Matrix3d dstress_drho = stress(primState)/rho;
	const Matrix3d drho_dF = - rho * (F.inverse()).transpose();
	const Matrix3d m2rho = -2. * rho * G;
	const vector<Eigen::Matrix3d> depsdF = dep_dF(primState.dI_dF(G,I), Eos.depsi_dI_dI(I, primState.S_()));
	const Vector3d de_dI = Eos.depsi_dI(I, primState.S_());
	const SquareTensor3 dsdeps = m2rho * primState.dI_dG(G,I); //cannot remember this derivation?
	double sigma_rho;

	vector<Matrix3d> A(3);
	const Matrix3d sigma = stress(primState);

	for (int k = 0; k < 3; k++)
  {
		sigma_rho = dstress_drho(0,k);
		Matrix3d ds_dG = primState.dsigma_dG(k, de_dI, G, rho); 
    for(int j = 0; j < 3; j++)
    {
      for(int m = 0; m < 3; m++){
				/* A(k, m, j) -- row, column, depth -- index using A[j](k,m) - this is a std array */
				Matrix3d dGdF = primState.dG_dF(G,F,j,m); //returns 2d slice of tensor at j, m (denominator constant)
				A[j](k,m) = sigma_rho * drho_dF(j,m) 
					 + depsdF[0](j,m) * dsdeps[0](0,k) 
           + depsdF[1](j,m) * dsdeps[1](0,k) 
           + depsdF[2](j,m) * dsdeps[2](0,k);
        /* + (dGdF * ds_dG).trace(); */
        //removed the trace 
        double sum = 0;
        for(int q = 0; q < 3; q++){
          for(int r = 0; r < 3; r++){
            sum += ds_dG(q,r) * dGdF(q,r);
          }
        }
        A[j](k,m) += sum;
			}
		}
	} 

	SquareTensor3 dstressdF(A);
	return dstressdF;

}

Vector3d System::B(const ElasticPrimState& pW) const
{
  const Matrix3d F = pW.F_();
	const Matrix3d G = strainTensor(F);
	const double rho = Density(pW);
  const Matrix3d m2rho = -2.*rho*G;
	const Vector3d Inv = getInvariants(F);
  const SquareTensor3 dsdeps = m2rho * pW.dI_dG(G,Inv); //cannot remember this derivation?
  Matrix3d dsdeps_ = dsdeps[2];
  Matrix3d dsdS = Eos.depsi_dI_dS(Inv, pW.S_())[2] * dsdeps_; //only non-zero component
  Vector3d B;
  for(int i = 0; i < 3; i++) {
    B[i] = dsdS(0, i); 
  }
  B *= 1./rho;

  return B;
  /* cout << "B\n" << endl; */
  /* cout << B << endl; */
}


/* struct eigen { */ 
/*   double lambda; */
/*   ElasticPrimState eigenvecR; */
/*   ElasticPrimState eigenvecL; */
  
/*   //default constructor */
/*   eigen(){} */

/*   eigen(const double val, Vector13d vecL, Vector13d vecR){ */
/*     lambda = val; */
/*     eigenvecL = convert(vecL); */
/*     eigenvecR = convert(vecR); //normalised */ 
/*   } */

/*   ElasticPrimState convert(Vector13d vecR) */
/*   { */
/*     ElasticPrimState temp; */
/*     for(int i = 0; i < 13; i++) */
/*     { */
/*       temp[i] = vecR[i]; */
/*     } */
/*     return temp; */
/*   } */

/*   //descending order */ 
/*   bool operator>(eigen const &other) const { */ 
/*     return lambda > other.lambda; */
/*   } */

/*   //ascending order */ 
/*   bool operator<(eigen const &other) const { */ 
/*     return lambda < other.lambda; */
/*   } */
/* }; */

/* vector<eigen> construct_Eigenvectors_A(System sys, ElasticPrimState pL, ElasticPrimState pR, ElasticEOS eos) */
/* { */


vector<eigen> System::stateEigenDecompose_A(const ElasticPrimState& pL, const ElasticPrimState& pR, const int dirn) const
{
  ElasticPrimState pW = 0.5*(pL + pR);
	const double rho = Density(pW);
  const Eigen::Matrix3d F = pW.F_();
  const Eigen::Matrix3d FT = pW.F_().transpose();
	const Matrix3d G = strainTensor(F);
	const Eigen::Vector3d Inv = getInvariants(F);
  SquareTensor3 dstressdF = dstress_dF(pW, G, Inv, dirn);
	const ElasticState consState = primitiveToConservative(pW);
  //construct curly A
  /* typedef Eigen::Matrix<double, 13, 13> Matrix13d; //each row corresponds to eigenvector */
  Eigen::Matrix3d zero = Eigen::Matrix3d::Zero();
  Eigen::Vector3d zeroV = Eigen::Vector3d::Zero();
  
  Eigen::Matrix3d A_11 = dstressdF[0];
  Eigen::Matrix3d A_12 = dstressdF[1];
  Eigen::Matrix3d A_13 = dstressdF[2];

#ifdef debug
  cout << "\ndstressdF\n" << dstressdF[0] << "\n" << dstressdF[1] << "\n" << dstressdF[2] << "\n" << endl;
#endif
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  double u1 = pW.u_()(dirn);
  Eigen::Matrix3d E1 = I.col(dirn)*I.row(0);
  Eigen::Matrix3d E2 = I.col(dirn)*I.row(1);
  Eigen::Matrix3d E3 = I.col(dirn)*I.row(2);
  Eigen::Vector3d B_ = B(pW);
  Eigen::MatrixXd AM(13,13);

  AM << u1*I, (-1./rho)*A_11, (-1./rho)*A_12, (-1./rho)*A_13, -1.0*B_,
     -1.*FT*E1, u1*I, zero, zero, zero.col(0),
     -1.*FT*E2, zero, u1*I, zero, zero.col(0),
     -1.*FT*E3, zero, zero, u1*I, zero.col(0),
    zero.row(0), zero.row(0), zero.row(0), zero.row(0), u1;

  Eigen::EigenSolver<Eigen::MatrixXd> es(AM);
  MatrixXd Re = es.eigenvectors().real();
  MatrixXd Le = es.eigenvectors().real().inverse(); //more efficient way of doing this 
  VectorXd lambda = es.eigenvalues().real();
  //check if eigenvalues are real
  /* if(eval[0]<=0.0 || eval[1]<=0.0 || eval[2]<=0.0) { */
  /*   cerr << "Error: non-positive eigenvalue in InterfaceDecompose\n"; */
  /*   exit(1); */
  /* } */
  vector<eigen> vecEig(13);

  for(int i = 0; i < 13; i++){
    //typedefs
    Vector13d tempR;
    Vector13d tempL;
    for(int j = 0; j < 13; j++)
    {
      tempR[j] = Re(j,i); 
      tempL[j] = Le(i,j); 
    }
    eigen etemp(lambda(i), tempL, tempR); 
    vecEig[i] = etemp;
  }
  std::sort(vecEig.begin(), vecEig.end(),less<eigen>()); //greater
  return vecEig;
}

VectorXd System::stateEigenDecompose_A_old(const ElasticPrimState& pW, 
		const int dirn,
		vector<ElasticPrimState>& LeE,
		vector<ElasticPrimState>& ReE) const 
{
	const double rho = Density(pW);
  const Eigen::Matrix3d F = pW.F_();
  const Eigen::Matrix3d FT = pW.F_().transpose();
	const Matrix3d G = strainTensor(F);
	const Eigen::Vector3d Inv = getInvariants(F);
  SquareTensor3 dstressdF = dstress_dF(pW, G, Inv, dirn);
	const ElasticState consState = primitiveToConservative(pW);
  //construct curly A
  /* typedef Eigen::Matrix<double, 13, 13> Matrix13d; //each row corresponds to eigenvector */
  Eigen::Matrix3d zero = Eigen::Matrix3d::Zero();
  Eigen::Vector3d zeroV = Eigen::Vector3d::Zero();
  
  Eigen::Matrix3d A_11 = dstressdF[0];
  Eigen::Matrix3d A_12 = dstressdF[1];
  Eigen::Matrix3d A_13 = dstressdF[2];

  
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  double u1 = pW.u_()(0);
  //might be worth creating this function in system

  Eigen::Matrix3d E1 = I.col(dirn)*I.row(0);
  Eigen::Matrix3d E2 = I.col(dirn)*I.row(1);
  Eigen::Matrix3d E3 = I.col(dirn)*I.row(2);

  Eigen::Vector3d B_ = B(pW);

  Eigen::MatrixXd AM(13,13);
  AM << u1*I, (-1./rho)*A_11, (-1./rho)*A_12, (-1./rho)*A_13, -1.0*B_,
     -1.*FT*E1, u1*I, zero, zero, zero.col(0),
     -1.*FT*E2, zero, u1*I, zero, zero.col(0),
     -1.*FT*E3, zero, zero, u1*I, zero.col(0),
    zero.row(0), zero.row(0), zero.row(0), zero.row(0), u1;

  Eigen::EigenSolver<Eigen::MatrixXd> es(AM);
  MatrixXd Re = es.eigenvectors().real();
  MatrixXd Le = es.eigenvectors().real().inverse(); //more efficient way of doing this 
  VectorXd lambda = es.eigenvalues().real();
  vector<eigen> vecEig(13);

  for(int i = 0; i < 13; i++){
    //typedefs
    Vector13d tempR;
    Vector13d tempL;
    for(int j = 0; j < 13; j++)
    {
      tempR[j] = Re(j,i); 
      tempL[j] = Le(i,j); 
    }
    eigen etemp(lambda(i), tempL, tempR); 
    vecEig[i] = etemp;
  }

  VectorXd eigs = es.eigenvalues().real(); 
  std::sort(vecEig.begin(), vecEig.end(),less<eigen>()); //greater
  
}

//construct and then decompose acoustic tensor: get Q, D, A and then construct Le, Re and 
VectorXd System::stateEigenDecompose_Omega(const ElasticPrimState& pW,
		const int dirn,
		vector<ElasticPrimState>& Le,
		vector<ElasticPrimState>& Re) const 
{
	const double rho = Density(pW);
	const Matrix3d G = strainTensor(pW.F_());
	/* const Invariants inv = G.getInvariants(); */
	const ElasticState consState = primitiveToConservative(pW);
	
	const Matrix3d omega = AcousticTensor(pW, dirn);

  VectorXd temp(3);
  
  return temp;
}

ElasticState System::godunovFlux(const ElasticState& qL, const ElasticState& qR) const
{
  return flux(godunovState(conservativeToPrimitive(qL), conservativeToPrimitive(qR), 0)); //add dirn to this 
}

ElasticPrimState System::godunovState(const ElasticPrimState& pL, const ElasticPrimState& pR, const int dirn) const
{
  if (pL == pR)
    return pL;
  vector<eigen> eigs = stateEigenDecompose_A(pL, pR, dirn);
  ElasticPrimState dp = pR - pL;
  ElasticPrimState god = pL;
  for(int w=0; w<13; ++w) {
    if(eigs[w].lambda <0.0) {
      /* cout << eigs[w].lambda << endl; */
      const double strength = dotState(eigs[w].eigenvecL, dp);
      if(strength!=0.0) 
        god+=strength*eigs[w].eigenvecR;
    }
  }

#ifdef debug_
  cout << "\npL\n" << pL << "\npR\n" << pR << "\nGodunov state\n" << god << "\nEigenvalues\n " << endl;  
  for(int w = 0; w < 13; ++w)
  {
    if(eigs[w].lambda <0.0) 
      cout << eigs[w].lambda << endl;
  }
#endif

  return god;
}

double System::dotState(ElasticPrimState& state1, ElasticPrimState& state2) const
{
  double temp = 0;
  for(int i = 0; i < 13; i++)
  {
    temp += state1[i]*state2[i]; 
  }
  return temp;
}

/* ElasticPrimState System::godunovState(ElasticPrimState pL, ElasticPrimState pR) const */
/* { */
/*   if (pL == pR) */
/*     return pL; */
/*   // compute the eigenvalues, L and R eigenvectors. */
/*   const ElasticPrimState pW = 0.5*(pL + pR); */
/*   vector<ElasticPrimState> Le(14), Re(14); //test if this has actually initialised correctly. */
  /* const VectorXd lambda = stateEigenDecompose(pW, 0, Le, Re); // within this construct the acoustic tensor */

  /* /1* VectorXd VectorXd(ElasticPrimState::e_size); *1/ */
  /* //compute the primitive variable jump across the Riemann problem: */
  /* const ElasticPrimState dp = pR - pL; */ 

  /* ElasticPrimState god = pL; */
  /* for(int w=0; w<14; ++w) { */
  /*   if(lambda[w]<0.0) { */
  /*     //do dot manually */

  /*     /1* double strength = 0; *1/ */
  /*     /1* for(int i = 0; i < 13; i++) *1/ */ 
  /*     /1* { *1/ */
  /*     /1*   strength += Le[w](i)(dp[i]); *1/ */
  /*     /1* } *1/ */
  /*     /1* /2* dot(Le[w], dp); //need to rewrite this line *2/ *1/ */ 
  /*     /1* if(strength!=0.0) *1/ */
  /*     /1*   god+=strength*Re[w]; *1/ */
  /*   } */
  /* } */
  /* return god; */
/* } */
	
/* TinyVector<double, 14> System::StateEigenDecompose(const ElasticPrimState& pW, */
/* 		const int dirn, */
/* 		TinyVector<ElasticPrimState, 14>& Le, */
/* 		TinyVector<ElasticPrimState, 14>& Re) const */ 
/* { */
//in equation of state
//deps_dI
//deps_dI_dI

//Primstate
/* dI_dF */
//calculate in dstress_dF function
/* depsi_dF */

// Aij = dsigma_dF


