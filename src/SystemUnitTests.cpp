#include "SolidSystem.h"
#include <vector> 
#include <iomanip>

//equation of state and system tests

void check_dI_dF(ElasticPrimState pri, System sys);
void check_depsi_dI_dI(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_depsi_dI_dS(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_deps_dF(ElasticPrimState pri, System sys, ElasticEOS eos);
void check_drho_dF(ElasticPrimState, System sys);
void check_dsigma_dG(ElasticPrimState pL, System sys, ElasticEOS eos, std::ofstream&);
void check_dsigma_deps(ElasticPrimState prim, System sys);
void check_B(ElasticPrimState prim, System sys, ElasticEOS eos);
void construct_Eigenvectors(System sys, ElasticPrimState pW, ElasticEOS eos);
double dot(ElasticPrimState state1, ElasticPrimState state2);
void callsToSystem(System sys, ElasticPrimState pL, ElasticPrimState pR);

typedef Eigen::Matrix<double, 13, 1> Vector13d;


void godunovState(System sys, ElasticPrimState pL, ElasticPrimState pR, ElasticEOS eos, std::vector<eigen>);
std::vector<eigen> construct_Eigenvectors_A(System sys, ElasticPrimState pL, ElasticPrimState pW, ElasticEOS eos);

using namespace Eigen;
using namespace std;

int main(void) 
{
  Eigen::Matrix3d FL, FR;
  Eigen::Vector3d uL, uR;
	double SL, SR;

  /* const ElasticPrimState priL(SquareMatrix(),DV3(5,500,1000),1000,0); */        

	/* FL <<	0.98,0.1,0.15,-0.02,1.1,0.1,0.15,0.015,1.2; */
	/* uL << 50., 500., 1000.; */
	/* SL = 1000; */

	FL <<	 0.98, 0, 0, 0.02, 1, 0.1, 0, 0, 1;
	uL << 0, 500, 1000;
	SL = 1000;

	uR << 0, 0, 0;
	FR << 1., 0, 0, 0, 1., 0.1, 0, 0, 1.;
	SR = 0.;

  ElasticPrimState primStateL(uL, FL, SL);  
  ElasticPrimState primStateR(uR, FR, SR);  
  ElasticEOS eos("copper"); // the romenski eos by default (to change)
  System sys(eos);
  ElasticState consStateL = sys.primitiveToConservative(primStateL);
  ElasticPrimState primStateCheck = sys.conservativeToPrimitive(consStateL);
  //std::cout << "Check that conversion works" << std::endl;
  assert(primStateL==primStateCheck);
  assert(primStateR!=primStateCheck);
  /* std::cout << "\nPRIMSTATE CHECK PASSED" << std::endl; */
  std::cout.precision(10); 
  ElasticPrimState pW = 0.5*(primStateL + primStateR);
  
	const Eigen::Matrix3d G = sys.strainTensor(pW.F_());
	const Eigen::Vector3d I = sys.getInvariants(pW.F_());
  int dirn = 0;
  SquareTensor3 dstressdF = sys.dstress_dF(pW, G, I, dirn);
  
  /* vector<eigen> egs = construct_Eigenvectors_A(sys, primStateL, primStateR, eos); */

  /* string file1 = "/home/raid/ma595/solid-1D/dstress_dFTest.out"; */
  /* ofstream dsdFout; */	
  /* dsdFout.open(file1.c_str(),ofstream::app); */


  /* dsdFout << "Matt's\ndstressdF\n " << dstressdF[0] << endl; */
  /* dsdFout << "\n " << dstressdF[1] << endl; */
  /* dsdFout << "\n " << dstressdF[2] << "\n" << endl; */

  /* check_dsigma_dG(pW, sys, eos, dsdFout); */
  /* dsdFout << "G\n" << G << endl; */
	/* const Eigen::Vector3d de_dI = eos.depsi_dI(I, pW.S_()); */
  /* dsdFout << "de_dI\n" << de_dI << endl; */
  /* dsdFout.close(); */
  //we can also do tests of the eigendecomposition

  //also experiment with single iteration, is this more efficient?
  callsToSystem(sys, primStateL, primStateR);
  construct_Eigenvectors(sys, pW, eos);
}


void callsToSystem(System sys, ElasticPrimState pL, ElasticPrimState pR)
{
  string file = "/home/raid/ma595/solid-1D/system.out";
  ofstream systemout;	
  systemout.open(file.c_str());
  ElasticPrimState god = sys.godunovState(pL, pR, 0);
  systemout << "Godunov State\n" << god << "\nLeft state\n" << pL << "\nRight state\n " << pR << endl;
  //output eigenvalues
}

//this is a test to sort objects  using the std::sort algorithm
//need to create objects of type eigen containing value and vector pairs
void godunovState(System sys, ElasticPrimState pL, ElasticPrimState pR, ElasticEOS eos, vector<eigen> eigs)
{
  /* vector<eigen> eigs = construct_Eigenvectors_A(sys, pL, pR, eos); */ 
  const ElasticPrimState dp = pR - pL;
  ElasticPrimState god = pL;
  for(int w=0; w<13; ++w) {
    if(eigs[w].lambda <0.0) {
      cout << eigs[w].lambda << endl;
      const double strength = dot(eigs[w].eigenvecL, dp);
      if(strength!=0.0) 
        god+=strength*eigs[w].eigenvecR;
    }
  }

  string file2 = "/home/raid/ma595/solid-1D/state_checkMatt.out";
  ofstream stateout;	
  stateout.open(file2.c_str());
  stateout << "pL\n" << pL << setprecision(5);
  stateout << "pR\n" << pR << setprecision(5);
  stateout << "Delta p\n" << dp << endl;
  stateout << "Godunov state\n" << god << setprecision(5);
  stateout.close();
   

  /* std::cout << "Godunov state \n" << god << std::endl; */
  /* std::cout << "pL state \n" << pL << std::endl; */
  /* std::cout << "pR state \n" << pR << std::endl; */

}

double dot(ElasticPrimState state1, ElasticPrimState state2)
{
  double temp = 0;
  for(int i = 0; i < 13; i++)
  {
    temp += state1[i]*state2[i]; 
  }
  return temp;
}

vector<eigen> construct_Eigenvectors_A(System sys, ElasticPrimState pL, ElasticPrimState pR, ElasticEOS eos)
{
  // construct curly A then decompose
  ElasticPrimState pW = 0.5*(pL + pR);
  int dirn = 0;
  const Eigen::Matrix3d F = pW.F_();
  const Eigen::Matrix3d FT = pW.F_().transpose();
	const Eigen::Matrix3d G = sys.strainTensor(F);
	const Eigen::Vector3d Inv = sys.getInvariants(F);
  SquareTensor3 dstressdF = sys.dstress_dF(pW, G, Inv, dirn);
  
  //construct curly A
  //is it possible to eigendecompose a typedef?
  //has to be an eigen type to perform the decomposition 
 
  /* typedef Eigen::Matrix<double, 13, 13> Matrix13d; //each row corresponds to eigenvector */
  
  Eigen::Matrix3d zero = Eigen::Matrix3d::Zero();
  Eigen::Vector3d zeroV = Eigen::Vector3d::Zero();
  
  double rho = sys.Density(pW);
  Eigen::Matrix3d A_11 = dstressdF[0];
  Eigen::Matrix3d A_12 = dstressdF[1];
  Eigen::Matrix3d A_13 = dstressdF[2];

  /* cout << "\ndstressdF\n" << dstressdF << endl; */
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  
  /* double u1 = pW.u_()(0); */
  double u1 = pW.u_()(0);

  const Eigen::Matrix3d m2rho = -2.*rho*G;
  const SquareTensor3 dsdeps = m2rho * pW.dI_dG(G,Inv); //cannot remember this derivation?
  Eigen::Matrix3d dsdeps_ = dsdeps[2];

  //might be worth creating this function in system
  Eigen::Matrix3d dsdS = eos.depsi_dI_dS(Inv, pW.S_())[2] * dsdeps_; //only non-zero component


  Eigen::Vector3d B;
  for(int i = 0; i < 3; i++) {
    B[i] = dsdS(0, i); 
  }
  B *= 1./rho;


  Eigen::Matrix3d E1 = I.col(dirn)*I.row(0);
  Eigen::Matrix3d E2 = I.col(dirn)*I.row(1);
  Eigen::Matrix3d E3 = I.col(dirn)*I.row(2);

  Eigen::MatrixXd AM(13,13);
  AM << u1*I, -(1./rho)*A_11, -(1./rho)*A_12, -(1./rho)*A_13, -B,
    -1.*FT*E1, u1*I, zero, zero, zero.col(0),
    -1.*FT*E2, zero, u1*I, zero, zero.col(0),
    -1.*FT*E3, zero, zero, u1*I, zero.col(0),
    zero.row(0), zero.row(0), zero.row(0), zero.row(0), u1;

  Eigen::EigenSolver<Eigen::MatrixXd> es(AM);
  /* cout << es.eigenvectors().col(0).norm(); */
  MatrixXd Re = es.eigenvectors().real();
  MatrixXd Le = es.eigenvectors().real().inverse();
  MatrixXd D = es.eigenvalues().real().asDiagonal();

  /* cout << "\nEigenvalues\n" << D << endl; */
  /* std::cout.precision(4); */ 
  vector<eigen> vecEig(13);
  /* typedef Eigen::Matrix<ElasticPrimState, 13, 1> Vector13Prim; //each row corresponds to eigenvector */

  for(int i = 0; i < 13; i++){
    Vector13d tempR;
    Vector13d tempL;
    for(int j = 0; j < 13; j++)
    {
      /* output << V(j,i) << "\t"; */
      tempR[j] = Re(j,i); 
      tempL[j] = Le(i,j); 
    }

    /* cout << temp.norm() << endl; */
    /* if(abs(D(i,i)) != 0) */
    /*   tempR /= D(i,i); */
    /* tempR /= 2.; */

    eigen etemp(D(i,i), tempL, tempR); 
    vecEig[i] = etemp;
  }

  VectorXd eigs = es.eigenvalues().real(); 
  std::sort(vecEig.begin(), vecEig.end(),less<eigen>()); //greater
  
  /* for(int i = 0; i < vecEig.size(); i++) */
  /* { */
  /*   std::cout << vecEig[i].lambda << endl; */
  /* } */


  string file = "/home/raid/ma595/solid-1D/eigenvectors.out";
  ofstream output;	
  output.open(file.c_str(), ofstream::app); //we append all results to the same file 
  output << "\nRight eigenvectors\n" << endl; 
  for(int i = 0; i < 13; i++)
  {
    ElasticPrimState eigvectemp = vecEig[i].eigenvecR;
    for(int j = 0; j < 13; j++)
    {
      output << eigvectemp[j] << "\t"; 
    }
    output << "\n";
  }

  //////////////////////////////////

  vector<double> sqsums(13);
  output << "\nLeft eigenvectors normalised\n" << endl; 
  for(int i = 0; i < 13; i++)
  {
    double sqsum = 0;
    ElasticPrimState eigvectemp = vecEig[i].eigenvecL;
    for(int j = 0; j < 13; j++)
    {
      output << eigvectemp[j] << "\t"; 
      sqsum += pow(eigvectemp[j], 2.);
    }
    sqsums[i] = pow(sqsum, 0.5);
    output << "\n";
  }
  output << "\n";

  for(int i = 0; i < 13; i++)
  {
    double sqsum = 0;
    ElasticPrimState eigvectemp = vecEig[i].eigenvecL;
    for(int j = 0; j < 13; j++)
    {
      output << eigvectemp[j]/sqsums[i] << "\t"; 
    }
    output << "\n";
  }


  for(int i = 0; i < 13; i++)
  {
    output.precision(7);
    output << vecEig[i].lambda << endl;
  }

  godunovState(sys, pL, pR, eos, vecEig);
  /* output << "AM\n " << AM << endl; */
  /* output << "B\n" << B << endl; */
  output.close();
  return vecEig;

  //these tests pass
  /* cout << Re*Le << endl; */
  /* cout << Le*Re << endl; */

  //as do these:
  /* for(int i = 0; i < 13; i++){ */
  /*   Vector13d vecEigL = vecEig[i].eigenvecL; */
  /*   Vector13d vecEigR = vecEig[i].eigenvecR; */
  /*   double val = 0; */
  /*   for(int j = 0; j < 13; j++) */
  /*   { */
  /*     val += vecEigL[j]*vecEigR[j]; */ 
  /*   } */
  /*   cout << val << endl; */
  /* } */

  /* cout << eos.depsi_dI_dS(Inv, pW.S_()) << endl; */
  /* cout << (1./rho) * dsdS << endl; */


  

  /* for(int i = 0; i < 13; i++) */
  /*   cout << "Length " << V.col(i).norm() << endl; */

  /* cout << "Reconstruct A from eigenvectors" << endl; */
  /* MatrixXd Arecon = V * D * V.inverse(); */
  /* cout << (Arecon - AM) << endl; */
  

  /* Eigen::MatrixXd eigs = es.eigenvectors(); */
  /* std::cout << "\n" << V << "\n" << std::endl; */
  /* std::cout << eigs << std::endl; */

  /* Eigen::MatrixXd A = MatrixXd::Random(6,6); */
  /* cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl; */
  /* EigenSolver<Eigen::MatrixXd> eis(A); */
  /* cout << "The eigenvalues of A are:" << endl << eis.eigenvalues() << endl; */
  /* cout << "The matrix of eigenvectors, V, is:" << endl << eis.eigenvectors() << endl << endl; */
  
}

void construct_Eigenvectors(System sys, ElasticPrimState pW, ElasticEOS eos)
{
  std::cout << "CONSTRUCT EIGENVECTORS\n" << std::endl;
  int dirn = 0;
  /* Eigen::MatrixXcd Q = es.eigenvectors(); */
  Eigen::Matrix3d omega = sys.AcousticTensor(pW, dirn);
  std::cout << "acoustic tensor\n" << omega << std::endl;
  Eigen::EigenSolver<Eigen::Matrix3d> es(omega);
  Eigen::Matrix3d V = es.pseudoEigenvectors();
  Eigen::Matrix3d eigs = es.pseudoEigenvalueMatrix();

  //we can check once they have been constructed
  //Q-1 = V
  //L = D^2
  Eigen::Matrix3d D;
  /* std::cout << L << std::endl; */
  D(0,0) = pow(eigs(0,0),0.5);
  D(1,1) = pow(eigs(1,1),0.5);
  D(2,2) = pow(eigs(2,2),0.5);
  Eigen::Matrix3d Dinv = D.inverse();
  /* std::cout << " D " << D << "\n" << std::endl; */
  Eigen::Matrix3d Q = V.inverse();
  Eigen::Matrix3d Qinv = V;
  /* std::cout << "Qinv\n " << Qinv << "\n " << std::endl; */
 
  Eigen::Matrix3d perm = Eigen::Matrix3d::Zero();
  //permutation matrix
  perm(2,0) = 1;
  perm(1,1) = 1;
  perm(0,2) = 1;

  const Eigen::Matrix3d F = pW.F_();
  const Eigen::Matrix3d Finv = F.inverse();
	const Eigen::Matrix3d G = sys.strainTensor(F);
	const Eigen::Vector3d I = sys.getInvariants(F);
  double rho = sys.Density(pW);
  const Eigen::Matrix3d m2rho = -2.*rho*G;
  
	const SquareTensor3 dsdeps = m2rho * pW.dI_dG(G,I); //cannot remember this derivation?
  Eigen::Matrix3d dsdeps_ = dsdeps[2];
  //might be worth creating this function in system
  Eigen::Matrix3d dsdS = eos.depsi_dI_dS(I, pW.S_())[2] * dsdeps_; //only non-zero component
  Eigen::Vector3d B;
  for(int i = 0; i < 3; i++) {
    B[i] = dsdS(dirn, i); 
  }
  B *= 1./rho;

  SquareTensor3 dstressdF = sys.dstress_dF(pW, G, I, dirn);
  Eigen::Matrix3d DQ = D*Q;
  Eigen::Matrix3d QA_11 = Q*(1./rho)*dstressdF[0];
  Eigen::Matrix3d QA_12 = Q*(1./rho)*dstressdF[1];
  Eigen::Matrix3d QA_13 = Q*(1./rho)*dstressdF[2];
  Eigen::Vector3d QB1 = Q*B;
  //the first row is B when dir = 0, self explanatory
 
  /* Eigen::Matrix3d mpiQA_11 = -perm*QA_11; */
  /* Eigen::Matrix3d mpiQA_12 = -perm*QA_12; */
  /* Eigen::Matrix3d mpiQA_13 = -perm*QA_13; */
  /* Eigen::Vector3d mpiQB1 = -perm*QB1; */

  /* ElasticPrimState primTest(pW.u_(), pW.F_(), pW.S_()); */
  /* std::cout << primTest << std::endl; */
  /* Eigen::Vector3d v(1,2,3); */
  /* primTest.u(v); */
  /* std::cout << primTest << std::endl; */
  /* Eigen::Matrix3d matTest = Eigen::Matrix3d::Ones(); */
  /* primTest.F(matTest); */
  /* std::cout << primTest << std::endl; */
  
  typedef Eigen::Matrix<ElasticPrimState, 13, 1> Vector13Prim; //each row corresponds to eigenvector
  
  Vector13Prim Le, Re;// initialise with zero components 
  ElasticPrimState la, law;
  Eigen::Matrix3d zero = Eigen::Matrix3d::Zero();
  Eigen::Vector3d zeroV(0,0,0);
  // it does make sense to wrap in elasticPrimState
  // a and aw = acoustic 
  // d and dp are linearly degenerate 

  int w, wr;
  for(w = 0; w < 3; w++)
  {
    wr = 12-w;

    Eigen::Matrix3d tl = zero; 
    tl << QA_11.row(w), QA_12.row(w), QA_13.row(w);

    la.u(DQ.row(w));
    la.F(tl);
    la.S(QB1(w));
    Le(w) = la;

    law.u(DQ.row(w));
    law.F(-tl);
    law.S(QB1(w));
    Le(wr) = law;   
  }

  w = 3;
  int pS = w;
  for(int i = 0; i < 3; i++)
  {
    ElasticPrimState ld, ldp;
    ld[pS] = F(0,1)/F(0,0);
    ld[pS+1] = -1;
    ldp[pS] = F(0,2)/F(0,0); 
    ldp[pS+2] = -1;
    pS+=3;
    Le(w) = ld;
    Le(w+1) = ldp;
    w+=2;
  }

  ElasticPrimState ld;
  ld.u(zeroV);
  ld.F(zero);
  ld.S(1.0);
  Le(9) = ld;
  //and linearly degenerate left eigenvectors
  //Right eigenvectors 
  //acoustic eigenvectors
  Eigen::Matrix3d QinvDinv = Qinv*Dinv;
  /* cout << "QINV*DINV " << QinvDinv <<  "\n " << Qinv << "\n " << Dinv << "\n" << endl; */
  Eigen::Matrix3d QinvDinvT = (Qinv*Dinv).transpose();
  Eigen::Matrix3d QinvDinv2 = QinvDinv*Dinv;
  /* Eigen::Matrix3d QinvDinv2 = (Qinv*Dinv*Dinv).transpose(); */
  Eigen::Matrix3d FT = F.transpose();
  
  for(w = 0; w < 3; w++)
  {
    wr = 12-w;
    ElasticPrimState ra, raw;
    //velocity components
    ra.u(QinvDinvT.row(w));
    raw.u(QinvDinvT.row(w));

    ra.F(QinvDinv2.col(w)*F.row(dirn));
    raw.F(-1.0*QinvDinv2.col(w)*F.row(dirn));

    //entropy part
    ra.S(0);
    raw.S(0);

    Re[w] = 0.5*ra;
    Re[wr] = 0.5*raw;
  }

  //indices 3 to 8
  //shear components
  Eigen::Matrix3d ominv = omega.inverse();
  //may be possible to vectorise this 
  Eigen::Matrix3d ominvA11 = ominv*(1./rho)*dstressdF[0];
  Eigen::Matrix3d ominvA12 = ominv*(1./rho)*dstressdF[1];
  Eigen::Matrix3d ominvA13 = ominv*(1./rho)*dstressdF[2];
  typedef Eigen::Matrix<Eigen::Matrix3d, 3, 1> Vector3M; //each row corresponds to eigenvector
  Vector3M ominvA(ominvA11, ominvA12, ominvA13);
  int s = 3;// these are the Re indices
  int sw = 4;

  for(int i = 0; i < 3; i++)
  {
    ElasticPrimState rs, rsw;
    rs.u(zeroV);
    rsw.u(zeroV);

    Eigen::Matrix3d Fs = ominvA[i].col(1)*F.row(dirn);
    Fs(i,1) -= 1.;
    rs.F(Fs);
    Eigen::Matrix3d Fsw = ominvA[i].col(2)*F.row(dirn);
    Fsw(i,2) -= 1.;
    rsw.F(Fsw);
    //need to subtract t
 
    rs.S(0);
    rsw.S(0);
    
    /* cout << s << endl; */
    /* cout << sw << endl; */
    
    Re[s] = rs;
    Re[sw] = rsw;

    s+=2;
    sw+=2;
  }

  //linearly degenerate eigenvector
  ElasticPrimState rd;
  Eigen::Vector3d ominvB = ominv * B;
  Eigen::Matrix3d FominvB = ominvB * F.row(dirn);
  rd.u(zeroV);
  rd.F(FominvB);
  rd.S(-1.0);
  Re[9] = -1.0*rd;

  std::cout << "Left eigenvectors " << std::endl;
  for(int i = 0; i < Le.size(); i++)
    std::cout << Le[i] << std::endl;

  std::cout << "Right eigenvectors " << std::endl;
  for(int i = 0; i < Re.size(); i++)
    std::cout << Re[i] << std::endl;

  /* std::cout << D << std::endl; */

  //check that LR = I 

  for(int i = 0; i < 13; i++)
  {
    std::cout << sys.dotState(Le[i],Re[i]) << endl;
  }




  //representation of unit dyads 
  //
  //
  //L = 13 * 13 

  //L = 13 * 13
  //This is quite cool: specify compile time sizes
  //do these need to be know at compile time?
  //we can make our own tinyvector this way
  /* Vector13d testtype = Vector13d::Zero(); */
  //may need to change this representation: think of Leveque's linearised Riemann solver
 
  /* typedef Eigen::Matrix<Vector13d, 13, 1> Vector13V; */
  /* Vector13V L, R; */
  /* L[0] = Vector13d::Ones(); */
  /* L[1] = Vector13d::Ones()*4; */
  /* std::cout << L[0].dot(L[1]) << std::endl; */
  /* //we can also get row vector and column vector using row and col functions */
  /* //we can also concatenate vectors using this approach */
  /* Eigen::Vector3d zeros = Eigen::Vector3d::Zero(); */
  /* Eigen::Vector3d testV1 = Eigen::Vector3d::Ones()*2; */
  /* Eigen::Vector3d testV2 = Eigen::Vector3d::Ones()*3; */
  /* Eigen::Matrix3d onesM = Eigen::Matrix3d::Ones(); */

  //could either use eigen or ElasticPrimState.
  //Eigen does give some flexibility

    /* L(i) << DQ.row(i), QA_11.row(i), QA_12.row(i), QA_13.row(i), 0; */
    /* L[12-i] << DQ.row(i), -QA_11.row(i), -QA_12.row(i), -QA_13.row(i), QB1(i); */
  /* } */


  /* L.transpose(); */
  /* std::cout << L[0] << std::endl; */

  // the r representation might not be correct here?
  
  /* std::cout << (Eigen::Matrix3d)perm << std::endl; */
  /* Eigen::Matrix3d permMatrix = (Eigen::Matrix3d)perm*Qinv; */
  /* std::cout << permMatrix << std::endl; */
  /* Matrix3d row1_3 = */ 
}
//problem with this approach is that I cannot assign by concatenating to a row or column of a matrix
//but it would be good to do it this way cause I could take the transpose at will.
  /* typedef Eigen::Matrix<double, 13, 13> Matrix13d; */ 
  /* Matrix13d L, R; */
  /* L = Matrix13d::Zero(); */
  /* std::cout << L.row(0) << std::endl; */
  /* typedef Eigen::Matrix<double, 13, 1> Vector13d; */
  /* Vector13d testtype = Vector13d::Zero(); */
  /* Vector13d::Ones(); */
  /* Eigen::Vector3d vec1; */
  /* vec1 = Eigen::Vector3d::Ones(); */

void compute_AcousticTensor(System sys, ElasticPrimState pW)
{
  int dirn = 0;
  Eigen::Matrix3d omega = sys.AcousticTensor(pW, dirn);
  std::cout << "Compute acoustic tensor test\n" << omega  << std::endl;
  printf("\nEigenvalues and orthogonal (eigenvectors) matrix Q\n"); 
  Eigen::EigenSolver<Eigen::Matrix3d> es(omega);
  std::cout << es.eigenvectors() << std::endl; 
  Eigen::MatrixXcd D = es.eigenvalues().asDiagonal();
  std::cout << D << std::endl;
  /* std::cout << "\n eigenvalues " << es.eigenvalues().asDiagonal() << std::endl; */ 
  /* std::cout <<  es.eigenvectors().col(0)[0] << std::endl; */ 
  printf("\nCheck that orthonormal\n");
  Eigen::MatrixXcd Q = es.eigenvectors();
  //how to only get non-trivial eigenvalues and eigenvectors from eigensolver
  Eigen::Matrix3d V = es.pseudoEigenvectors();
  Eigen::Matrix3d L = es.pseudoEigenvalueMatrix();

  printf("\nCheck that we can reconstruct the acoustic tensor: omega = Q*L*Q^-1");
  //this checks out
  Eigen::Matrix3d omegaR = V * L * V.inverse();
  std::cout << "\n" << omegaR  << std::endl;
}

void check_B(ElasticPrimState prim, System sys, ElasticEOS eos)
{
  std::cout << "Compute numerical result B " << std::endl;
  Eigen::Matrix3d sigma = sys.stress(prim);
  double h = 1e-7;
  ElasticPrimState prim_T(prim.u_(), prim.F_(), prim.S_() + h);
  Eigen::Matrix3d sigma_T = sys.stress(prim_T);

  std::cout << (sigma_T - sigma)/h << std::endl;
  printf("\n compute analytic result\n");

	const Eigen::Matrix3d G = sys.strainTensor(prim.F_());
	const Eigen::Vector3d I = sys.getInvariants(prim.F_());
  double rho = sys.Density(prim);
	const Eigen::Matrix3d m2rho = -2.*rho*G;
	const SquareTensor3 dsdeps = m2rho * prim.dI_dG(G,I); //cannot remember this derivation?
  /* std::cout << eos.depsi_dI_dS(I, prim.S_()) << std::endl; */
  check_depsi_dI_dS(prim, sys, eos);


  Eigen::Matrix3d dsdeps_ = dsdeps[2];
  Eigen::Matrix3d dsdS = eos.depsi_dI_dS(I, prim.S_())[2] * dsdeps_;
  //the first row is B when dir = 1, self explanatory
  std::cout << dsdS << std::endl;
  int dirn = 0;
  Eigen::Vector3d B;
  std::cout << " B " << std::endl;
  for(int i = 0; i < 3; i++)
  {
    B[i] = dsdS(dirn, i); 
  }
  B *= 1./rho;
  std::cout << B << std::endl;
}

//general matrix by matrix derivative
//pass in function pointers 
//try lambda functions C++11
void fd_dscalar_dmatrix(double scalar, Eigen::Matrix3d matrix)
{
  
}

//likely to be incorrect:
void check_dsigma_deps(ElasticPrimState prim, System sys)
{
  double rho = sys.Density(prim);
  Eigen::Matrix3d F = prim.F_();
	const Eigen::Vector3d I = sys.getInvariants(F);
	const Eigen::Matrix3d G = sys.strainTensor(prim.F_());
	const Eigen::Matrix3d m2rho = -2.*rho*G;
  const SquareTensor3 dsdeps = m2rho * prim.dI_dG(G,I);
  std::cout << "dsigma_deps analytical result: \n " <<  dsdeps << std::endl;
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
       
    }
}

void check_drho_dF(ElasticPrimState pL, System sys)
{
  double rhoL = sys.Density(pL);
  Eigen::Matrix3d FL = pL.F_();
  std::cout << "drho_dF analytical result: \n " << -rhoL*FL.transpose().inverse() << std::endl;
  Eigen::Vector3d uL = pL.u_();
  double SL = pL.S_();
  Eigen::Matrix3d FD_approx;
  double h = 1e-8; // for this test h = 1e-08 is the best, why is this? 
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
    {
      Eigen::Matrix3d FL_temp = FL;
      FL_temp(i,j) += h;
      ElasticPrimState primStateLTemp(uL, FL_temp, SL);
      double rhoph = sys.Density(primStateLTemp);
      FD_approx(i,j) = (rhoph-rhoL)/h;
    }
  std::cout << "drho_dF numerical approximation: \n" << FD_approx << std::endl;
}

//checked against kevin
void check_dsigma_dG(ElasticPrimState pL, System sys, ElasticEOS eos, std::ofstream& output)
{
  double rhoL = sys.Density(pL);
  Eigen::Matrix3d FL = pL.F_();
	const Eigen::Vector3d IL = sys.getInvariants(FL);
	const Eigen::Matrix3d GL = sys.strainTensor(pL.F_());
	const Eigen::Vector3d de_dI = eos.depsi_dI(IL, pL.S_());
  std::cout << "Analytical result: \n " << std::endl;

  output << "dsigma_dG\n" << endl;
  for(int i = 0; i < 3; i++)
  {
    const Eigen::Matrix3d ds_dG = pL.dsigma_dG(i, de_dI, GL, rhoL);
    output << ds_dG << "\n " << std::endl;
  }
}

void check_dG_dF(ElasticPrimState pL, System sys)
{
  std::cout << "Analytical: dG_dF " << std::endl;
  Eigen::Matrix3d FL = pL.F_();
	const Eigen::Matrix3d GL = sys.strainTensor(pL.F_());
  for (int j = 0; j < 3; j++)
    for (int m = 0; m < 3; m++){
      Eigen::Matrix3d dGdF = pL.dG_dF(GL,FL,j,m);
      int k1 = j*3 + m;
      std::cout << "\n" << k1 << ".\n" << dGdF << std::endl;
    }

}

//from Barton et al
//checked against kevin
void check_dI_dF(ElasticPrimState pri, System sys)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Matrix3d G = sys.strainTensor(F);
  const Eigen::Vector3d I = sys.getInvariants(F);
  std::cout << "Analytical: dI_dF\n" << pri.dI_dF(G,I) << std::endl;
}

void check_deps_dF(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Matrix3d GO = sys.strainTensor(F);
  const Eigen::Vector3d IO = sys.getInvariants(F);
	const std::vector<Eigen::Matrix3d> depsdF = sys.dep_dF(pri.dI_dF(GO,IO), eos.depsi_dI_dI(IO, pri.S_())); // this is the analytical solution
  std::cout << "Analytical: deps_dF\n" << depsdF << std::endl;

  std::cout << "Finite-difference approximation (correct): deps_dF\n" << std::endl;
  double h = 1e-7;
  ElasticPrimState priT = pri;
  for(int s = 0; s < 3; s++)
  {
    for(int i = 0; i < 3; i++){
      for(int j = 0; j < 3; j++) {
        Eigen::Matrix3d FTemp = pri.F_();
        FTemp(i,j) += h;
        const Eigen::Matrix3d G = sys.strainTensor(FTemp);
        const Eigen::Vector3d I = sys.getInvariants(FTemp);
        //wrap in primstate to calculate entropy
        Eigen::Vector3d Id = eos.depsi_dI(I, pri.S_());
        std::cout << (Id[s] - eos.depsi_dI(sys.getInvariants(F), pri.S_())[s])/h << "\t";
      }   
      std::cout << "\n";
    } 
    std::cout << "\n";
  }
}

//this checks out
void check_depsi_dI_dI(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Vector3d I = sys.getInvariants(F);
  Eigen::Matrix3d depsi_dI_dI = eos.depsi_dI_dI(I, pri.S_());
  printf("\ndepsi_dI_dI\n");
  std::cout << depsi_dI_dI << std::endl;
  //numerical bit
  double h = 1e-08;
  double S = pri.S_();
  Eigen::Vector3d eps_dI = eos.depsi_dI(I, S);
  //change each individual invariant:
  printf("\nFinite difference approximation: Correct \n");
  for(int i = 0; i < 3; i++)
  {
    Eigen::Vector3d IC_temp = I;
    IC_temp[i] += h;
    std::cout << (eos.depsi_dI(IC_temp, S) - eos.depsi_dI(I, S))/h << '\n' <<  std::endl;
  }

}

void check_depsi_dI_dS(ElasticPrimState pri, System sys, ElasticEOS eos)
{
  Eigen::Matrix3d F = pri.F_();
  const Eigen::Vector3d I = sys.getInvariants(F);
  Eigen::Vector3d depsi_dI_dS = eos.depsi_dI_dS(I, pri.S_());
  printf("\ndepsi_dI_dS\n");
  std::cout << depsi_dI_dS << std::endl;
  //numerical bit
  double h = 1e-08;
  double S = pri.S_();
  Eigen::Vector3d eps_dI = eos.depsi_dI(I, S);
  //change each individual invariant:
  printf("\nFinite difference approximation: Correct \n");
  std::cout << (eos.depsi_dI(I, S+=h) - eos.depsi_dI(I, S))/h << '\n' <<  std::endl;
}

