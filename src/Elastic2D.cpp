#include <libconfig.h++>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <fenv.h>
#include <ctime>

#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "SolidSystem.h"
#include "InputSolid.h"

using namespace std;
using namespace libconfig;
using namespace Eigen;

/* #define debug_ */
/* #define debugrun_ */

//use kevin's domain and Material structs as a base to modify later
//work by solving over each material in turn
//(How does 1D ghost fluid work?)
//two or more separate domains which we solve then repopulate 'solution vector'
//system, dom and EOS associated with each material

enum slope_lim{superbee, minbee, vanleer, albalda};//0,1,2,3
slope_lim sl;

struct Domain {
	int Ni; //number of x cells
	int GNi; // number of x cells plus ghost cells
	int GC; // number of ghost cells
	int starti; //start of domain
	int endi; //end of domain
	double Lx; // length of domain
	double dx; // cell size

	Domain(const int a_Ni, const int a_GC, const double a_Lx) :
		Ni(a_Ni), GNi(a_Ni+2*a_GC), GC(a_GC),	starti(a_GC), endi(a_Ni+a_GC),
		Lx(a_Lx), dx(a_Lx/a_Ni){}
  Domain(){}
	~Domain(){
	}
};

struct Material {
  const System sys;
  vector<ElasticState> sol;  
  const Domain dom;
  string name;

  Material(const string& a_name, const Domain& a_dom, const ElasticEOS Eos) :
    sys(Eos), 
    sol(a_dom.GNi), dom(a_dom), name(a_name) {}
};

double slopelim(double);
double ksi_r(double);
ElasticState grad(const Material&, int);
void outputAll(string, const vector<ElasticState>);

//ppm functions
ElasticState theta(const Material& mat, int j);

void BCs(Material& mat) {  

//transmissive boundary conditions 
//start
#pragma omp parallel for schedule(dynamic)
  for(int i=0; i<mat.dom.starti; ++i) {
    const int imagei = 2*mat.dom.starti-1-i;
    mat.sol[i] = mat.sol[imagei];
  }

//end
#pragma omp parallel for schedule(dynamic)
  for(int i=mat.dom.endi; i<mat.dom.GNi; ++i) {
    const int imagei = 2*mat.dom.endi-1-i;
    mat.sol[i] = mat.sol[imagei];
  }  
}

void ICInterface(Material& mat,
    const double iface,
    const ElasticPrimState& left, const ElasticPrimState& right) 
{
  //convert to conserved states
  ElasticState stateL = mat.sys.primitiveToConservative(left);
  ElasticState stateR = mat.sys.primitiveToConservative(right);

  for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
    const double x = (i-mat.dom.starti+0.5)*mat.dom.dx;
    if(x<iface) 
    {      
      mat.sol[i] = stateL;
    }
    else 
    {
      mat.sol[i] = stateR;
    } 
  }
}

//supply left state 
//this will affect the calculation of the timestep
void solveXGodunov(Material& mat, const double dt)
{
	const double dt_dX = dt/mat.dom.dx;
  const vector<ElasticState> soln = mat.sol;
  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; i++)
	{
	  mat.sol[i] += dt_dX*(mat.sys.godunovFlux(soln[i-1], soln[i]) - mat.sys.godunovFlux(soln[i], soln[i+1]));
	}
  /* exit(1); */
	//printArray(U);
	BCs(mat);
}

ElasticState forceFlux(System& sys, vector<ElasticState>& left, vector<ElasticState>& right, double dt_dX, int i)
{
	const ElasticState F_lf = 0.5*(sys.flux(right[i]) + sys.flux(left[i+1])) + (0.5/dt_dX)*(right[i]-left[i+1]);
	const ElasticState C_r = 0.5*(right[i] + left[i+1]) + 0.5*(dt_dX)*(sys.flux(right[i]) - sys.flux(left[i+1]));
	ElasticState F_force = 0.5*(F_lf + sys.flux(C_r));
	return F_force;
}

double alpha_r(double d_r, double beta_r)
{
  double denom = pow(1e-6+beta_r, 2);
  return d_r/denom;
}

vector<ElasticState> beta(vector<ElasticState>& soln, int k, int i)
{
  //j should never change here
  // repeatedly calculating at least the denominator?
  // think this could be vectorised and stored 
  vector<ElasticState> b(k);
  for(unsigned int j = 0; j < ElasticState::e_size; j++)
  {
    b[0][j] = 13./12.*pow((soln[i][j] -   2.*soln[i+1][j] + soln[i+2][j]), 2.) + 1./4.*pow((3.*soln[i][j] - 4.*soln[i+1][j] + soln[i+2][j]), 2);
    b[1][j] = 13./12.*pow((soln[i-1][j] - 2.*soln[i][j]   + soln[i+1][j]), 2.) + 1./4.*pow(( soln[i-1][j] +    soln[i+1][j]), 2.);
    b[2][j] = 13./12.*pow((soln[i-2][j] - 2.*soln[i-1][j] + soln[i][j]), 2.)   + 1./4.*pow(( soln[i-2][j] - 4.*soln[i-1][j] + 3.*soln[i][j]), 2);
  }
  return b;
}

vector<ElasticState> omega(Material& mat, vector<ElasticState> vr, int soln_ind, bool right)
{
  const double d[3] = {3./10., 3./5., 1./10.}; // and dbar_r = d_{k-1-r}
  int k = vr.size();
  vector<ElasticState> w(k);  // r .. k
  double d_r;
  //find smoothness indicators for all r = 0,...,k-1 before looping over all states? This can be vectorised.
  //3xElasticState
  vector<ElasticState> beta_r = beta(mat.sol, w.size(), soln_ind); // store 1, 2 and 3.

  for(unsigned int r = 0; r < vr.size(); r++)
  {
    // if statement on right boolean
    if(right)
      d_r = d[r];
    else
      d_r = d[k-1-r];

    for(unsigned int j = 0; j < ElasticState::e_size; j++)
    {
      double alpha_den;
      double d_s;
      for (int s = 0; s < k; s++)
      {
        if(right)
          d_s = d[s];
        else
          d_s = d[k-1-s];

        alpha_den += alpha_r(d[s], beta_r[s][j]); // right or left?
        w[r][j] = alpha_r(d_r, beta_r[r][j])/alpha_den;
      }
    }
  }
  return w;
}

ElasticState WENO_recon(Material& mat, int i, bool right)
{
  const int r = 4, j = 3; // r = -1, 0, 1, 2 :: j = 0, 1, 2, 3
  const double C[r][j] = {
    {11./6., -7./6.,  1./3.} , 
    { 1./3.,  5./6., -1./6.} ,
    {-1./6.,  5./6.,  1./3.} ,
    { 1./3., -7./6., 11./6.} ,
  };
  ElasticState weno;

  unsigned short int k = 3;
  vector<ElasticState> vr(k);
  vector<ElasticState> q = mat.sol;

  unsigned short int Crinc(0);
  if(right)
   Crinc = 1; 

  //(2k - 1) order k = 3 :: 5th order 
  for(unsigned int r = 0; r < k; r++) // r is the row 0 ... k - 1
  {
    for(unsigned int j = 0; j < k; j++)
    {
      //k different reconstructions of v+-1/2
      vr[r] += C[r + Crinc][j] * q[i-r+j]; // check this!
      /* std::cout << "r " <<  r << " j " << j << " " << std::endl; */
    }
  }


  //calculate left and right weights form reconstructions on each variable separately 
  //need r sets of weights (each applied to individual variables)
  vector<ElasticState> w = omega(mat, vr, i, right);

  for(unsigned int r = 0; r < k; r++)
  {
    for(unsigned int j = 0; j < ElasticState::e_size; j++)
    {
      weno[j] += w[r][j]*vr[r][j];
    }
  }
  /* std::cout << weno << std::endl; */
  return weno;
}

ElasticState flux(Material& mat, ElasticState& ql, ElasticState& qr, double dx, double dt);
ElasticState rk3(ElasticState dq_dt, ElasticState q, double dt); //dq_dt, state, dt

void solveXWENO(Material& mat, const double dt)
{
  //From Shu (1998) pg 353 
  // weno constants C_rj
  //  k = r + s + 1 = 3, 
  // j = 0 to j = k - 1
  //k = 3 order is 5 : 2k - 1 
  // r has 4 possible values, but we iterate from 0 to 2
  //for k = 3 
  double dx = mat.dom.dx; //bug
  vector<ElasticState> vl(mat.dom.GNi); // left reconstruction v_i-1/2(+) : 
  vector<ElasticState> vr(mat.dom.GNi); // right reconstruction v_i+1/2(-) : 
  vector<ElasticState> dv_dt(mat.dom.GNi);

  //assuming boundary conditions already applied.
  for (int i = mat.dom.starti-1; i < mat.dom.endi+1; i++)
  {
    //compute reconstructions
    vr[i] = WENO_recon(mat, i, true);
    /* vl[i] = WENO_recon(mat, i, false); */
  }

  /* outputAll("/home/raid/ma595/solid-1D/output/vl", vl); */
  outputAll("/home/raid/ma595/solid-1D/output/vr", vr);
  exit(1);

  ElasticState vlr, vml, vmr, vrl;
  ElasticState fl, fr;
 
  for (int i = mat.dom.starti; i < mat.dom.endi; i++)
  {
    vlr = vr[i-1];
    vml = vl[i];
    vmr = vr[i];
    vrl = vl[i+1];
    /* std::cout << qlr << " " << qml << " " << qmr << " " << qrl << std::endl; */

    /* std::cout << mat.sys.flux(qlr) << std::endl; */
    /* std::cout << mat.sys.flux(qml) << std::endl; */
    /* double a = 5500; */
    /* ElasticState fluxState =  mat.sys.flux(qlr) - mat.sys.flux(qml) + a * (qlr - qrl); */ 
    fl = flux(mat, vlr, vml, dx, dt);
    fr = flux(mat, vmr, vrl, dx, dt);

    dv_dt[i] = (fl - fr)/dx; // rearranged to avoid implementing additional operators // this is the derivative l in J's code.
  }


  // rk3 time-stepping
 
  /* #pragma omp parallel for schedule(dynamic) */
  for (int i = mat.dom.starti; i < mat.dom.endi; i++)
  {
    // test with normal euler time integration // what is the flux?
    mat.sol[i] = rk3(dv_dt[i], mat.sol[i], dt); 
  }


  BCs(mat);
  
}
  //need to do a similar thing as ppm: (RK timestepping)
  /* for (int i = mat.dom.starti; i < mat.dom.endi; i++) */
  /* { */
  /*   double dt_dx = dt/dx; */
  /*   ElasticState f_lf = (mat.sys.flux(vl) + mat.sys.flux(vr) + 1./dt_dx * (vl - vr))/2.; */
  /*   ElasticState v_ri = 0.5*(vl + vr) + 0.5*dt_dx*(mat.sys.flux(vl) - mat.sys.flux(vr)); */
  /*   return 0.5 * (f_lf + mat.sys.flux(v_ri)); */

  /* } */
  //compute force flux using WENO reconstructions 

//need following functions:
//ppm
//--Q (Quadratic reconstruction)
//--P (new reconstruction
//--theta
//--minmod
//SolveXPPM
//flux
//step

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T> double minmod(T a, T b)
{
  return (sgn(a) - sgn(b))*min(abs(a), abs(b))/2.;
}

//we don't need k here, all that's necessary is to use the state that we created 
/* double L(k,j,x) */ 
//If I was more intelligent, I could do this vectorwise instead of componentwise

/* for(int j = 0; j < ElasticState::e_size; j++) */

//do it component wise first:
//limiter 
//i is the integer we use to index the state
/* double theta(Material& mat, const int i) */
/* { */
  
/* } */

/* //reconstruction */
/* ElasticState P(Material) */
/* { */

/* } */

/* double L(int k, int j,) */  
/* { */
  
/* } */

/* //k is the state index */
/* //j is the solution index */
// this can eventually be extended to be vectorised
//
// if using pointers, would have to dereference all the time
// we can potentially speed all this up

double Q(const Material& mat, int k, int j, double x)
{
  vector<ElasticState> q = mat.sol;
  double dx = mat.dom.dx;
  //partial derivatives:
  double dd = (q[j+1][k] - 2.*q[j][k] + q[j-1][k])/pow(dx,2.); //second order derivative
  double d2 = (q[j+1][k] - q[j-1][k])/(2.*dx); //
  return q[j][k] - pow(dx,2.)*dd/24. + d2*x + dd*pow(x,2.)/2.; //this x can be +/- dx/2
}

double L(const Material& mat, int k, int j, double x)
{
  vector<ElasticState> q = mat.sol;
  const double dx = mat.dom.dx;
  return q[j][k] + x*minmod((q[j][k] - q[j-1][k])/dx, (q[j+1][k] - q[j][k])/dx);
}

//this is the reconstruction
ElasticState P(const Material& mat, int j, double x)
{
  ElasticState th = theta(mat, j);
  ElasticState qn;
  for(unsigned int k = 0; k < ElasticState::e_size; k++)
  {
    qn[k] = (1.0 - th[k])*L(mat,k,j,x) + th[k]*Q(mat,k,j,x);
  }
  return qn; 
}

// could eventually turn this function into something that is more lengthy 
ElasticState theta(const Material& mat, int j)
{
  vector<ElasticState> q = mat.sol;
  double xl = -mat.dom.dx/2.;
  double xr = mat.dom.dx/2.;

  ElasticState theta;

  for(unsigned int k = 0; k < ElasticState::e_size; k++)
  {
    if(q[j-1][k] < q[j][k] && q[j][k] < q[j+1][k]) //a < b && b < c
    {
      double Qkjr = Q(mat,k,j,xr); double Lkjr = L(mat,k,j,xr);
      double Qkjl = Q(mat,k,j,xl); double Lkjl = L(mat,k,j,xl);

      double b1 = max(Qkjr, Qkjl) - Lkjr;

      if(b1 == 0.0)
        theta[k] = 1.0;

      double b2 = max(Qkjr, Qkjl) - Lkjl;

      if(b2 == 0.0)
        theta[k] = 1.0;

      if(b1 != 0 || b2 != 0){
        double t1 = max((L(mat,k,j,xl) + L(mat,k,j+1,xl))/2., Q(mat,k,j+1,xl)) - Lkjr;
        double t2 = min((L(mat,k,j,xr) + L(mat,k,j-1,xr))/2., Q(mat,k,j-1,xr)) - Lkjl;
        theta[k] = min(t1/b1,min(t2/b2,1.0));
      }
    }
    else if(q[j-1][k] > q[j][k] && q[j][k] > q[j+1][k])
    {
      double Qkjr = Q(mat,k,j,xr); double Lkjr = L(mat,k,j,xr);
      double Qkjl = Q(mat,k,j,xl); double Lkjl = L(mat,k,j,xl);

      double b1 = max(Qkjr, Qkjl) - Lkjl;

      if(b1 == 0.0)
        theta[k] = 1.0;

      double b2 = max(Qkjr, Qkjl) - Lkjr;

      if(b2 == 0.0)
        theta[k] = 1.0;

      if(b1 != 0 || b2 != 0){
        double t1 = max((L(mat,k,j,xr) + L(mat,k,j-1,xr))/2., Q(mat,k,j-1,xr)) - Lkjl;
        double t2 = min((L(mat,k,j,xl) + L(mat,k,j+1,xl))/2., Q(mat,k,j+1,xl)) - Lkjr;
        theta[k] = min(t1/b1,min(t2/b2,1.0));
      }
    }
    else 
      theta[k] = 0.0;
  }
  return theta;
}

ElasticState flux(Material& mat, ElasticState& ql, ElasticState& qr, double dx, double dt)
{
  /* double a = 6000; // get this from material information */
  /* std::cout << dx/dt << std::endl; */
  /* return (mat.sys.flux(ql) + mat.sys.flux(qr) + a * (ql - qr))/2.; */

  //compute force flux 
  double dt_dx = dt/dx;
  ElasticState f_lf = (mat.sys.flux(ql) + mat.sys.flux(qr) + 1./dt_dx * (ql - qr))/2.;
  ElasticState q_ri = 0.5*(ql + qr) + 0.5*dt_dx*(mat.sys.flux(ql) - mat.sys.flux(qr));
  return 0.5 * (f_lf + mat.sys.flux(q_ri));
}

ElasticState rk1(ElasticState dq_dt, ElasticState q, double dt)
{
  ElasticState q1 = q + dt*dq_dt;
  return q1;
}

ElasticState rk3(ElasticState dq_dt, ElasticState q, double dt) //dq_dt, state, dt
{
  ElasticState q1 =    q                +    dt*dq_dt;
  ElasticState q2 = 3.*q/4. +     q1/4. +    dt*dq_dt/4.;
  ElasticState q3 =    q/3. +  2.*q2/3. + 2.*dt*dq_dt/3.;
  return q3;
}

//runge kutta time stepping?
void solveXPPM(Material& mat, const double dt) //evolve function
{
  double dx = mat.dom.dx; //bug
  double xl = -dx/2.;
  double xr = dx/2.;

  vector<ElasticState> ql(mat.dom.GNi);
  vector<ElasticState> qr(mat.dom.GNi);
  vector<ElasticState> dq_dt(mat.dom.GNi);

  //apply boundaries
  //need to get the number of ghost cells correct here.
  //3:104
  /* #pragma omp parallel for schedule(dynamic) */
  for (int i = mat.dom.starti-1; i < mat.dom.endi+1; i++) 
  {
    ql[i] = P(mat,i,xl); //left reconstruction
    qr[i] = P(mat,i,xr); //left reconstruction
  }


#ifdef debugrun_
  outputAll("/home/raid/ma595/solid-1D/output/ql", ql);
  outputAll("/home/raid/ma595/solid-1D/output/qr", qr);
  outputAll("/home/raid/ma595/solid-1D/output/barton1D", mat.sol);
#endif

  //compute fluxes between cells
  //4:103 (3 ghost cells either side
  ElasticState qlr, qml, qmr, qrl;
  ElasticState fl, fr;
 
  for (int i = mat.dom.starti; i < mat.dom.endi; i++)
  {
    qlr = qr[i-1];
    qml = ql[i];
    qmr = qr[i];
    qrl = ql[i+1];
    /* std::cout << qlr << " " << qml << " " << qmr << " " << qrl << std::endl; */

    /* std::cout << mat.sys.flux(qlr) << std::endl; */
    /* std::cout << mat.sys.flux(qml) << std::endl; */
    /* double a = 5500; */
    /* ElasticState fluxState =  mat.sys.flux(qlr) - mat.sys.flux(qml) + a * (qlr - qrl); */ 
    fl = flux(mat, qlr, qml, dx, dt);
    fr = flux(mat, qmr, qrl, dx, dt);

    dq_dt[i] = (fl - fr)/dx; // rearranged to avoid implementing additional operators // this is the derivative l in J's code.
  }

  // rk3 time-stepping
 
  /* #pragma omp parallel for schedule(dynamic) */
  for (int i = mat.dom.starti; i < mat.dom.endi; i++)
  {
    // test with normal euler time integration // what is the flux?
    mat.sol[i] = rk3(dq_dt[i], mat.sol[i], dt); 
  }

  BCs(mat);
  
}

void solveXSLIC(Material& mat, const double dt)
{
	int N = mat.dom.GNi;
	vector<ElasticState> left(N);
	vector<ElasticState> right(N);
	const double dt_dX = dt/mat.dom.dx;

  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; ++i) 
	{
		//extrapolated values (bar) pg 514 and slope pg 506.
		const ElasticState slopebar = grad(mat, i);
		ElasticState Cleft = mat.sol[i] - 0.5 * slopebar;
		ElasticState Cright = mat.sol[i] + 0.5 * slopebar;
		const ElasticState Cbar = 0.5 * (dt_dX) * (mat.sys.flux(Cleft) - mat.sys.flux(Cright));
/* Cleft.flux() - Cright.flux()); //change these flux functions */
		left[i] = Cleft + Cbar;
		right[i] = Cright + Cbar;
	}		
	
	left[0] = mat.sol[0];
	right[0] = mat.sol[0]; //is required
	left[1] = mat.sol[1];
	right[1] = mat.sol[1]; //is required
	left[N - 1] = mat.sol[N - 1];
	right[N - 1] = mat.sol[N - 1];	
	left[N - 2] = mat.sol[N - 2];
	right[N - 2] = mat.sol[N - 2];	

	System sys = mat.sys;	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; i++)
	{
	  mat.sol[i] += dt_dX*(forceFlux(sys, left, right, dt_dX, i - 1) - forceFlux(sys, left, right, dt_dX, i));
	}
	//printArray(U);
	BCs(mat);
}

//Force flux calculation pg 512
ElasticState grad(const Material& mat, int i)
{
	//1. calculate r
	const double w = 1; // check this 
	ElasticState num = mat.sol[i] - mat.sol[i-1];
	ElasticState den = mat.sol[i+1] - mat.sol[i];
	ElasticState r;
	ElasticState delta = 0.5*(1+w)*num;// + 0.5*(1-w)*denom;
	//aNew = Ui[i].soundSpeed();	
	//cout << i << endl;	
	for(unsigned int j = 0; j < ElasticState::e_size; j++)
	{	
		if(num[j] == 0 && den[j] == 0)
		{
			r[j] = 1;
		}
		else if(den[j] == 0)
		{
			den[j] = 1.e-10;
			r[j] = num[j]/den[j];			
		}
		else 
			r[j] = num[j]/den[j];		
	}
	//ElasticState r = num/denom;
	ElasticState deltaold = delta;	
	for(unsigned int j = 0; j < ElasticState::e_size; j++)
	{
		const double ksi = slopelim(r[j]);
		//ksi = 1;
		delta[j] *= ksi; 	 
	}
	return delta; 
} 

//need to parallelise
double getMinDt(const Material& mat)
{
  System sys = mat.sys; 
  double sharedmindt = std::numeric_limits<double>::max(); 
#pragma omp parallel
  {
    double mindt = std::numeric_limits<double>::max();
#pragma omp for nowait schedule(dynamic)
    for (int i = mat.dom.starti; i < mat.dom.endi; i++) 
    {
      double Smax = sys.getMaxWaveSpeed(sys.conservativeToPrimitive(mat.sol[i]),0);
      mindt = min(mat.dom.dx/Smax, mindt);
    }	
#pragma omp critical
    {
      sharedmindt = std::min(sharedmindt, mindt);
    }
  }
  return sharedmindt;
}

void printArray(Material mat)
{
	cout.precision(4);
	for(int i = 0; i < 140; i++)
	{
		cout << '-';
	}
	cout << endl;
	for(int i = 0; i < mat.dom.GNi; i++)
	{
		cout << setw(7) << left << i << ' ' << mat.sol[i] << endl;
	}
}


/**
 * Calculate the slope limiter and slope delta_i.
 * Update delta_i.
 * @param Ui Pointer to solution domain Ui
 * @param i Domain index of conserved vector  
 * @return updated slope vector delta_i
 */
// delta ibar slope on pg 509
// slope limiters pg 509 and 510
double slopelim(double r)
{
	//calculate ksi_r
	double ksi_sl = 0;
	//slope limiters:
  /* cout << "Slope limiter " << sl << endl; */ 
	switch(sl)
	{		
		case superbee:
			if(r <= 0)
				ksi_sl = 0;
			else if(r > 0 && r <= 0.5)
				ksi_sl = 2*r;
			else if(r >= 0.5 && r <= 1)
				ksi_sl = 1;
			else if (r >= 1)
				ksi_sl = min(r, min(ksi_r(r), 2.));
			break;
		case vanleer:		
			if(r <= 0)
				ksi_sl = 0;
			else
				ksi_sl = min(2*r/(1+r), ksi_r(r));
			break;
		case albalda:
			if(r <= 0)
				ksi_sl = 0;
			else 
				ksi_sl = min(r*(1+r)/(1+r*r), ksi_r(r));			
			break;
		case minbee:	
			//cout << "minbee"	<< endl;
			if(r <= 0)
				ksi_sl = 0;
			else if(r >= 0 && r <= 0.5)
				ksi_sl = r;
			else 
				ksi_sl = min(1., ksi_r(r));
			break;
	}
	return ksi_sl;			
}

double ksi_r(double r)
{
	/* double cNew = aNew*dt/dX; */
	//cout << "courant Number	" << cNew << endl;
	//return 2./(1. + r);
	double beta_fw = 1.;//2/(0.1);
	double w = 1.;
	return 2.*beta_fw/(1. - w + (1. + w)*r);
}

//this function outputs the conservative state and all boundary conditions
//may be better to include 
void outputAll(string file, const vector<ElasticState> vec)
{
	ofstream output;	
  output.open(file.c_str()); //we append all results to the same file 
  int GCs = 3;
  double dx = 1./(vec.size()-2*GCs);
  double state;
  for (int i = GCs; i < vec.size()-GCs; i++)
  {
	  output << (double)dx*((i-GCs)+0.5);
    ElasticState consState = vec[i];	
    for(unsigned int i = 0 ; i < ElasticState::e_size; i++)
    {
      if(consState[i] < 1e-20)
        state = 0;
      else
        state = consState[i];
      output << '\t' << state;
    }
    output << endl;
  }
  output.close();
}

void outputGnu(string file, Material mat, int outStep, double t)
{
	/* int ret = system(("mkdir -p " + fileName).c_str()); */
	/* if(ret!=0) cerr << "Error creating directory?\n"; */
	cerr << "Writing results to \"" << file << "\"\n";
	ofstream output;	
	/* output.precision(3); */
  if(outStep == 0)
  {
    output.open(file.c_str(), ios::app); //we append all results to the same file 
    output << "# t = 0 " << '\n'
      << "# Column 1: x-coordinate" << '\n'
      << "# Column 2: density" << '\n'
      << "# Column 3: U" << '\n'
      << "# Column 4: V" << '\n'
      << "# Column 5: W" << '\n'
      << "# Column 6: sigma_11" << '\n'
      << "# Column 7: sigma_12" << '\n'
      << "# Column 8: sigma_13" << '\n'
      << "# Column 9: sigma_21" << '\n'
      << "# Column 10: sigma_22" << '\n'
      << "# Column 11: sigma_23" << '\n'
      << "# Column 12: sigma_31" << '\n'
      << "# Column 13: sigma_32" << '\n'
      << "# Column 14: sigma_33" << '\n'
      << "# Column 15: S " << '\n'
      << "# Column 16: I1" << '\n'
      << "# Column 17: I2" << '\n' 
      << "# Column 18: I3" << '\n'
      << "# Column 19: dI_1" << '\n'
      << "# Column 20: dI_2" << '\n'
      << "# Column 21: dI_3" << '\n' 
	    << "# Column 22: G11" << '\n'  
	    << "# Column 22: G12" << '\n'  
	    << "# Column 22: G13" << '\n'  
	    << "# Column 22: G21" << '\n'  
	    << "# Column 22: G22" << '\n'  
	    << "# Column 22: G23" << '\n'  
	    << "# Column 22: G31" << '\n'  
	    << "# Column 22: G32" << '\n'  
	    << "# Column 22: G33" << '\n' << endl;
    output.close();
  }
  //some work needs doing here
	vector <double> out;
  output.open(file.c_str(), ios::app);
  if(outStep != 0)
    output << "# t = " << t << endl;
	for(int i = mat.dom.starti; i < mat.dom.endi; i++)
	{
		//stress (system), velocity (primState), entropy (EOS), invariants
		const ElasticState consState = mat.sol[i];	
		double rho = mat.sys.Density(consState);
		const ElasticPrimState primState = mat.sys.conservativeToPrimitive(consState);
		double entropy = primState.S_();
		Vector3d u = primState.u_();
		Vector3d inv = mat.sys.getInvariants(primState.F_());
		Matrix3d sigma = mat.sys.stress(primState);
  	const Matrix3d G = mat.sys.strainTensor(primState.F_());
    // convert cell into 
    // phys domain is mat.dom.Lx/cell 
    // (i-2)*dx + dx/2
    // dx*((i-mat.dom.starti) + 0.5)
    //get invariants
    
    
		output << (double)mat.dom.dx*((i-mat.dom.starti)+0.5) << '\t'
			<< rho/1e3 << '\t' 
			<< u(0)/1000.<< '\t' 
			<< u(1)/1000.<< '\t' 
			<< u(2)/1000.<< '\t' 
			<< sigma(0,0)/1e9 << '\t'
			<< sigma(0,1)/1e9 << '\t'
			<< sigma(0,2)/1e9 << '\t'
			<< sigma(1,0)/1e9 << '\t'
			<< sigma(1,1)/1e9 << '\t'
			<< sigma(1,2)/1e9 << '\t'
			<< sigma(2,0)/1e9 << '\t'
			<< sigma(2,1)/1e9 << '\t'
			<< sigma(2,2)/1e9 << '\t'
			<< entropy/1e6 << '\t'
			<< inv[0] << '\t'
			<< inv[1] << '\t'
			<< inv[2] << '\t' 
      << G(0,0) << '\t'
      << G(0,1) << '\t'
      << G(0,2) << '\t'
      << G(1,0) << '\t'
      << G(1,1) << '\t'
      << G(1,2) << '\t'
      << G(2,0) << '\t'
      << G(2,1) << '\t'
      << G(2,2) << '\t'
      << endl;

		/* for(unsigned int j = 0; j < out.size(); j++) */
		/* { */
		/* 	output << '\t' << out[j]; */
		/* } */
		/* out.clear(); */
	}
  output << '\n' << std::endl;
	output.close();	
}

void advance(Material& mat, const double dt) //or evolve?
{
  //add the plasticity bit
// apply curl constraint (2D) // geometric -
// cylindrical //spherical bcs //plasticity
  //series of advance functions: levelset, geometric bcs, 
  /* solveXWENO(mat, dt); */
  solveXGodunov(mat, dt);
}


int solveSystem(InputSolid inputSolid, Material* mat){
  
  double t(0), dt(0), tend(inputSolid.end_time);
  const double CFL(inputSolid.input_CFL);
  //output names 
  double outFreq(tend/inputSolid.frequency);
	int step(0), outStep(0); 
  std::string outDir(inputSolid.filePath), outName(inputSolid.fileName);
  std::string outFile(outDir + outName);
  //delete existing output
  ofstream myfile;
  myfile.open(outFile, ios::trunc); // this needs to be kept here (or have as an if statement)
  myfile.close();

	while(t < tend)
	{
		BCs(*mat);
    #ifdef debugrun_
    std::cout << "writing initial condition + BCs" << std::endl;
    outputAll(outFile, (*mat).sol);
    outputAll("/home/raid/ma595/solid-1D/output/initBCs", (*mat).sol);
    #endif
		//1. calculate time step: CFL and boundary conditions pg 495
		dt = getMinDt(*mat);
		/* if(step < 20)	{dt /= 10;} */
		dt *= CFL;

#ifdef debugrun_
    cout << "The calculated time-step is: " << dt << endl;
#endif

    if(t < outStep * outFreq && t + dt > outStep * outFreq) {
      dt = outStep * outFreq - t;
    }
    
    if(t >= outStep * outFreq) 
    {
      outputGnu(outFile, *mat, outStep, outStep*outFreq);
      cerr << "Saved results to output file " << outStep << "\n";
      ++outStep;
    }
		/* if(t > tf) */
		/* { */				
		/* 	t -= dt; */
		/* 	dt = tf - t; */
		/* 	t += dt; // should equal tf */
		/* } */
		cout.precision(4);
		cout << " [" << outStep << "] " << setw(6) << step << setw(6) << "Time " << setw(6) << t << " dt " << setw(6) << dt << 
						setw(15) << " Remaining " << setw(6) <<tend-t<< endl;

    advance(*mat, dt);

		t += dt;
    /* if(step == 1) */
    /* { */
    /*   outputGnu(outFile, *mat, outStep, t); */
    /*   break; */
    /* } */
		++step;
	}

	outputGnu(outFile, *mat, outStep, tend); //final output

  //total time
#ifdef debug_
	printArray(*mat);
#endif
  return 0;
}

void unitTests(Material mat, ElasticPrimState iPrimL, ElasticPrimState iPrimR)
{
	cout << "perform simple solution array tests " << endl;
	/* printArray(mat); */

#ifdef debug_
	cout.precision(6);

	mat.sys.Eos.checkEosConstants();

	cout << "Prim stateL" << endl;
	cout << iPrimL << endl;
	cout << "Density" << mat.sys.Density(iPrimL) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimL.F_()) << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	

	cout << "Prim stateR" << endl;
	cout << iPrimR << endl;
	cout << "Density" << mat.sys.Density(iPrimR) << endl;
	cout << "Strain tensor L " << mat.sys.strainTensor(iPrimR.F_()) << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	

	ElasticState iStateL = mat.sys.primitiveToConservative(iPrimL);
	ElasticState iStateR = mat.sys.primitiveToConservative(iPrimR);

	cout <<"Cons stateL" << endl;
	cout << iStateL << endl;
	cout << "Density" << mat.sys.Density(iStateL) << endl;
	cout << endl;

	cout <<"Cons stateR" << endl;
	cout << iStateR << endl;
	cout << "Density" << mat.sys.Density(iStateR) << endl;
	cout << endl;

	ElasticPrimState iPrimL2 = mat.sys.conservativeToPrimitive(iStateL);
	ElasticPrimState iPrimR2 = mat.sys.conservativeToPrimitive(iStateR);

	cout << "Prim stateL" << endl;
	cout << iPrimL2 << endl;
	cout << "Density" << mat.sys.Density(iPrimL2) << endl;
 	cout << endl;
	cout << "INVARIANTS L " << endl;
	cout << mat.sys.getInvariants(iPrimL.F_()) << endl;	

	cout << "Prim stateR" << endl;
	cout << iPrimR2 << endl;
	cout << "Density" << mat.sys.Density(iPrimR2) << endl;
	cout << endl;
	cout << "INVARIANTS R " << endl;
	cout << mat.sys.getInvariants(iPrimR.F_()) << endl;	

//---------------------
	cout << "Check EOS constants" << endl;
	mat.sys.Eos.checkEosConstants();
#endif
}

void getLimiter(InputSolid inputSolid)
{
  std::string limStr = inputSolid.limiter;
  if(limStr == "superbee") 
  {
    sl = superbee;
  }
  else if(limStr == "minbee") 
  {
    sl = minbee;
  }
  else if(limStr == "vanleer")
  {
    sl = vanleer;
  }
  else if(limStr == "albalda")
  {
    sl = albalda;
  }
  else 
  {
    cerr << "Unknown limiter \'" << sl << "\' - options are superbee, minbee, vanleer, albalda\n";
    exit(1);
  }
}

int main(int argc, char ** argv)
{
  // Enable floating point error checking
  /* std::vector<ElasticState> solnVector; */
  /* int i = 0; */
  /* int j = 0; */
  /* double output = solnVector[i+1][j] - 1.0; */

  /* std::cout << output << std::endl; */
  feenableexcept(FE_INVALID);                                                                                                                                                                                    
  feenableexcept(FE_DIVBYZERO);

  char* icStr = argv[1];
  printf("Settings file: %s\n", icStr);
  InputSolid inputSolid;
  inputSolid.readConfigFile(icStr);

  //create and initialize domain
	Domain dom = Domain(inputSolid.cellCountX, 3, inputSolid.xMax);  //changed to 3 ghost cells for PPM
  ElasticEOS Eos(inputSolid.matL);
  Material* mat = new Material(inputSolid.matL, dom, Eos);

  //create left and right states and initialise
  ElasticPrimState primStateL(inputSolid.uL, inputSolid.FL, inputSolid.SL);
  ElasticPrimState primStateR(inputSolid.uR, inputSolid.FR, inputSolid.SR);
  double iface = inputSolid.iface;
  ICInterface(*mat, iface, primStateL, primStateR);
  outputAll("/home/raid/ma595/solid-1D/output/initial", (*mat).sol);
  //get limiter
  getLimiter(inputSolid);

  #ifdef debug_
  unitTests(*mat, primStateL, primStateR);
  exit(1); //return 0
  #endif
  //solvesystem
  double begin = omp_get_wtime();
  solveSystem(inputSolid, mat);
  double end = omp_get_wtime();
  double elapsed_secs = double(end - begin);
  std::cout << "TIME " << elapsed_secs << std::endl;
  delete mat;
  return 0;
}



