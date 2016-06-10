#include <libconfig.h++>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <omp.h>
#include <fenv.h>
#include <ctime>

#include "ElasticState.h"
#include "ElasticPrimState.h"
#include "SolidSystem.h"
#include "InputSolid.h"
#include "Utils.h"
/* #include "Domain.hpp" */

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

typedef Eigen::Matrix<ElasticPrimState, Eigen::Dynamic, Eigen::Dynamic> SolnArryPrim;
typedef Eigen::Matrix<ElasticState, Eigen::Dynamic, Eigen::Dynamic> SolnArry;

struct Domain {
	int Ni; //number of y cells
	int GNi; // number of y cells plus ghost cells
	int Nj; //number of x cells
	int GNj; // number of x cells plus ghost cells
	int GC; // number of ghost cells
	int starti; //start of domain
	int endi; //end of domain
	int startj; //start of domain
	int endj; //end of domain
	double Lx; // length of domain
	double Ly; // length of domain
	double dx; // cell size
	double dy; // cell size

  Domain(const int a_Ni, const int a_Nj, const int a_GC,
      const double a_Lx) :
    Ni(a_Ni), Nj(a_Nj), GNi(a_Ni+2*a_GC), GNj(a_Nj+2*a_GC), GC(a_GC),
    starti(a_GC), startj(a_GC), endi(a_Ni+a_GC), endj(a_Nj+a_GC), 
    Lx(a_Lx), Ly((a_Lx*a_Nj)/a_Ni), dx(a_Lx/a_Ni), dy(a_Lx/a_Ni) {}
  Domain(){}
	~Domain(){
	}
};

struct Material {
  const System sys;
  /* const Domain<ElasticState> dom; */
  Domain dom;
  SolnArry sol;
  string name;

  Material(const string& a_name, const Domain a_dom, const ElasticEOS Eos):
    sys(Eos), dom(a_dom), sol(a_dom.GNi, a_dom.GNj), name(a_name) {}
};

double slopelim(double);
double ksi_r(double);
ElasticState grad(const ElasticState, const ElasticState, const ElasticState);
/* void outputAll(string, const ); */

void BCs(Material& mat) {  
  //transmissive boundary conditions 
  //apply boundary conditions row wise first 
  for(int j=mat.dom.startj; j<mat.dom.endj; ++j){
    //start
#pragma omp parallel for schedule(dynamic)
    for(int i=0; i<mat.dom.starti; ++i) {
      const int imagei = 2*mat.dom.starti-1-i;
      mat.sol(i,j) = mat.sol(imagei,j);
    }

    //end
#pragma omp parallel for schedule(dynamic)
    for(int i=mat.dom.endi; i<mat.dom.GNi; ++i) {
      const int imagei = 2*mat.dom.endi-1-i;
      mat.sol(i,j) = mat.sol(imagei,j);
    }  
  }

  //loop over rows
  for(int i=mat.dom.starti; i<mat.dom.endi; ++i){
    //start
#pragma omp parallel for schedule(dynamic)
    for(int j=0; j<mat.dom.startj; ++j) {
      const int imagej = 2*mat.dom.startj-1-j;
      mat.sol(i,j) = mat.sol(i,imagej);
    }

    //end
#pragma omp parallel for schedule(dynamic)
    for(int j=mat.dom.endj; j<mat.dom.GNj; ++j) {
      const int imagej = 2*mat.dom.endj-1-j;
      mat.sol(i,j) = mat.sol(i,imagej);
    }  
  }
}

//rotate states
void ICRiemann2DTrivial(Material& mat,
    const double iface,
    const ElasticPrimState& left, const ElasticPrimState& right) 
{
  //convert to conserved states
  ElasticState stateL = mat.sys.primitiveToConservative(left);
  ElasticState stateR = mat.sys.primitiveToConservative(right);

  //trivial Riemann problem
  for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
    for(int j=mat.dom.startj; j<mat.dom.endj; ++j) {
      const double x = (j-mat.dom.startj+0.5)*mat.dom.dx;
      if(x<iface) 
      {      
        mat.sol(i,j) = stateL;
      }
      else 
      {
        mat.sol(i,j) = stateR;
      } 
    }
  }
  cout << "Initialised 2D domain " << endl;
}

//general 2D Riemann initialisation 
void ICRiemann2DAngle(Material& mat, double iface, double theta, const ElasticPrimState& left, const ElasticPrimState& right) 
{
  //  
  ElasticState stateL = mat.sys.primitiveToConservative(left);
  ElasticState stateR = mat.sys.primitiveToConservative(right);

  cout << mat.sys.Density(left) << endl;

  cout << mat.sys.Density(right) << endl;

  //change in i alters the position of the interface 
  //
  /* cout <<  << endl; */
  double ifaceInit = 0.01 + mat.dom.dx/2.;

  for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
    double y = (i-mat.dom.starti+0.5)*mat.dom.dy;
    iface = ifaceInit - y*tan(theta*M_PI/180.);
    for(int j=mat.dom.startj; j<mat.dom.endj; ++j) {
      double x = (j-mat.dom.startj+0.5)*mat.dom.dx;
      if(x < iface)
        mat.sol(i,j) = stateL;
      else
        mat.sol(i,j) = stateR;
    }
  }

}

//45 degree angle
void ICRiemann2D(Material& mat, const ElasticPrimState& left, const ElasticPrimState& right) 
{
  //convert to conserved states
  ElasticState stateL = mat.sys.primitiveToConservative(left);
  ElasticState stateR = mat.sys.primitiveToConservative(right);

  cout << mat.sys.Density(left) << endl;

  cout << mat.sys.Density(right) << endl;

  double iface = 0;
  //trivial Riemann problem
  for(int i=mat.dom.starti; i<mat.dom.endi; ++i) {
    for(int j=mat.dom.startj; j<mat.dom.endj; ++j) {
      const double x = (j-mat.dom.startj+0.5)*mat.dom.dx;
      if(i<=mat.dom.endj - j) 
      {      
        mat.sol(i,j) = stateL;
      }
      else 
      {
        mat.sol(i,j) = stateR;
      } 
    }
  }
  cout << "Initialised 2D domain " << endl;
}

//supply left state 
//this will affect the calculation of the timestep
/* void solveXGodunov(Material& mat, const double dt) */
/* { */
/* 	const double dt_dX = dt/mat.dom.dx; */
/*   const vector<ElasticState> soln = mat.sol; */
/*   #pragma omp parallel for schedule(dynamic) */
/* 	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; i++) */
/* 	{ */
/* 	  mat.sol[i] += dt_dX*(mat.sys.godunovFlux(soln[i-1], soln[i]) - mat.sys.godunovFlux(soln[i], soln[i+1])); */
/* 	} */
/*   /1* exit(1); *1/ */
/* 	//printArray(U); */
/* 	BCs(mat); */
/* } */

ElasticState forceFlux(System& sys, vector<ElasticState>& left, vector<ElasticState>& right, double dt_dX, int i,  const int dirn)
{
	const ElasticState F_lf = 0.5*(sys.flux(right[i], dirn) + sys.flux(left[i+1], dirn)) + (0.5/dt_dX)*(right[i]-left[i+1]);
	const ElasticState C_r = 0.5*(right[i] + left[i+1]) + 0.5*(dt_dX)*(sys.flux(right[i], dirn) - sys.flux(left[i+1], dirn));
	ElasticState F_force = 0.5*(F_lf + sys.flux(C_r, dirn));
	return F_force;
}

template <typename T> int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T> double minmod(T a, T b)
{
  return (sgn(a) - sgn(b))*min(abs(a), abs(b))/2.;
}


//row or col(could pass the entire matrix), dt, Nx, Ny, dirn
//matrix, dt, dir, row/col val
void solveXSLIC(Material& mat, const double dt, const int row)
{
	int N = mat.dom.GNj;
	vector<ElasticState> left(N);
	vector<ElasticState> right(N);
	const double dt_dX = dt/mat.dom.dx;

  #pragma omp parallel for schedule(dynamic)
	for(int j = mat.dom.startj - 1; j < mat.dom.endj + 1; ++j)
	{
		//extrapolated values (bar) pg 514 and slope pg 506.
		const ElasticState slopebar = grad(mat.sol(row,j), mat.sol(row,j+1), mat.sol(row,j-1));
    /* cout << slopebar << endl; */
    /* double slopebarval = slopebar[0]; */
    /* cout << slopebarval << endl; */
		ElasticState Cleft = mat.sol(row,j) - 0.5*slopebar;
		ElasticState Cright = mat.sol(row,j) + 0.5*slopebar;
    /* if(Cleft != Cright) */
    /*   cout << Cleft << "\n" << Cright << endl; */
    /* cout << mat.sol(row,j+1) << endl; */
		const ElasticState Cbar = 0.5*(dt_dX)*(mat.sys.flux(Cleft,0) - mat.sys.flux(Cright,0));
/* Cleft.flux() - Cright.flux()); //change these flux functions */
		left[j] = Cleft + Cbar;
		right[j] = Cright + Cbar;
	}		
	
	left[0] = mat.sol(row,0);
	right[0] = mat.sol(row,0); //is required
	left[1] = mat.sol(row,1);
	right[1] = mat.sol(row,1); //is required
	left[N-1] = mat.sol(row,N-1);
	right[N-1] = mat.sol(row,N-1);	
	left[N-2] = mat.sol(row,N-2);
	right[N-2] = mat.sol(row,N-2);	

	System sys = mat.sys;	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
  #pragma omp parallel for schedule(dynamic)
	for(int j = mat.dom.startj-1; j < mat.dom.endj+1; j++)
	{
	  mat.sol(row,j) += dt_dX*(forceFlux(sys, left, right, dt_dX, j-1, 0) - forceFlux(sys, left, right, dt_dX, j, 0));
	}
	//printArray(U);
	BCs(mat);
}

void solveYSLIC(Material& mat, const double dt, const int col)
{
	int N = mat.dom.GNi;
	vector<ElasticState> left(N);
	vector<ElasticState> right(N);
	const double dt_dX = dt/mat.dom.dx;

  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti - 1; i < mat.dom.endi + 1; ++i)
	{
		//extrapolated values (bar) pg 514 and slope pg 506.
		const ElasticState slopebar = grad(mat.sol(i,col), mat.sol(i+1,col), mat.sol(i-1,col));
		ElasticState Cleft = mat.sol(i,col) - 0.5*slopebar;
		ElasticState Cright = mat.sol(i,col) + 0.5*slopebar;
		const ElasticState Cbar = 0.5*(dt_dX)*(mat.sys.flux(Cleft,1) - mat.sys.flux(Cright,1));
		left[i] = Cleft + Cbar;
		right[i] = Cright + Cbar;
	}		
	
	left[0] = mat.sol(0,col);
	right[0] = mat.sol(0,col); //is required
	left[1] = mat.sol(1,col);
	right[1] = mat.sol(1,col); //is required
	left[N-1] = mat.sol(N-1,col);
	right[N-1] = mat.sol(N-1,col);	
	left[N-2] = mat.sol(N-2,col);
	right[N-2] = mat.sol(N-2,col);	

	System sys = mat.sys;	
	//3. calculate force flux using LF and RI, and calculate new cell averaged Ui pg 494
  #pragma omp parallel for schedule(dynamic)
	for(int i = mat.dom.starti-1; i < mat.dom.endi + 1; i++)
	{
	  mat.sol(i,col) += dt_dX*(forceFlux(sys, left, right, dt_dX, i-1, 0) - forceFlux(sys, left, right, dt_dX, i, 0));
	}
	//printArray(U);
	BCs(mat);
}

//Force flux calculation pg 512
ElasticState grad(const ElasticState cons_i, const ElasticState cons_ip1, const ElasticState cons_im1)
{
  //opportunities to speed up here (check for equality)
  //i,i+1,i-1
	//1. calculate r
	const double w = 1; // check this 
	ElasticState num = cons_i - cons_im1;
	ElasticState den = cons_ip1 - cons_i;
	ElasticState r;
	ElasticState delta = 0.5*(1+w)*num;// + 0.5*(1-w)*denom;
	
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

double getMinDt(const Material& mat)
{
  System sys = mat.sys; 
  double mindt = std::numeric_limits<double>::max();
  for (int i = mat.dom.starti; i < mat.dom.endi; i++) 
  {
    for(int j = mat.dom.startj; j < mat.dom.endj; j++)
    {
      ElasticPrimState prim = sys.conservativeToPrimitive(mat.sol(i,j));
      double Smax_x = sys.getMaxWaveSpeed(prim,0);
      double Smax_y = sys.getMaxWaveSpeed(prim,1);
      mindt = min(mat.dom.dx/Smax_x, mindt);
      mindt = min(mat.dom.dy/Smax_y, mindt);
    }	
  }
  return mindt;
}

/* void printArray(Material mat) */
/* { */
/* 	cout.precision(4); */
/* 	for(int i = 0; i < 140; i++) */
/* 	{ */
/* 		cout << '-'; */
/* 	} */
/* 	cout << endl; */
/* 	for(int i = 0; i < mat.dom.GNi; i++) */
/* 	{ */
/* 		cout << setw(7) << left << i << ' ' << mat.sol[i] << endl; */
/* 	} */
/* } */


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
/* void outputAll(string file, const Material mat) */
/* { */
/* 	ofstream output; */	
/*   output.open(file.c_str()); //we append all results to the same file */ 
/*   int GCs = 2; */


/*   //this all needs to be changed */ 

/*   double dx = 1./(vec.size()-2*GCs); */
/*   double state; */
/*   for (int i = GCs; i < vec.size()-GCs; i++) */
/*   { */
/* 	  output << (double)dx*((i-GCs)+0.5); */
/*     ElasticState consState = matrix(i, j); */	
/*     for(unsigned int i = 0 ; i < ElasticState::e_size; i++) */
/*     { */
/*       if(consState[i] < 1e-20) */
/*         state = 0; */
/*       else */
/*         state = consState[i]; */
/*       output << '\t' << state; */
/*     } */
/*     output << endl; */
/*   } */
/*   output.close(); */
/* } */

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
	/* for(int i = mat.dom.starti; i < mat.dom.endi; i++) */
	/* { */
    /* for(int j = mat.dom.startj; j < mat.dom.endj; j++) */
    /* { */
  for(int i = 0; i < mat.dom.GNi; i++)
  {
    for(int j = 0; j < mat.dom.GNj; j++)
    {
      //stress (system), velocity (primState), entropy (EOS), invariants
      const ElasticState consState = mat.sol(i,j);	
      double rho = mat.sys.Density(consState);
      /* const ElasticPrimState primState = mat.sys.conservativeToPrimitive(consState); */
      /* double entropy = primState.S_(); */
      /* Vector3d u = primState.u_(); */
      /* Vector3d inv = mat.sys.getInvariants(primState.F_()); */
      /* Matrix3d sigma = mat.sys.stress(primState); */
      /* const Matrix3d G = mat.sys.strainTensor(primState.F_()); */
    
      //output row wise
      output << (double)mat.dom.dx*((j-mat.dom.startj)+0.5) << '\t' << (double)mat.dom.dy*((i-mat.dom.starti)+0.5) << '\t' << rho/1e3 << endl;

       

        /* << u(0)/1000.<< '\t' */ 
        /* << u(1)/1000.<< '\t' */ 
        /* << u(2)/1000.<< '\t' */ 
        /* << sigma(0,0)/1e9 << '\t' */
        /* << sigma(0,1)/1e9 << '\t' */
        /* << sigma(0,2)/1e9 << '\t' */
        /* << sigma(1,0)/1e9 << '\t' */
        /* << sigma(1,1)/1e9 << '\t' */
        /* << sigma(1,2)/1e9 << '\t' */
        /* << sigma(2,0)/1e9 << '\t' */
        /* << sigma(2,1)/1e9 << '\t' */
        /* << sigma(2,2)/1e9 << '\t' */
        /* << entropy/1e6 << '\t' */
        /* << inv[0] << '\t' */
        /* << inv[1] << '\t' */
        /* << inv[2] << '\t' */ 
        /* << G(0,0) << '\t' */
        /* << G(0,1) << '\t' */
        /* << G(0,2) << '\t' */
        /* << G(1,0) << '\t' */
        /* << G(1,1) << '\t' */
        /* << G(1,2) << '\t' */
        /* << G(2,0) << '\t' */
        /* << G(2,1) << '\t' */
        /* << G(2,2) << '\t' */
        /* << endl; */

      /* for(unsigned int j = 0; j < out.size(); j++) */
      /* { */
      /* 	output << '\t' << out[j]; */
      /* } */
      /* out.clear(); */
    }
    output << "\n";
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
  /* solveXGodunov(mat, dt); */

  //x sweeps
  for(int i = mat.dom.starti; i <mat.dom.endi; i++)
  {
    solveXSLIC(mat, dt, i);
  }
  //y sweeps
  for(int j = mat.dom.startj; j <mat.dom.endj; j++)
  {
    solveYSLIC(mat, dt, j);
  }

}


int solveSystem(InputSolid inputSolid, Material* mat)
{
  double t(0), dt(0), tend(inputSolid.end_time);
  const double CFL(inputSolid.input_CFL);
  //output names 
  double outFreq(tend/inputSolid.frequency);
	int step(0), outStep(0); 
  std::string outDir(inputSolid.filePath), outName(inputSolid.fileName);
  std::string outFile(outDir + outName);
  //delete existing output
  ofstream myfile;
  cout << "Output file " << outFile << endl;
  myfile.open(outFile, ios::trunc); // this needs to be kept here (or have as an if statement)
  myfile.close();
  //initial output
  outputGnu(outFile, *mat, 0, 0);

	while(t < tend)
	{
		BCs(*mat);
    #ifdef debugrun_
    /* std::cout << "writing initial condition + BCs" << std::endl; */
    /* outputAll(outFile, (*mat).sol); */
    /* outputAll("/home/raid/ma595/solid-1D/output/initBCs", (*mat).sol); */
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
      /* std::string out(outFile + "_" + convertToStr(outStep)); */
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
    /* std::string out(outFile + "_" + convertToStr(outStep)); */
    /* outputGnu(outFile, *mat, 1, 10.); */
    /* exit(1); */
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
  //2D domain: x corresponds to the j direction, y corresponds to the i direction 
	Domain dom = Domain(inputSolid.cellCountY, inputSolid.cellCountX, 2, inputSolid.xMax);  //rows, cols, GC, 
  cout << "input material " << inputSolid.matL << endl;
  ElasticEOS Eos(inputSolid.matL);
  Material* mat = new Material(inputSolid.matL, dom, Eos);
  //create left and right states and initialise
  ElasticPrimState primStateL(inputSolid.uL, inputSolid.FL, inputSolid.SL);
  ElasticPrimState primStateR(inputSolid.uR, inputSolid.FR, inputSolid.SR);
  double iface = inputSolid.iface;
  double theta = 30;
  ICRiemann2DAngle(*mat, iface, theta, primStateL, primStateR);
  /* cout << (*mat).sol(100,100) << endl; */
  std::string outDir(inputSolid.filePath), outName(inputSolid.fileName);
  std::string outFile(outDir + outName);
  outputGnu(outFile, *mat, 0, 0);
  //get limiter
  getLimiter(inputSolid);

  #ifdef debug_
  unitTests(*mat, primStateL, primStateR);
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



