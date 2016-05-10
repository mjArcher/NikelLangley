/* #include <config4cpp/Configuration.h> */
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <libconfig.h++>

using namespace libconfig;
/* using namespace config4cpp; */
using namespace Eigen;
using namespace std;

void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting&);

//originally in main
/* string getDirName(const string& prefix, const Domain& dom) */ 
/* { */
/*   stringstream filenameSS; */
/*   filenameSS << "out/" << prefix << "_" << dom.Ni << "/"; */
/*   return filenameSS.str(); */
/* } */    




int readFromFile()
{
	Config cfg;
	Matrix3d F_t;
	Vector3d u_t;
	double S_t, rho_t, e_t;
	try
	{
		/* cfg.readFile("/lsc/zeushome/ma595/Code/computing/solid/solid-1D/input/barton1D.cfg"); */
		cfg.readFile("./barton1D.cfg");
	}
	catch(const FileIOException &fioex)
	{
		std::cerr << "I/O error while reading file." << std::endl;
		return(EXIT_FAILURE);
	}
	catch(const ParseException &pex)
	{
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
			<< " - " << pex.getError() << std::endl;
		return(EXIT_FAILURE);
	}
	try
	{
    std::cout << "Reading from configuration file" << std::endl;
    /* int slLimChoice; */
		/* cfg.lookupValue("slChoice", slLimChoice); */
    /* std::cout << "slChoice " << slLimChoice << std::endl; */
    /* cfg.lookupValue("base", baseDir); */
    // can either use lookup here or index it as shown below 
		Setting &settingState = cfg.lookup("simulation.riemannstate.left");
    Setting& root = cfg.getRoot();
    Setting& simulation = root["simulation"];
		/* cout << "Length " << settingState.getLength() << endl; */
		/* int temp = settingState["rho"].getLength(); */
    int max = 3;
    int count = 3;
    Matrix3d F;
    for (unsigned int i = 0; i < max; i++)
    {
      for (unsigned int j = 0; j < max; j++)
      {
        count = i*max + j;
        F(i,j) = simulation["riemannstate"]["left"]["F"][count];
        std::cout << F(i,j) << std::endl;
      }
    }

		/* /1* /2* //get leftstate *2/ *1/ */
		/* getState(u_t, F_t, rho_t, S_t, settingState); */
  



		/* /1* /2* e_t = internalEnergy(F_t, S_t); *2/ *1/ */	
    /* /1* double e_t = 1; *1/ */
		/* /1* /2* ElasticState ElasticStateL(u_t, F_t, rho_t, S_t, e_t); *2/ *1/ */
		/* /1* /2* cout << ElasticStateL << endl; *2/ *1/ */	
		/* /1* settingState.lookupValue("slLimChoice", rho_t); *1/ */
		/* /1* cout << rho_t << endl; *1/ */

		/* cout << "slChoice" <<slLimChoice << endl; */
	}
	catch(const SettingNotFoundException &nfex)
	{
		cerr << "No 'name' setting in configuration file." << endl;
	}
}

// get settings from configuration file
void getState(Vector3d &u, Matrix3d &F, double &rho, double& S, Setting& setting)
{
	setting.lookupValue("rho", rho);
	setting.lookupValue("S", S);
	int max = setting["u"].getLength();
	double temp = 0; 
	for (unsigned int i = 0; i < max; i++)
	{
		u(i) = setting["u"][i];
    std::cout << u(i) << std::endl;
	}
	int count = 0;
	for (unsigned int i = 0; i < max; i++)
	{
		for (unsigned int j = 0; j < max; j++)
		{
			count = i*max + j;
			F(i,j) = setting["F"][count];
      std::cout << F(i,j) << std::endl;
		}
	}
}

int main(int argc, char ** argv)
{
  std::cout << "Hello world" << std::endl;
  readFromFile();
}

