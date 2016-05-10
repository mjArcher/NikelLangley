#include <libconfig.h++>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace libconfig;

//read in and store simulation parameters from cfg file
class InputSolid{
  public:
    int test;

    bool cartesian;
    bool curvilinear; 
    bool cut_cell;

    bool WAF;
    bool MUSCL;
    bool SLIC;
    bool WENO;
    bool FORCE;

    bool start_output;
    bool end_output;

    bool time;

    double interval;
    int frequency;

    // currently only 2D 
    double xMin;
    double xMax;
    double yMin;
    double yMax;

    double deltaX;
    double deltaY;

    int cellCountX;
    int cellCountY;

    double end_time;

    // solid parameters
    std::string matL;
    Eigen::Vector3d uL;
    Eigen::Matrix3d FL;
    double SL;

    std::string matR;
    Eigen::Vector3d uR;
    Eigen::Matrix3d FR;
    double SR;

    //cfl
    double input_CFL;
    double iface;

    std::string filePath;
    std::string fileName;
    std::string geometryPath;

    std::string limiter;

    InputSolid(){
      time = false;
      test = 0;

      cartesian = true;
      curvilinear = false;
      cut_cell = false;

      WAF = false;
      MUSCL = false;
      SLIC = true;
      WENO = false;

      xMin = 0.f;
      xMax = 1.f;

      cellCountX = 100;

      deltaX = 1.0/200;    

      end_time = 0.0;
      start_output = false;
      end_output = false;

      FL(0);
      uL(0);
      SL = 0;

      FR(0);
      uR(0);
      SR = 0;

      iface = 0;

      limiter = "vanleer";
    }

    int readConfigFile(const char* filepath){
      printf("Performing file read\n");
      std::cout << filepath << std::endl;
      Config cfg;
      try
      {
        cfg.readFile(filepath);
      }
      catch(const FileIOException &fioex)
      {
        std::cerr << "I/O error while reading file." << std::endl;
        return(EXIT_FAILURE);
      }
      catch(const ParseException &pex)
      {
        std::cerr << "Parse error at " << pex.getFile() << ":" << 
          pex.getLine() << " - " << pex.getError() << std::endl;
        return(EXIT_FAILURE);
      }
      try
      {
        Setting& root = cfg.getRoot();

        Setting& simulation = root["simulation"];

        // get domain information 
        cellCountX = simulation["domain"]["cells"]["x"];

        xMin = simulation["domain"]["dimensions"]["x"][0];
        xMax = simulation["domain"]["dimensions"]["x"][1];

        //get boundary information

        if(simulation.exists("space")){
          Setting& space = simulation["space"];
          if(std::string(space.c_str()) == "cartesian"){
            cartesian = true;
          }
          if(std::string(space.c_str()) == "curvilinear"){
            curvilinear = true;
            cartesian = false;
          }
          if(std::string(space.c_str()) == "cut-cell"){
            cut_cell = true;
            cartesian = false;
          }
        }

        if(simulation.exists("time")){
          time = true;
          end_time = simulation["time"]["end"];
        }

        if(simulation.exists("output")){
          Setting& output = simulation["output"];
          filePath = output["directory"].c_str();
          fileName = output["name"].c_str();
          frequency = output["frequency"];
          start_output = output["start_output"];
          end_output = output["end_output"];
        }

        if(simulation.exists("method")){
          Setting& method = simulation["method"];
          if(std::string(method.c_str()) == "WAF"){
            WAF = true;
            MUSCL = false;
            SLIC = false;
            WENO = false;
            FORCE = false;
          }
          else if(std::string(method.c_str()) == "MUSCL"){
            WAF = false;
            MUSCL = true;
            SLIC = false;
            WENO = false;
            FORCE = false;
          }
          else if(std::string(method.c_str()) == "SLIC"){
            WAF = false;
            MUSCL = false;
            SLIC = true;
            WENO = false; 
            FORCE = false;
          }
          else if(std::string(method.c_str()) == "WENO"){
            WAF = false;
            MUSCL = false;
            SLIC = false;
            WENO = true; 
            FORCE = false;
          }
          else if(std::string(method.c_str()) == "FORCE"){
            WAF = false;
            MUSCL = false;
            SLIC = false;
            WENO = false;
            FORCE = true;
          }
        }

        // get elastic states: (this could be used as well)
        // delegate this to a function?
        if(simulation.exists("riemannstate"))
        {
          Setting& settingState = cfg.lookup("simulation.riemannstate.left");
          // get deformation gradient:
          //
          // can we get the left Riemann state?
          Setting& riemannStateL = simulation["riemannstate"]["left"];
          Setting& riemannStateR = simulation["riemannstate"]["right"];
          unsigned int max = 3, count = 0;
          for (unsigned int i = 0; i < max; i++) {
            for (unsigned int j = 0; j < max; j++) {
              count = i*max + j;
              FL(i,j) = riemannStateL["F"][count];
              FR(i,j) = riemannStateR["F"][count];
              double val = settingState["F"][count]; //alternate way
              /* printf("F(%d,%d) = %4.2f\n", i, j, val); */
            }
          }
          // get velocity components:
          for(int i = 0; i < 3; i++)
          {
            uL(i) = riemannStateL["u"][i];
            uR(i) = riemannStateR["u"][i];
          }

          // get entropy components:
          SR = riemannStateR["S"];
          SL = riemannStateL["S"];
          // in future it will be better to adopt cns approach. 
          // get material: does the material exist?? 
          matL = riemannStateL["mat"].c_str();
          matR = riemannStateR["mat"].c_str();
          // do not create prim state here; wait.
          // get e.o.s properties?
          // get the interface
          iface = simulation["riemannstate"]["discont"];
        } 

       if(simulation.exists("CFL")){
        input_CFL = simulation["CFL"];
       }

       if(simulation.exists("Limiter"))
        limiter = simulation["Limiter"].c_str();

        std::cout << "input_CFL " << input_CFL << std::endl;
        std::cout << "limiter " << limiter << std::endl;

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
        std::cerr << "No 'name' setting in configuration file." << std::endl;
      }

      //we can initialise the domain by calling a function in Elastic1D.cpp
      //and then pass this object to init all states etc and domain throwing
      //errors if they haven't been supplied here.
      //need a DEBUG macro here to check all inputs and output to screen 
      return 0;
    }


    // MPI
    // extend to 2D, WENO etc 
    // Julia, CUDA
    //
   

};
