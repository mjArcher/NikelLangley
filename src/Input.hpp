#include <libconfig.h++>

using namespace libconfig;

//read in and store simulation parameters from cfg file
class Input{
  public:
    int device;

    bool toro;
    int test;

    bool cartesian;
    bool curvilinear; 
    bool cut_cell;

    bool TVDWAF;
    bool MUSCL;

    bool wedge;
    bool cylinder;
    bool NACA;
    bool bubble;
    bool quadrant;
    bool shock;
    bool ambient;

    bool start_output;
    bool end_output;

    bool time;

    double interval;
    int frequency;

    double angle; 
    double shockDiaph_x;
    double shockDiaph_y;
    double shockNorm_x;
    double shockNorm_y;

    double xMin;
    double xMax;
    double yMin;
    double yMax;
    double zMin;
    double zMax;

    double deltaX;
    double deltaY;
    double deltaZ;

    int cellCountX;
    int cellCountY;
    int cellCountZ;
    double end_time;

    double ambient_rho;
    double ambient_u;
    double ambient_v;
    double ambient_p;

    double bubble_rho;
    double bubble_u;
    double bubble_v;
    double bubble_p;
    double bubble_x;
    double bubble_y;
    double radius;

    double rhoL;
    double uL;
    double vL;
    double pL;

    double rhoR;
    double uR;
    double vR;
    double pR;

    double mach;
    double rho_ahead;
    double u_ahead;
    double v_ahead;
    double p_ahead;

    double wedgeStart;
    double wedgeEnd;
    double wedgeAngle;


    std::string filePath;
    std::string fileName;
    std::string geometryPath;

    Input(){
      device = 0;

      time = false;

      toro = false;
      test = 0;
      angle = 0.0;

      cartesian = true;
      curvilinear = false;
      cut_cell = false;

      TVDWAF = false;
      MUSCL = false;

      wedge = false;
      cylinder = false;
      NACA = false;
      bubble = false;
      quadrant = false;
      shock = false;
      ambient = false;

      xMin = 0.f;
      xMax = 1.f;
      yMin = 0.f;
      yMax = 1.f;

      cellCountX = 100;
      cellCountY = 100;

      deltaX = 1.0/200;    
      deltaY = 1.0/200;    

      end_time = 0.0;
      start_output = false;
      end_output = false;

      ambient_rho = 0.0;
      ambient_u = 0.0;
      ambient_v = 0.0;
      ambient_p = 0.0;

      bubble_rho = 0.0;
      bubble_u = 0.0;
      bubble_v = 0.0;
      bubble_p = 0.0;

      bubble_x = 1.0;
      bubble_y = 1.0;
      radius = 0.5;

      rhoL = 0.0;
      uL = 0.0;
      vL = 0.0;
      pL = 0.0;

      rhoR = 0.0;
      uR = 0.0;
      vR = 0.0;
      pR = 0.0;
    }

    void readConfigFile(const char* filepath){

      Config config;

      config.readFile(filepath);

      Setting& root = config.getRoot(); 

      Setting& simulation = root["simulation"];

      cellCountX = simulation["domain"]["cells"]["x"];
      cellCountY = simulation["domain"]["cells"]["y"];
      cellCountZ = simulation["domain"]["cells"]["z"];

      xMin = simulation["domain"]["dimensions"]["x"][0];
      xMax = simulation["domain"]["dimensions"]["x"][1];

      yMin = simulation["domain"]["dimensions"]["y"][0];
      yMax = simulation["domain"]["dimensions"]["y"][1];

      zMin = simulation["domain"]["dimensions"]["z"][0];
      zMax = simulation["domain"]["dimensions"]["z"][1];

      deltaX = (xMax - xMin)/cellCountX;
      deltaY = (yMax - yMin)/cellCountY;
      deltaZ = (zMax - zMin)/cellCountZ;

      if(simulation.exists("device")){
        device = simulation["device"];
      }

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

      if(simulation.exists("toro")){
        toro = true;
        test = simulation["toro"]["test"];
        angle = simulation["toro"]["angle"];
      }

      if(simulation.exists("shock")){
        shock = true;
        mach = simulation["shock"]["mach"];
        rho_ahead = simulation["shock"]["rho_ahead"];	
        u_ahead = simulation["shock"]["u_ahead"];	
        v_ahead = simulation["shock"]["v_ahead"];	
        p_ahead = simulation["shock"]["p_ahead"];	
        shockDiaph_x = simulation["shock"]["diaph"]["x"];	
        shockDiaph_y = simulation["shock"]["diaph"]["y"];	
        shockNorm_x = simulation["shock"]["norm"]["x"];	
        shockNorm_y = simulation["shock"]["norm"]["y"];	
        angle = simulation["shock"]["angle"];
      }

      if(simulation.exists("wedge")){
        wedge = true;
        wedgeStart = simulation["wedge"]["start"];	
        wedgeEnd = simulation["wedge"]["end"];	
        wedgeAngle = simulation["wedge"]["angle"];	
      }

      if(simulation.exists("bubble")){
        bubble = true;
        bubble_rho = simulation["bubble"]["rho"];	
        bubble_u = simulation["bubble"]["u"];	
        bubble_v = simulation["bubble"]["v"];	
        bubble_p = simulation["bubble"]["p"];	
        bubble_x = simulation["bubble"]["centre"]["x"];
        bubble_y = simulation["bubble"]["centre"]["y"];
        radius = simulation["bubble"]["radius"];	
      }

      if(simulation.exists("ambient")){
        ambient = true;
        ambient_rho = simulation["ambient"]["rho"];	
        ambient_u = simulation["ambient"]["u"];	
        ambient_v = simulation["ambient"]["v"];	
        ambient_p = simulation["ambient"]["p"];	
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
        if(std::string(method.c_str()) == "TVDWAF"){
          TVDWAF = true;
          MUSCL = false;
        }
        if(std::string(method.c_str()) == "MUSCL"){
          MUSCL = true;
          TVDWAF = false;
        }
      }

      if(cut_cell){
        if(simulation.exists("geometry")){
          Setting& geometry = simulation["geometry"];
          geometryPath = geometry["directory"].c_str();
        }
      }

    }
};
