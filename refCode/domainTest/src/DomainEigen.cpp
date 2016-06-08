#include "DomainEigen.h"

using namespace std;
using namespace Eigen;

//initialise the 2D matrix 
void initialise(Domain& dom)
{

  Vector3d vec(2.,2.,2.);  
  cout << vec << endl;
  Matrix3d F;
  F << 3,0,0,0,3,0,0,0,3;
  double S(1.);

  //first state
  ElasticPrimState prim1(vec,F,S);

  //second state
  ElasticPrimState prim2;
  prim2 = prim1*2;

  //perform initialisation using a simple block function
  std::cout << dom.rows() << endl;
  /* for(int i = 0; i < dom.rows()/2.; i++) */
  /* { */
  /* dom.block(49,49,50,50) << prim1, prim1; */
  for (int i = 0; i < dom.rows(); i++)
  {
    for(int j = 0; j < dom.cols(); j++)
    {
      if(i < dom.rows()/2)
        dom(i,j) = prim1;
      else
        dom(i,j) = prim2;
    }
  }
 
  cout << "position 1" << endl;
  std::cout << dom(0,10) << std::endl;
  cout << "position 2" << endl;
  std::cout << dom(60,10) << std::endl;
}

int eigenMatrix(int rows, int cols)
{
  Domain dom(rows, cols);
  initialise(dom);
  Domain domT = dom.transpose();
}



int main(int argc, char** argv)
{
  std::cout << "Hello world " << std::endl;
  eigenMatrix(100,100);
}
