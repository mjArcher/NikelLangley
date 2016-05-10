#include "Domain.hpp"
#include "stdio.h"
using namespace std;


//problem with this approach is that we need to pass this class everywhere we need to update the domain
//pass by reference or pass by pointer?

//operators
//functions to get domain information, i.e. dx, rows, cols

Domain::Domain(int rows, int cols):rows(rows),cols(cols){ 
  vector<vector<double> > A(rows, vector<double>(cols));
  dom = A;
}

int main(int argc, char ** argv)
{
  // this all works
  Domain dom(100, 50);
  cout << dom.getRows() << endl;
  printf("\nGet element from the array\n");
  cout << dom.returnDomain()[20][30] << endl;


}
