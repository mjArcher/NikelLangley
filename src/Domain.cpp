#include "Domain.hpp"
#include "stdio.h"
#include <sstream>
using namespace std;


//problem with this approach is that we need to pass this class everywhere we need to update the domain
//pass by reference or pass by pointer?

//operators
//functions to get domain information, i.e. dx, rows, cols
//wrap domain in material?
//
//Requirements:
//
//Want to be able to create any problem type, multiple materials etc
//implement mask.
//currently have a Material class which has vectors of ElasticStates
//could eventually have a material class which has also a mask - boolean array (indicates where on the domain the material is)
//might be quite tricky: since we have states wrapped up in material structs.
//we can get the y-length from the ratio between the number of x and y cells - always assume square cells!
Domain::Domain(int Ni, int Nj, int GC, double Lx):Ni(Ni), Nj(Nj), GC(GC), Lx(Lx) { 
  vector<vector<double> > A(Ni + 2*GC, vector<double>(Nj + 2*GC));
  starti = GC, startj = GC, endi = Ni+GC, endj = Nj+GC;
  GNi = Ni+2*GC, GNj = Nj+2*GC;
  dom = A;
}

Domain::~Domain()
{
}

vector<double> Domain::getRowi(int row) const
{
  vector<double> rowv(Nj);
  for(int j = startj; j < endj; j++)
  {
    rowv[j-GC] = dom[row][j];
  }
  return rowv;
}

vector<double> Domain::getColj(int col) const
{
  vector<double> colv(Ni);
  for(int i = starti; i < endi; i++)
  {
    colv[i-GC] = dom[i][col];
  }
  return colv;
}

//outputs the entire row plus boundary cells
vector<double> Domain::getGRowi(int row) const
{
  vector<double> rowv(GNj);
  for(int j = 0; j < GNj; j++)
  {
    rowv[j] = dom[row][j];
  }
  return rowv;
}

//entire col plus boundary cells
vector<double> Domain::getGColj(int col) const
{
  vector<double> colv(GNi);
  for(int i = 0; i < GNi; i++)
  {
    colv[i] = dom[i][col];
  }
  return colv;
}

double Domain::operator() (unsigned row, unsigned col) const
{
  return dom[row][col];
}

std::ostream& operator<<(std::ostream& os, const vector<double>& param) 
{
	os.precision(3);
  cout << "Hello world" << endl;
  for (int i = 0; i < param.size(); i++)
  {
    os << param[i] << endl;
  }
  return os;
}

int main(int argc, char ** argv)
{
  //this all works
  Domain dom(10, 30, 2, 1.0);
  cout << dom.getNi() << endl;
  printf("\nGet element from the array\n");
  cout << dom.returnDomain()[5][20] << endl;

  //concatenate string and double with iostreams
  stringstream ss;
  cout << "output col to screen" << endl;
  cout << dom.getColj(0) << endl;
}

