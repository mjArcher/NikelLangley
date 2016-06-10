#include "Domain.hpp"
#include "stdio.h"
#include <sstream>
#include "ElasticPrimState.h" 
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
////lx is in the j direction
////ly is in the i direction

template <class T> Domain<T>::Domain(int Ni, int Nj, int GC, double Lx):Ni(Ni), Nj(Nj), GC(GC), Lx(Lx) { 
  vector<vector<T> > A(Ni + 2*GC, vector<T>(Nj + 2*GC));
  starti = GC, startj = GC, endi = Ni+GC, endj = Nj+GC;
  GNi = Ni+2*GC, GNj = Nj+2*GC;
  dom = A;
  Ly = Lx * (double)Ni/Nj;
  dx = Lx/Nj;
  //
}

template <class T> Domain<T>::~Domain()
{
}

template <class T>
vector<T> Domain<T>::getRowi(int row) const
{
  vector<T> rowv(Nj);
  for(int j = startj; j < endj; j++)
  {
    rowv[j-GC] = dom[row][j];
  }
  return rowv;
}

template <class T>
vector<T> Domain<T>::getColj(int col) const
{
  vector<T> colv(Ni);
  for(int i = starti; i < endi; i++)
  {
    colv[i-GC] = dom[i][col];
  }
  return colv;
}

//outputs the entire row plus boundary cells
template <class T>
vector<T> Domain<T>::getGRowi(int row) const
{
  vector<T> rowv(GNj);
  for(int j = 0; j < GNj; j++)
  {
    rowv[j] = dom[row][j];
  }
  return rowv;
}

//entire col plus boundary cells
template <class T>
vector<T> Domain<T>::getGColj(int col) const
{
  vector<T> colv(GNi);
  for(int i = 0; i < GNi; i++)
  {
    colv[i] = dom[i][col];
  }
  return colv;
}


/* template <class T> */
/* void Domain<T>::assignRowi(unsigned int rowi, vector<T> row) */ 
/* { */
/*   for(int j = 0; j < GNj; j++) */
/*   { */
/*     dom[rowi][j] = row[j]; */
/*   } */
/* } */

/* template <class T> */
/* void Domain<T>::assignColj(unsigned int colj, vector<T> col) */ 
/* { */
/*   for(int i = 0; i < GNi; i++) */
/*   { */
/*     dom[i][colj] = col[i]; */
/*   } */
/* } */

template <class T>
T Domain<T>::operator()(unsigned int row, unsigned int col) const
{
  return dom[row][col];
}


//get the cartesian position indexed by row and col
//think about what output this should be i.e. a map? tuple?
template <class T>
T Domain<T>::getCellCartPos(unsigned int row, unsigned int col)
{
  
}

template <class T>
std::ostream& operator<<(std::ostream& os, const vector<T>& param) 
{
	os.precision(3);
  cout << "Hello world" << endl;
  for (int i = 0; i < param.size(); i++)
  {
    os << param[i] << endl;
  }
  return os;
}

/* int main(int argc, char ** argv) */
/* { */
/*   Domain<ElasticPrimState> dom(10, 40, 2, 1.0); */
/*   cout << dom.getNi() << endl; */
/*   printf("\nGet element from the array\n"); */
/*   cout << dom.returnDomain()[5][20] << endl; */
/*   cout << "output col to screen" << endl; */
/*   //use more sophisticated getter */
/* } */

