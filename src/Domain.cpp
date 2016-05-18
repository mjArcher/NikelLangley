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
//

Domain::Domain2D(int rows, int cols):rows(rows),cols(cols){ 
  vector<vector<double> > A(rows, vector<double>(cols));
  dom = A;
}

double Domain2D::operator() (unsigned row, unsigned col) const
{
  return dom[row][col];
}

/* int main(int argc, char ** argv) */
/* { */
/*   // this all works */
/*   Domain dom(100, 50); */
/*   cout << dom.getRows() << endl; */
/*   printf("\nGet element from the array\n"); */
/*   cout << dom.returnDomain()[20][30] << endl; */

/*   //concatenate string and double with iostreams */
/*   stringstream ss; */

/*   cout << "output col to screen (eventually implement function"  << endl; */
/*   for(int i = 0; i < 10; ++i) */
/*   { */
/*     /1* iost << dom(i,0) *1/ */
/*     /1* //then flush *1/ */
/*     cout << i << "\t" << dom(i,0) << endl; */
/*   } */

}
