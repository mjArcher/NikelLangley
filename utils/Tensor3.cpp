#include "Tensor3.h"
#include <iomanip>

const int Symmetric23Tensor3::Lookup_[3][3]={{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};

Tensor3 operator-(const Tensor3& T1, const Tensor3& T2) 
{
  Tensor3 D(T1);
  D-=T2;
  return D;
}


ostream& operator<<(ostream& os, const Tensor3& T) 
{
  for(int j=0; j<3; ++j) 
  {
    os << "{ ";
    for(int i=0; i<3; ++i) 
    {
      os << "{ ";
      for(int k=0; k<3; ++k) 
      {
	os << setw(9) << setprecision(2) << scientific << T(i,j,k) << " ";
      }
      os << "} ";
    }
    os << "}\n";
  }
  return os;
}
