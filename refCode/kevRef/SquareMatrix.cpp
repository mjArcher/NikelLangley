#include "SquareMatrix.h"

/** Identity matrix */
const SquareMatrix SquareMatrix::Identity(1,0);

double dot(const SquareMatrix& A, const SquareMatrix& B) {
  double tot=0;
  for ( int i=0; i<3; ++i )
    {
      for ( int j=0; j<3; ++j )
        {
	  tot+=A(i,j)*B(i,j);
        }
    }
  return tot;
}

SquareMatrix operator* ( const double scal, const SquareMatrix& M )
{
    SquareMatrix C = M;
    return C*=scal;
}

SquareMatrix operator* ( const SquareMatrix& M, const double scal )
{
    SquareMatrix C = M;
    return C*=scal;
}

SquareMatrix operator/ ( const SquareMatrix& M, const double scal )
{
    SquareMatrix C = M;
    return C/=scal;
}

ostream& operator<< ( ostream& os, const SquareMatrix& M )
{
    os << setprecision(10);
    for ( int i=0; i<3; ++i )
    {
        os << "\t[ ";
        for ( int j=0; j<3; ++j )
        {
            os << setw(16) << M(i,j) << " ";
        }
        os << "]\n";
    }
    return os;
}

SquareMatrix operator+ (const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix C=A;
    return C+=B;
}

SquareMatrix operator- (const SquareMatrix& A, const SquareMatrix& B) {
    SquareMatrix C=A;
    return C-=B;
}

SquareMatrix rotateToX(const TinyVector<double, 2>& reference)
{  
  const double cos_theta = reference[0];  
  const double sin_theta = reference[1];
  return SquareMatrix(cos_theta, sin_theta, 0.0, -sin_theta, cos_theta, 0.0, 0.0, 0.0, 1.0);
}

TinyVector<double, 3> operator* (const SquareMatrix& A, const TinyVector<double, 3>& v) 
{
  return TinyVector<double, 3>(A(0,0)*v[0]+A(0,1)*v[1]+A(0,2)*v[2],
			       A(1,0)*v[0]+A(1,1)*v[1]+A(1,2)*v[2],
			       A(2,0)*v[0]+A(2,1)*v[1]+A(2,2)*v[2]);
}  

bool operator==(const SquareMatrix& M1, const SquareMatrix& M2) 
{
  for(int i=0; i<3; ++i) 
  {
    for(int j=0; j<3; ++j) 
    {
      if(M1(i,j)!=M2(i,j)) return false;
    }
  }
  return true;
}

bool operator>(const SquareMatrix& M1, const SquareMatrix& M2) 
{
  for(int i=0; i<3; ++i) 
  {
    for(int j=0; j<3; ++j) 
    {
      if(M1(i,j)>M2(i,j)) return true;
    }
  }
  return false;
}

bool operator<(const SquareMatrix& M1, const SquareMatrix& M2) 
{
  for(int i=0; i<3; ++i) 
  {
    for(int j=0; j<3; ++j) 
    {
      if(M1(i,j)<M2(i,j)) return true;
    }
  }
  return false;
}

bool operator!=(const SquareMatrix& M1, const SquareMatrix& M2)
{
  return !(M1==M2);
}

SquareMatrix abs(const SquareMatrix& M) 
{
  SquareMatrix A;
  for(int i=0; i<3; ++i) 
  {
    for(int j=0; j<3; ++j) 
    {
      A(i,j) = fabs(M(i,j));
    }
  }
  return A;
}
