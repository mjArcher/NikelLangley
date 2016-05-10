#ifndef __SQUAREMATRIX_H
#define __SQUAREMATRIX_H

#include "SymmetricMatrix.h"
#include "blitz/tinyvec.h"
#include "blitz/tinyvec-et.h"
using namespace blitz;

class SquareMatrix
{

    double Data_[3][3];

public:

    /** Identity matrix */
    static const SquareMatrix Identity;

    SquareMatrix() {}
    explicit SquareMatrix ( const double scal )  {
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                Data_[i][j] = scal;
            }
        }
    }

    SquareMatrix ( const double diag, const double off_diag ) {
        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                if (i==j)
                    Data_[i][j] = diag;
                else
                    Data_[i][j] = off_diag;
            }
        }
    }

    explicit SquareMatrix (const double * const data)  {
        const double * datap = data;
        for (int j=0; j<3; ++j) {
            for (int i=0; i<3; ++i, ++datap) {
                Data_[i][j] = *datap;
            }
        }
    }

    SquareMatrix ( const double A,const double B, const double C, const double D, const double E, const double F, const double G, const double H, const double I) {
        Data_[0][0]=A;
        Data_[0][1]=B;
        Data_[0][2]=C;
        Data_[1][0]=D;
        Data_[1][1]=E;
        Data_[1][2]=F;
        Data_[2][0]=G;
        Data_[2][1]=H;
        Data_[2][2]=I;
    }

  TinyVector<double,3> col(const int j) const
    {
      TinyVector<double, 3> c;
      for(int i=0; i<3; ++i) 
      {
	c[i] = (*this)(i,j);
      }
      return c;
    }

    TinyVector<double,3> row(const int i) const
    {
      TinyVector<double, 3> r;
      for(int j=0; j<3; ++j) 
      {
	r[j] = (*this)(i,j);
      }
      return r;
    }

    SquareMatrix ( const SquareMatrix& B ) {
      for (int j=0; j<3; ++j)
            for (int i=0; i<3; ++i)
                Data_[i][j] = B.Data_[i][j];
    }

    explicit SquareMatrix ( const SymmetricMatrix& B ) {
        for (int i=0; i<3; ++i) 
            for (int j=0; j<3; ++j) 
                Data_[i][j] = B(i,j);
    }

    SquareMatrix& operator*= ( const double s )  {
        for ( int i=0; i<3; ++i )
            for ( int j=0; j<3; ++j )
                Data_[i][j]*=s;
        return *this;
    }

    inline double& operator() (const int i, const int j) {
      assert(i>=0 && i<3 && j>=0 && j<3);
      return Data_[i][j];
    }
    
    inline const double& operator() (const int i, const int j) const {
      assert(i>=0 && i<3 && j>=0 && j<3);
      return Data_[i][j];
    }
    
    SquareMatrix& operator/= ( const double s ) {
        return (*this)*=1./s;
    }

    SquareMatrix& SetZero(void) {
        for ( int i=0; i<3; ++i )
            for ( int j=0; j<3; ++j )
                Data_[i][j]=0.;
        return *this;
    }
      
    double Invariant2(void) const {
      return Data_[0][0]*Data_[1][1]+Data_[1][1]*Data_[2][2]+Data_[0][0]*Data_[2][2]-Data_[0][1]*Data_[0][1]-Data_[1][2]*Data_[1][2]-Data_[0][2]*Data_[0][2];      
    }
    
    SquareMatrix& operator+=( const SquareMatrix& B ) {
        for ( int i=0; i<3; ++i )
            for ( int j=0; j<3; ++j )
                Data_[i][j]+=B.Data_[i][j];
        return *this;
    }

    SquareMatrix& operator-=( const SquareMatrix& B ) {
        for ( int i=0; i<3; ++i )
            for ( int j=0; j<3; ++j )
                Data_[i][j]-=B.Data_[i][j];
        return *this;
    }

    SquareMatrix operator* ( const SquareMatrix& B ) const {
        SquareMatrix container ( 0. );
        for ( int i=0; i<3; ++i )
        {
            for ( int j=0; j<3; ++j )
            {
                for ( int k=0; k<3; ++k )
                {
                    container.Data_[i][j] += Data_[i][k] * B.Data_[k][j];
                }
            }
        }
        return container;
    }

    SquareMatrix& SetIdentity ( void )
    {
        *this = Identity;
        return *this;
    }

    SquareMatrix& SetAll ( const double scal )
    {
        for (int i=0; i<3; ++i) {
            for (int j=0; j<3; ++j) {
                Data_[i][j]=scal;
            }
        }
        return *this;
    }

    SquareMatrix Transpose ( void ) const
    {
        SquareMatrix container;

        for ( int i=0; i<3; ++i )
        {
            for ( int j=0; j<3; ++j )
            {
                container.Data_[j][i] = Data_[i][j];
            }
        }
        return container;
    }

    SquareMatrix Inverse ( void ) const
    {
      // Minors
        const double m11 = Data_[1][1]*Data_[2][2]-Data_[1][2]*Data_[2][1];
        const double m12 = Data_[1][2]*Data_[2][0]-Data_[1][0]*Data_[2][2];
        const double m13 = Data_[1][0]*Data_[2][1]-Data_[1][1]*Data_[2][0];
        const double m21 = Data_[0][2]*Data_[2][1]-Data_[0][1]*Data_[2][2];
        const double m22 = Data_[0][0]*Data_[2][2]-Data_[0][2]*Data_[2][0];
        const double m23 = Data_[0][1]*Data_[2][0]-Data_[0][0]*Data_[2][1];
        const double m31 = Data_[0][1]*Data_[1][2]-Data_[0][2]*Data_[1][1];
        const double m32 = Data_[0][2]*Data_[1][0]-Data_[0][0]*Data_[1][2];
        const double m33 = Data_[0][0]*Data_[1][1]-Data_[0][1]*Data_[1][0];
	
	// Determinant
        const double det = Data_[0][0]*m11+Data_[0][1]*m12+Data_[0][2]*m13;
        const double rdet = 1./det;
	
	// Inverse
        return SquareMatrix(rdet*m11, rdet*m21, rdet*m31, rdet*m12, rdet*m22, rdet*m32, rdet*m13, rdet*m23, rdet*m33);
    }

    SquareMatrix DetInverse ( void ) const
    {
      // Minors
        const double m11 = Data_[1][1]*Data_[2][2]-Data_[1][2]*Data_[2][1];
        const double m12 = Data_[1][2]*Data_[2][0]-Data_[1][0]*Data_[2][2];
        const double m13 = Data_[1][0]*Data_[2][1]-Data_[1][1]*Data_[2][0];
        const double m21 = Data_[0][2]*Data_[2][1]-Data_[0][1]*Data_[2][2];
        const double m22 = Data_[0][0]*Data_[2][2]-Data_[0][2]*Data_[2][0];
        const double m23 = Data_[0][1]*Data_[2][0]-Data_[0][0]*Data_[2][1];
        const double m31 = Data_[0][1]*Data_[1][2]-Data_[0][2]*Data_[1][1];
        const double m32 = Data_[0][2]*Data_[1][0]-Data_[0][0]*Data_[1][2];
        const double m33 = Data_[0][0]*Data_[1][1]-Data_[0][1]*Data_[1][0];
		
	// Inverse
        return SquareMatrix(m11, m21, m31, m12, m22, m32, m13, m23, m33);
    }

    double Trace ( void ) const
    {
        return Data_[0][0]+Data_[1][1]+Data_[2][2];
    }
    
    double Determinant ( void ) const {
      return Data_[0][0] * ( Data_[1][1] * Data_[2][2] - Data_[1][2]*Data_[2][1] )
           - Data_[1][0] * ( Data_[0][1] * Data_[2][2] - Data_[2][1]*Data_[0][2] )
           + Data_[2][0] * ( Data_[0][1] * Data_[1][2] - Data_[1][1]*Data_[0][2] );
}

    double Schur2 (void) const 
    {
      double schur=0.;
      for(int i=0; i<3; ++i) 
      {
	for(int j=0; j<3; ++j) 
	{
	  schur+=Data_[i][j]*Data_[i][j];
	}
      }
      return schur;
    }
    
    double Schur(void) const 
    {
      return sqrt(Schur2());
    }
};

SquareMatrix operator+ (const SquareMatrix& A, const SquareMatrix& B);
SquareMatrix operator- (const SquareMatrix& A, const SquareMatrix& B);
SquareMatrix operator* ( const double scal, const SquareMatrix& M );
SquareMatrix operator* ( const SquareMatrix& M, const double scal );
SquareMatrix operator/ ( const SquareMatrix& M, const double scal );

bool operator==(const SquareMatrix& M1, const SquareMatrix& M2);
bool operator!=(const SquareMatrix& M1, const SquareMatrix& M2);
bool operator<(const SquareMatrix& M1, const SquareMatrix& M2);
bool operator>(const SquareMatrix& M1, const SquareMatrix& M2);
SquareMatrix abs(const SquareMatrix& M);

TinyVector<double, 3> operator* (const SquareMatrix& A, const TinyVector<double, 3>& v);

double dot(const SquareMatrix& A, const SquareMatrix& B);

ostream& operator<< ( ostream& os, const SquareMatrix& M );

SquareMatrix rotateToX(const TinyVector<double, 2>& reference);
SquareMatrix rotateToX(const TinyVector<double, 3>& reference);

#endif
