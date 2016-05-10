#include "SymmetricMatrix.h"
#include "SquareMatrix.h"

const SymmetricMatrix SymmetricMatrix::Identity(1,0);
const int SymmetricMatrix::Lookup_[3][3]={{0, 3, 5}, {3, 1, 4}, {5, 4, 2}};

SymmetricMatrix Multiply_AT_A(const SquareMatrix& A)
{
    SymmetricMatrix result;
    for (int i=0; i<3; ++i)
    {
        for (int j=i; j<3; ++j)
        {
            result(i,j) = 0.;
            for (int k=0; k<3; ++k)
            {
                result(i,j) += A(k,i) * A(k,j);
            }
        }
    }
    return result;
}

SymmetricMatrix Multiply_A_AT(const SquareMatrix& A)
{
    SymmetricMatrix result;
    for (int i=0; i<3; ++i)
    {
        for (int j=i; j<3; ++j)
        {
            result(i,j) = 0.;
            for (int k=0; k<3; ++k)
            {
                result(i,j) += A(i,k) * A(j,k);
            }
        }
    }
    return result;
}

SymmetricMatrix Multiply_AT_S_A(const SquareMatrix& A, const SymmetricMatrix& S)
{
    SymmetricMatrix result;
    for (int i=0; i<3; ++i)
    {
        for (int j=i; j<3; ++j)
        {
            result(i,j) = 0.;
            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    result(i,j) += A(k,i) * S(k,l) * A(l,j);
                }
            }
        }
    }
    return result;
}

SymmetricMatrix Multiply_A_S_AT(const SquareMatrix& A, const SymmetricMatrix& S)
{
    SymmetricMatrix result;
    for (int i=0; i<3; ++i)
    {
        for (int j=i; j<3; ++j)
        {
            result(i,j) = 0;
            for (int k=0; k<3; ++k)
            {
                for (int l=0; l<3; ++l)
                {
                    result(i,j) += A(i,k) * S(k,l) * A(j,l);
                }
            }
        }
    }
    return result;
}

TinyVector<double, 3> operator* (const SymmetricMatrix& A, const TinyVector<double, 3>& v) 
{
  return TinyVector<double, 3>(A(0,0)*v[0]+A(0,1)*v[1]+A(0,2)*v[2],
			       A(1,0)*v[0]+A(1,1)*v[1]+A(1,2)*v[2],
			       A(2,0)*v[0]+A(2,1)*v[1]+A(2,2)*v[2]);
}  

SymmetricMatrix SymmetricMatrix::Inverse ( void ) const
{      
  const SymmetricMatrix detInv = DetInverse();
  
  // Determinant
  const double det = Data_[0]*detInv(0,0)+Data_[3]*detInv(0,1)+Data_[5]*detInv(0,2);
  
  // Inverse
  return detInv/det;
}


/** Stream output of symmetric matrix */
ostream& operator<< ( ostream& os, const SymmetricMatrix& M )
{
    for ( int i=0; i<3; ++i )
    {
        os << setprecision(10);
        os << "\t[ ";
        for ( int j=0; j<3; ++j )
        {
            os << setw (14) << M( i,j ) << " ";
        }
        os << "]\n";
    }
    return os;
}

SymmetricMatrix operator*(const double s, const SymmetricMatrix& M) {
  SymmetricMatrix C=M;
  return C*=s;
}

SymmetricMatrix operator*(const SymmetricMatrix& M, const double s) {
  SymmetricMatrix C=M;
  return C*=s;
}

SymmetricMatrix operator/(const SymmetricMatrix& M, const double s) {
  SymmetricMatrix C=M;
  return C*=(1./s);
}

SquareMatrix operator*(const SymmetricMatrix& A, const SymmetricMatrix& B) {
  SquareMatrix C;
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      C(i,j) = 0.;
      for(int k=0; k<3; ++k) {
	C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return C;
}

SquareMatrix operator*(const SquareMatrix& A, const SymmetricMatrix& B) {
  SquareMatrix C;
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      C(i,j) = 0.;
      for(int k=0; k<3; ++k) {
	C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return C;
}

SquareMatrix operator*(const SymmetricMatrix& A, const SquareMatrix& B) {
  SquareMatrix C;
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      C(i,j) = 0.;
      for(int k=0; k<3; ++k) {
	C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return C;
}

SymmetricMatrix operator-(const SymmetricMatrix& A, const SymmetricMatrix& B) {
  SymmetricMatrix C=A;
  return C-=B;
}

SquareMatrix operator-(const SymmetricMatrix& A, const SquareMatrix& B) {
  SquareMatrix C(A);
  return C-=B;
}

SquareMatrix operator-(const SquareMatrix& A, const SymmetricMatrix& B) {
  SquareMatrix C(A);
  return C-=SquareMatrix(B);
}

SymmetricMatrix operator+(const SymmetricMatrix& A, const SymmetricMatrix& B) {
  SymmetricMatrix C=A;
  return C+=B;
}

SymmetricMatrix CommutativeMult(const SymmetricMatrix& A, const SymmetricMatrix& B) {
  SymmetricMatrix C;
  for(int i=0; i<3; ++i) {
    for(int j=i; j<3; ++j) {
      C(i,j) = 0.;
      for(int k=0; k<3; ++k) {
				C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return C;
}


void EigenDecomposeRotate(SymmetricMatrix& A,
                          const int i, const int j, const int k, const int l,
                          const double s, const double tau) {
  const double g=A(i,j);
  const double h=A(k,l);
  A(i,j) = g - s * (h + g * tau);
  A(k,l) = h + s * (g - h * tau);
}

void EigenDecomposeRotate(SquareMatrix& A,
                          const int i, const int j, const int k, const int l,
                          const double s, const double tau) {
  const double g=A(i,j);
  const double h=A(k,l);
  A(i,j) = g - s * (h + g * tau);
  A(k,l) = h + s * (g - h * tau);
}

blitz::TinyVector<double, 3> EigenDecompose(const SymmetricMatrix &A1,
                                     SquareMatrix& EigenVecs) {
  blitz::TinyVector<double, 3> b, z, EigenVals;
  EigenVecs.SetIdentity();
  SquareMatrix A = SquareMatrix(A1);
  
  const double r32 = 0.2/9.0;

  /* Initialize b and EigenVals to the diagonal of A */
  for (unsigned i=0; i<3; ++i) {
    EigenVals[i] = A(i,i);
  }

  b=EigenVals;
  z = 0.0;

  for (unsigned iter=0; iter<50; ++iter) {
    /* Sum off-diagonal elements */
    const double sm=fabs(A(0,1)) + fabs(A(0,2)) + fabs(A(1,2));

    if (sm == 0.0) return EigenVals;

    const double tresh = (iter<4 ? r32*sm : 0.0);

    for (unsigned i=0; i<3; ++i) {
      for (unsigned j=i+1; j<3; ++j) {
        double& Aij   = A(i,j);
        const double absAij = fabs(Aij);
        const double absdi  = fabs(EigenVals[i]);
        const double absdj  = fabs(EigenVals[j]);

        const double g=100.0*absAij;
	
        /* After four sweeps, skip the rotation if the off-diagonal element is small. */
        if (iter > 4 && double(absdi+g) == absdi && double(absdj+g) == absdj) {
          Aij=0.0;
        } else if (absAij > tresh) {
          const double h = EigenVals[j]-EigenVals[i];
          const double absh = fabs(h);
          double t;

          if (double(absh+g) == absh)
            t = Aij / h;
          else {
            const double theta = 0.5 * h / Aij;
            t = 1.0 / (fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }

          const double c = 1.0 / sqrt(1.0+t*t);
          const double s = t*c;
          const double tau = s / (1.0+c);
          const double hh = t*Aij;
          z[i] -= hh;
          z[j] += hh;
          EigenVals[i] -= hh;
          EigenVals[j] += hh;
          Aij = 0.0;

          /* Case of rotations 0 ≤ j < p. */
          for (unsigned k=0;k<i;++k) {
            EigenDecomposeRotate(A,k,i,k,j, s , tau);
          }

          /* Case of rotations p < j < q. */
          for (unsigned k=i+1;k<j;++k) {
            EigenDecomposeRotate(A,i,k,k,j, s, tau);
          }

          /* Case of rotations q < j ≤ n. */
          for (unsigned k=j+1;k<3;++k) {
            EigenDecomposeRotate(A,i,k,j,k, s, tau);
          }

          for (unsigned k=0;k<3;++k) {
            EigenDecomposeRotate(EigenVecs,k,i,k,j, s, tau);
          }
        }
      }
    }

    b += z;
    EigenVals = b;
    z = 0.;

  }

  cerr << "Too many iterations in EigenDecompose(SymmetricMatrix)\n";
	cerr << "Matrix = " << A1 << endl;
  exit(1);

  return EigenVals;
}

double dot(const SymmetricMatrix& A, const SymmetricMatrix& B) {
  double tot=0.;
  for(int i=0; i<3; ++i) {
    for(int j=0; j<3; ++j) {
      tot+=A(i,j)*B(i,j);
    }
  }
  return tot;
}
