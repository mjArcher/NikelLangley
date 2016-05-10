#ifndef SYMMETRICMATRIX_H
#define SYMMETRICMATRIX_H

#include <iostream>
#include <blitz/tinyvec.h>
using namespace std;
using namespace blitz;

/** 3D Symmetric matrix stored in Voigt-style vector form */
class SymmetricMatrix {

	private:

		double Data_[6];
		static const int Lookup_[3][3];

	public:

		/** Identity matrix */
		static const SymmetricMatrix Identity;

		/** Default constructor - all components uninitialized */
		SymmetricMatrix ( void ) {}

		/** Constant component constructor - sets all components of vector to specified value
		 *  @param val Value to be assigned to all components */
		explicit SymmetricMatrix ( const double val ) {
			for (int i=0; i<6; ++i) Data_[i] = val;
		}

		/** Diagonal constructor - sets all diagonal components of vector to single value and all off-diagonal components to
		 *  another value.
		 *  @param diag_val Component value on the diagonal
		 *  @param offdiag_val Component value for all non-diagonal entries */
		SymmetricMatrix ( const double diag_val, const double offdiag_val )
		{
			for (int i=0; i<3; ++i)
			{
				Data_[i] = diag_val;
			}
			for (int i=3; i<6; ++i)
			{
				Data_[i] = offdiag_val;
			}
		}

		SymmetricMatrix(const double A, const double B, const double C, const double D, const double E, const double F) {
			Data_[0]=A;
			Data_[1]=B;
			Data_[2]=C;
			Data_[3]=D;
			Data_[4]=E;
			Data_[5]=F;
		}

		/** Copy constructor
		 *  @param B Matrix to copy from */
		SymmetricMatrix ( const SymmetricMatrix& B )  {
			for(int i=0; i<6; ++i) Data_[i]=B.Data_[i];
		}

		SymmetricMatrix& operator=(const SymmetricMatrix& B) {
			for(int i=0; i<6; ++i) Data_[i]=B.Data_[i];
			return *this;
		}

		SymmetricMatrix& SetZero(void) {
			for(int i=0; i<6; ++i) Data_[i]=0.;
			return *this;
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

		const SymmetricMatrix& Transpose(void) const 
		{
			return *this;
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

		/** Component access - non-const version
		 * @param i Matrix row (0 = first row)
		 * @param j Matrix column (0 = first column)
		 * @return Reference to component at i, j */
		inline double& operator() ( const int i, const int j )
		{
#ifdef BOUNDS_CHECK
			if (i<0 || i>2 || j<0 || j>2)
			{
				cerr << "Bad bounds in SymmetricMatrix::operator()\n";
				exit(1);
			}
#endif
			return Data_[Lookup_[i][j]];
		}

		/** Component access - const version
		 * @param i Matrix row (0 = first row)
		 * @param j Matrix column (0 = first column)
		 * @return Value of component at i, j */
		inline const double& operator() ( const int i, const int j ) const
		{
#ifdef BOUNDS_CHECK
			if (i<0 || i>2 || j<0 || j>2)
			{
				cerr << "Bad bounds in SymmetricMatrix::operator()\n";
				exit(1);
			}
#endif
			return Data_[Lookup_[i][j]];
		}

		double Determinant ( void ) const
		{
			return Data_[0] * ( Data_[1]*Data_[2] - Data_[4]*Data_[4] )
				- Data_[3] * ( Data_[3]*Data_[2] - Data_[4]*Data_[5] )
				+ Data_[5] * ( Data_[3]*Data_[4] - Data_[1]*Data_[5] );
		}

		double Trace(void) const {
			return Data_[0] + Data_[1] + Data_[2];
		}

		double Schur2(void) const {
			return Data_[0]*Data_[0]+Data_[1]*Data_[1]+Data_[2]*Data_[2] 
				+ 2.*(Data_[3]*Data_[3]+Data_[4]*Data_[4]+Data_[5]*Data_[5]);
		}

		double Schur(void) const {
			return sqrt(Schur2());
		}

		SymmetricMatrix Deviator(void) const
		{
			const double tr = Trace() / 3.0;
			SymmetricMatrix dev;
			dev.Data_[0] = Data_[0] - tr;
			dev.Data_[1] = Data_[1] - tr;
			dev.Data_[2] = Data_[2] - tr;
			dev.Data_[3] = Data_[3];
			dev.Data_[4] = Data_[4];
			dev.Data_[5] = Data_[5];
			return dev;
		}    


		SymmetricMatrix Inverse ( void ) const;

		SymmetricMatrix DetInverse ( void ) const
		{      
			// Minors
			const double m11 = Data_[1]*Data_[2]-Data_[4]*Data_[4];
			const double m12 = Data_[4]*Data_[5]-Data_[3]*Data_[2];
			const double m13 = Data_[3]*Data_[4]-Data_[1]*Data_[5];
			const double m22 = Data_[0]*Data_[2]-Data_[5]*Data_[5];
			const double m23 = Data_[3]*Data_[5]-Data_[0]*Data_[4];
			const double m33 = Data_[0]*Data_[1]-Data_[3]*Data_[3];

			// Inverse
			return SymmetricMatrix(m11, m22, m33, m12, m23, m13);
		}

		SymmetricMatrix& operator*=(const double s) {
			for(int i=0; i<6; ++i) Data_[i]*=s;
			return *this;
		}

		SymmetricMatrix& operator-=(const SymmetricMatrix& B) {
			for(int i=0; i<6; ++i) Data_[i]-=B.Data_[i];
			return *this;
		}

		SymmetricMatrix& operator+=(const SymmetricMatrix& B) {
			for(int i=0; i<6; ++i) Data_[i]+=B.Data_[i];
			return *this;
		}

};

/** Stream output of symmetric matrix */
ostream& operator<< ( ostream& os, const SymmetricMatrix& M );

class SquareMatrix;
SymmetricMatrix operator*(const double s, const SymmetricMatrix& M);
SymmetricMatrix operator*(const SymmetricMatrix& M, const double s);
SymmetricMatrix operator/(const SymmetricMatrix& M, const double s);
SquareMatrix operator-(const SymmetricMatrix& A, const SquareMatrix& B);
SquareMatrix operator-(const SquareMatrix& A, const SymmetricMatrix& B);
SymmetricMatrix operator-(const SymmetricMatrix& A, const SymmetricMatrix& B);
SymmetricMatrix operator+(const SymmetricMatrix& A, const SymmetricMatrix& B);
SquareMatrix operator*(const SymmetricMatrix& A, const SymmetricMatrix& B);
SquareMatrix operator*(const SquareMatrix& A, const SymmetricMatrix& B);
SquareMatrix operator*(const SymmetricMatrix& A, const SquareMatrix& B);

TinyVector<double, 3> operator* (const SymmetricMatrix& A, const TinyVector<double, 3>& v); 

bool operator==(const SymmetricMatrix& M1, const SymmetricMatrix& M2);
bool operator!=(const SymmetricMatrix& M1, const SymmetricMatrix& M2);
bool operator==(const SquareMatrix& M1, const SymmetricMatrix& M2);
bool operator==(const SymmetricMatrix& M1, const SquareMatrix& M2);
bool operator!=(const SquareMatrix& M1, const SymmetricMatrix& M2);
bool operator!=(const SymmetricMatrix& M1, const SquareMatrix& M2);

double dot(const SymmetricMatrix& A, const SymmetricMatrix& B);

SymmetricMatrix CommutativeMult(const SymmetricMatrix& A, const SymmetricMatrix& B);

SymmetricMatrix Multiply_AT_A(const SquareMatrix& A);
SymmetricMatrix Multiply_A_AT(const SquareMatrix& A);
SymmetricMatrix Multiply_AT_S_A(const SquareMatrix& A, const SymmetricMatrix& S);
SymmetricMatrix Multiply_A_S_AT(const SquareMatrix& A, const SymmetricMatrix& S);

blitz::TinyVector<double, 3> EigenDecompose(const SymmetricMatrix &A1,
		SquareMatrix& EigenVecs);
#endif /* SYMMETRICMATRIX_H */

