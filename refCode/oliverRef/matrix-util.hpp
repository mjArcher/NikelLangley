#ifndef MATRIX_UTIL_HPP
#define MATRIX_UTIL_HPP

#include <tvmet/Vector.h>
#include <tvmet/Matrix.h>
#include <math.h>
#include <utility> // pair
#include <algorithm> // max

// Find Givens rotation (cos, sin) pair.  Used by jacobi_eig_sym.
template <typename T> std::pair<T, T> approx_givens_angle(T a11, T a12, T a22)
{
	bool b = a12*a12 < (a11-a22)*(a11-a22);
	T omega = 1.0/sqrt(a12*a12 + (a11-a22)*(a11-a22));
	T s = b?(omega*a12):(sqrt(0.5));
	T c = b?(omega*(a11-a22)):(sqrt(0.5));
	return std::make_pair(c,s);
}

// Jacobi Eigendecomposition of the symmetric matrix A.  The
// eigenvalues are output to 'lambda' and the right eigenvectors to
// 'right'.  Returns the maximum off-diagonal element of the
// eigenvalue matrix, as a measure of the error (i.e., smaller is
// better).
template <class T> T jacobi_eig_sym(tvmet::Matrix<T,3,3> A,
				    tvmet::Vector<T,3>& lambda,
				    tvmet::Matrix<T,3,3>& right)
{
	const int num_sweeps = 3;
	
	right = tvmet::identity<T,3,3>();
	
	tvmet::Matrix<T,3,3> Q;
	std::pair<double, double> givens;

	for (int n=0; n<num_sweeps; n++) {
		//////////////
		// p=0; q=1;
		givens = approx_givens_angle(A(0,0), A(0,1), A(1,1));
		double c = givens.first;
		double s = givens.second;

		Q = c,-s, 0,
		    s, c, 0,
		    0, 0, 1;

		// tvmet cannot update A in place - need a temporary here 
		tvmet::Matrix<T,3,3> A_next;
		A_next = tvmet::trans(Q) * A * Q;
		A = A_next;

		tvmet::Matrix<T,3,3> right_next;
		right_next = right * Q;
		right = right_next;

		//////////////
		// p=0; q=2;
		givens = approx_givens_angle(A(0,0), A(0,2), A(2,2));
		c = givens.first;
		s = givens.second;

		Q = c, 0,-s,
		    0, 1, 0,
		    s, 0, c;
		
		A_next = tvmet::trans(Q) * A * Q;
		A = A_next;
		right_next = right * Q;
		right = right_next;

		//////////////
		// p=1; q=2;
		givens = approx_givens_angle(A(1,1), A(1,2), A(2,2));
		c = givens.first;
		s = givens.second;

		Q = 1, 0, 0,
	 	    0, c,-s,
		    0, s, c;

		A_next = tvmet::trans(Q) * A * Q;
		A = A_next;
		right_next = right * Q;
		right = right_next;

	}
	lambda[0] = A(0,0);
	lambda[1] = A(1,1);
	lambda[2] = A(2,2);

	// the maximum off diagonal element;
	return std::max(A(0,1),std::max(A(1,0),A(1,2))); 
}

#endif
