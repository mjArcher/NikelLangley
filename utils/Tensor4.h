#ifndef __TENSOR4_H
#define __TENSOR4_H

/* #include "blitz/tinyvec.h" */
/* #include "blitz/tinyvec-et.h" */
#include "Invariants.h"
#include "SquareMatrix.h"
#include "SymmetricMatrix.h"

using namespace blitz;

//! Fourth order tensor
class Tensor4
{

    double Data_[3][3][3][3];

public:

    Tensor4() {}

    explicit Tensor4 ( const double scal )  {
      for (int i=0; i<3; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<3; ++k) {
            for (int l=0; l<3; ++l) {
	      Data_[i][j][k][l] = scal;
            }
	  }
	}
      }      
    }

    const double& operator()(int i, int j, int k, int l) const
    {
      return Data_[i][j][k][l];
    }

    double& operator()(int i, int j, int k, int l)
    {
      return Data_[i][j][k][l];
    }    

    SquareMatrix operator()(int i, int j) const
    {
      SquareMatrix M;
      for(int k=0;k<3;++k) 
      {
	for(int l=0;l<3;++l) 
	{
	  M(k,l) = Data_[i][j][k][l];
	}
      }
      return M;
    }

    const Tensor4& operator*=(const double s) 
    {
      for (int i=0; i<3; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<3; ++k) {
	    for (int l=0; l<3; ++l) {
	      Data_[i][j][k][l] *= s;	    
	    }
	  }
	}	
      }      
      return *this;
    }
    
};

//! Fourth order tensor
class Symmetric12Tensor4
{

    double Data_[6][3][3];
    static const int Lookup_[3][3];

public:

    Symmetric12Tensor4() {}

    explicit Symmetric12Tensor4 ( const double scal )  {
      for (int i=0; i<6; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<3; ++k) {
	    Data_[i][j][k] = scal;	    
	  }
	}
      }      
    }

    inline double& operator() ( const int i, const int j, const int k, const int l )
    {
      return Data_[Lookup_[i][j]][k][l];
    }

    inline const double& operator() ( const int i, const int j, const int k, const int l ) const
    {
      return Data_[Lookup_[i][j]][k][l];
    }

    SquareMatrix operator()(int i, int j) const
    {
      SquareMatrix M;
      for(int k=0;k<3;++k) 
      {
	for(int l=0;l<3;++l) 
	{
	  M(k,l) = Data_[Lookup_[i][j]][k][l];
	}
      }
      return M;
    }

    const Symmetric12Tensor4& operator*=(const double s) 
    {
      for (int i=0; i<3; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<6; ++k) {
	      Data_[k][i][j] *= s;	    
	  }
	}	
      }      
      return *this;
    }

};

//! Fourth order tensor
class Symmetric34Tensor4
{

    double Data_[3][3][6];
    static const int Lookup_[3][3];

public:

    Symmetric34Tensor4() {}

    explicit Symmetric34Tensor4 ( const double scal )  {
      for (int i=0; i<3; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<6; ++k) {
	    Data_[i][j][k] = scal;	    
	  }
	}
      }      
    }

    inline double& operator() ( const int i, const int j, const int k, const int l )
    {
      return Data_[i][j][Lookup_[k][l]];
    }

    inline const double& operator() ( const int i, const int j, const int k, const int l ) const
    {
      return Data_[i][j][Lookup_[k][l]];
    }

    SymmetricMatrix operator()(int i, int j) const
    {
      SymmetricMatrix M;
      for(int k=0;k<3;++k) 
      {
	for(int l=k;l<3;++l) 
	{
	  M(k,l) = Data_[i][j][Lookup_[k][l]];
	}
      }
      return M;
    }

    const Symmetric34Tensor4& operator*=(const double s) 
    {
      for (int i=0; i<3; ++i) {
	for (int j=0; j<3; ++j) {
	  for (int k=0; k<6; ++k) {
	      Data_[i][j][k] *= s;	    
	  }
	}	
      }      
      return *this;
    }

};

//! Fourth order tensor
class Symmetric1234Tensor4
{

    double Data_[6][6];
    static const int Lookup_[3][3];

public:

    Symmetric1234Tensor4() {}

    explicit Symmetric1234Tensor4 ( const double scal )  {
      for (int i=0; i<6; ++i) {
	for (int j=0; j<6; ++j) {
	  Data_[i][j] = scal;	    
	}
      }      
    }

    inline double& operator() ( const int i, const int j, const int k, const int l )
    {
      return Data_[Lookup_[i][j]][Lookup_[k][l]];
    }

    inline const double& operator() ( const int i, const int j, const int k, const int l ) const
    {
      return Data_[Lookup_[i][j]][Lookup_[k][l]];
    }

    SymmetricMatrix operator()(int i, int j) const
    {
      SymmetricMatrix M;
      for(int k=0;k<3;++k) 
      {
	for(int l=k;l<3;++l) 
	{
	  M(k,l) = Data_[Lookup_[i][j]][Lookup_[k][l]];
	}
      }
      return M;
    }

    const Symmetric1234Tensor4& operator*=(const double s) 
    {
      for (int i=0; i<6; ++i) {
	for (int j=0; j<6; ++j) {
	  Data_[i][j] *= s;	    
	}
      }      
      return *this;
    }
    
};

#endif
