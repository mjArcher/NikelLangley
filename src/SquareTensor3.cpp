#include "SquareTensor3.h"

using namespace std;
using namespace Eigen;

//multiply 3rd order tensor by a mtrix
SquareTensor3 operator*(const SquareTensor3& lhs, const Matrix3d& s)
{
	return s*lhs;
}
//const safety needs to be added
SquareTensor3 operator*(const Matrix3d& s, const SquareTensor3& rhs)
{
	SquareTensor3 C(rhs);
	C *= s;
	return C;
}

ostream& operator<<(ostream& os, const SquareTensor3& t) 
{
	for(int i =0; i < 3; i++)
	{
		std::cout << t[i] << std::endl;
    std::cout << "\n";
	}
 /* std::cout << "done " << std::endl; */
	return os;
};




