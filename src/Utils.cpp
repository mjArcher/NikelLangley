#include "Utils.h"

using namespace std;

void printDoubleArray(double * array, int N)
{
	for(int i = 0; i < N; i++)
	{
		cout << array[i] << '\t' << i << endl;
	}
}

string convertToStr(int num)
{
	ostringstream convert;
	convert << num; 
	string number = convert.str();
	convert.str("");
	convert.clear();
	return number;
}

double kron(int j, int k)
{
	if (j == k)
		return 1.;
	else
		return 0;
}

