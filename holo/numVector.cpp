// numVector class definition

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
//#include "Define.h"
#include "numMatrix.h"
#include "numVector.h"
using namespace std;

numVector::numVector()
{
	_n_rows = 0;
	_coeff = NULL;
}

numVector::numVector(int nr)
{
	_n_rows = nr;

	// Allocate pointer to rows
	_coeff = new double[_n_rows];
	if(_coeff == NULL) {
		cout << "allocation failure in numVector";
		exit(0);
	}
	zero();
}

//Destructor
numVector::~numVector()
{
	delete [] _coeff;
}		

// Functions
void numVector::zero()
{
	for(int i=0; i<_n_rows; i++){
		_coeff[i] = 0.;
	}
}

int numVector::getNumRows() const
{
	return _n_rows;
}

void numVector::setCoeff(int i, double value)
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in SetCoeff";
		exit(0);
	}
	_coeff[i] = value;
}

void numVector::addCoeff(int i, double value)
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in AddCoeff";
		exit(0);
	}
	_coeff[i] += value;
}

double numVector::getCoeff(int i)
{
	if((i<0) || (i>_n_rows)){
		cout << "Error in GetCoeff";
		exit(0);
	}
	return _coeff[i];
}

void numVector::print(ostream &out)
{
	for(int i=0; i<_n_rows; i++){
		out << "Row " << i << "  ";
		out.setf(ios::scientific);
		out.precision(3);
		out.width(12);
		out << _coeff[i] << endl;	
	}
    out << endl;
}
