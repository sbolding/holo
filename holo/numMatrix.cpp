// numMatrix class definition

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
//#include "Define.h"
#include "numMatrix.h"

using std::ostream;
using std::endl;
using std::cout;
using std::ios;

// Constructors
numMatrix::numMatrix()
{
	_n_rows = 0;
	_n_cols = 0;
	_coeff = NULL;
}

//Destructor
numMatrix::~numMatrix()
{
	for(int i=0; i<_n_rows; i++) {
		delete [] _coeff[i];
	}
	delete [] _coeff;
}		

// Functions
void numMatrix::zero()
{
	for(int i=0; i<_n_rows; i++){
		for(int j=0; j<_n_cols; j++){
			_coeff[i][j] = 0.;
		}
	}
}

void numMatrix::scale(double value)
{
	for(int i=0; i<_n_rows; i++){
		for(int j=0; j<_n_cols; j++){
			_coeff[i][j] *= value;
		}
	}
}


int numMatrix::getNumRows() const
{
	return _n_rows;
}

int numMatrix::getNumCols() const
{
	return _n_cols;
}

void numMatrix::print(ostream &out){
	for(int i=0; i<_n_rows; i++){
		out << "Row " << i << endl;
		for(int j=0; j<_n_cols; j++){
			out.setf(ios::scientific);
			out.precision(4);
			out.width(13);
			out << _coeff[i][j];
			if((j+1)%6 == 0){
				out << endl;
			}
		}
		out << endl;	
		out << endl;	
	}
}

void numMatrix::printCompact(ostream &out){
	for(int i=0; i<_n_rows; i++){
		for(int j=0; j<_n_cols; j++){
			out.setf(ios::scientific);
			out.precision(3);
			out.width(12);
			out << _coeff[i][j];
		}
		out << endl;		
	}
    out << endl;
}


