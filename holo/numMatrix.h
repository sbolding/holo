/* 
 numMatrix

 This class is the interface to the derived matrix classes.
*/

#ifndef numMatrix_h
#define numMatrix_h

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "numVector.h"

class numMatrix{
  protected:
    double **_coeff;		// Coefficients of matrix
    int _n_rows, _n_cols;        // Number of columns and rows
	// Never to be used copy constructor
	numMatrix(const numMatrix &matrix);
	// Never to be used Constructor
	numMatrix();
  public:
    //Destructor
    virtual ~numMatrix();
    // Functions
	void zero();
	int getNumRows() const;
	int getNumCols() const;
	void print(std::ostream &out);
    void printCompact(std::ostream &out);
    void scale(double value);
    virtual void setCoeff(int i, int j, double value) = 0;
	virtual void addCoeff(int i, int j, double value) = 0;
	virtual double getCoeff(int i, int j) = 0;
	virtual void mult(numMatrix *b, numMatrix *c) = 0;
	virtual void mult(numVector *b, numVector *c) = 0;
	virtual void trans(numMatrix *a) = 0;
	virtual void gauss(numVector *b, numVector *x) = 0;
};

#endif
