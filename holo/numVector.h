/* 
 numVector

 This class implements a vector.
*/

#ifndef numVector_h
#define numVector_h

#include <stdio.h>
#include <stdlib.h>
#include <ostream>

class numVector{
  protected:
    double *_coeff;		// Coefficients of matrix
    int _n_rows;       		// Number of rows
	//never to be used copy constructor
	numVector(const numVector &vector);
  public:
	// Constructors
	numVector();
	numVector(int nr);

    //Destructor
    virtual ~numVector();
    // Functions
	void zero();
	int getNumRows() const;
	void setCoeff(int i, double value);
	void addCoeff(int i, double value);
	double getCoeff(int i);
	void print(std::ostream &out);
};

#endif
