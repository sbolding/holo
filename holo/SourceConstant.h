#ifndef _SOURCECONSTANT_H
#define _SOURCECONSTANT_H

#include "Source.h"

class SourceConstant : public Source
{
private:

	double _value; //constant source value
	SourceConstant(); //default constructor, never use

public:

	SourceConstant(double value); //Constant source value initializer
	virtual void sampleSource(Particle1D * p, Mesh* mesh);  //overwrite sample source function.  Sample is done to particle p

};


#endif