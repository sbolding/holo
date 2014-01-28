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
	virtual void sampleSource();  //overwrite sample source function


};


#endif