#include "SourceConstant.h"

SourceConstant::SourceConstant()
{
	//should never be called
}


SourceConstant::SourceConstant(double value) : Source()
{
	_value = value;
}

void SourceConstant::sampleSource()
{
	//do nothing yet;
}