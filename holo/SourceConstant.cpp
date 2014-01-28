#include "SourceConstant.h"

SourceConstant::SourceConstant()
{
	//should never be called
}


SourceConstant::SourceConstant(double value) : Source()
{
	_value = value;
}

void SourceConstant::sampleSource(Particle1D * particle, Mesh* mesh)
{
	//For a constant source just need to sample which cell it is 
	part


}