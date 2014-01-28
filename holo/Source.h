#ifndef _SOURCE_H
#define _SOURCE_H

#include <vector>
#include "Particle1D.h"

class Source
{
private:
	
public:
	Source();
	virtual void sampleSource(Particle1D *p, Mesh * mesh) = 0;

};








#endif