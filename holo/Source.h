#ifndef _SOURCE_H
#define _SOURCE_H

#include <vector>
#include "RNG.h"

class Particle1D; //forward declaration

class Source
{
protected:

	Particle1D* _particle; //pointer to the particle class, source is a friend of particle class
	Source(); //Never use the default constructor
	RNG* _rng;
	unsigned int _sampling_method; //0 for alias sampling, 1 for stratified sampling
	double _vol_src_total; //The total volumetric source (p/sec)
	double _BC_src_total; //The total BC source (p/sec)

	//Angle sampling methods
	void sampleAngleIsotropic();
	void sampleAngleIsotropic(double min_cosine, double max_cosine, bool directional_bias=false);  //will adjust weight to produce unbiased particle if source distribution is uniform over sphere, i.e., 1/4PI (default)

	//Spatial sampling methods
	void sampleLinDiscontSource(std::vector<double>); //method used to sample lin discontinuous sources
	//void sampleConstExtSource();	//Don't know if this is actually ever needed? probably not
	void initializeSamplingSource(); //will create the total source vector, as well as initilize sampling routines, is a unique function because source can be reset between cycles

public:
	Source(Particle1D* particle, string sampling_method);
	virtual void sampleSourceParticle() = 0; //pure virtual class, samples a source particle direction and 
	double getAreaLinDiscFunction(std::vector<double> nodal_values, double element_volume) //for computing total source magnitudes
	{
		return 0.5*(nodal_values[0] + nodal_values[1])*element_volume;
	};

};








#endif