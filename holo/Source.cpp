#include "Source.h"
#include "Particle1D.h"
#include <iostream>

Source::Source()
{
	//do nothing;
}

Source::Source(Particle1D* particle, string sampling_method)
{
	_particle = particle;
	_rng = particle->_rng;
	_sampling_method = HoMethods::sampling_map.at(sampling_method);
}

void Source::sampleAngleIsotropic()
{
	_particle->_mu = _rng->rand_num()*2.0 - 1.;
}

void Source::sampleAngleIsotropic(double min_cosine, double max_cosine, bool directional_bias)
{
	if (min_cosine > max_cosine)
	{
		std::cerr << "Passed a min cosine greater than max cosine to sampling routine" << std::endl;
		exit(1);
	}
	_particle->_mu = _rng->rand_num()*(max_cosine - min_cosine) + min_cosine;

	//Adjust weight to be unbiased if necessary
	if (!directional_bias)
	{
		_particle->_weight *= (max_cosine - min_cosine)*0.5;
	}

}


double Source::sampleLinDiscFunc(std::vector<double> nodal_values, double left_node_coor, double right_node_coor)
{
	std::cout << "This function has not been checked yet" << std::endl;
	system("pause");
	exit(1); 

	//Sample the position based on the nodal values, should write a function to get the area of the source from the element somehow
	double coordinate;
	double width = (right_node_coor - left_node_coor);
	if (width > 0.0)
	{
		std::cerr << "Passed in coordinates in reverse order to sampleLinDiscFunc, in source.cpp" << std::endl;
		exit(1);
	}

	//If this routine is too slow, do a soft check to see if they are different first, then do the check below
	if (abs(nodal_values[0] - nodal_values[1]) / nodal_values[0] < 1.E-10) //then effectively a constant source, sampling is uniform across the cell
	{
		coordinate = _rng->rand_num()* width + left_node_coor;
	}
	else //need to sample from lin discontinuous source //THIS ROUTINE WORKS
	{
		double left_hat, right_hat; //Normalized nodal values, such that CDF is normalized
		left_hat = 2.0*nodal_values[0] / (nodal_values[1] + nodal_values[0]);
		right_hat = 2.0 - left_hat;
		//use direct inversion of CDF to sample position, based on quadratic formula
		coordinate = -left_hat + sqrt(left_hat*left_hat + 2 * _rng->rand_num()*(right_hat - left_hat));
		coordinate /= (right_hat - left_hat);
		coordinate = coordinate*width + left_node_coor; //convert to mfp
	}
	return coordinate;
}