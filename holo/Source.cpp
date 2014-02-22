#include "Source.h"
#include "Particle1D.h"
#include <iostream>

Source::Source()
{
	//do nothing;
}

Source::~Source()
{
	//do nothing
}

Source::Source(Particle1D* particle, string sampling_method)
{
	_particle = particle;
	_rng = particle->_rng; 
	_sampling_method = HoMethods::sampling_map.at(sampling_method);
	_vol_src_total = 0.0;
	_BC_src_total = 0.0;
}

void Source::sampleAngleIsotropic()
{
	_particle->_mu = _rng->rand_num()*2.0 - 1.;
}

void Source::sampleAngleIsotropic(double min_cosine, double max_cosine, bool directional_bias)
{
	//This function works correctly for - or + numbers, as long as min_cosine < max_cosine, e.g., -1 < -0.5
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

void Source::sampleAngleCosineLaw(int n, double min_cosine, double max_cosine)
{
	if (min_cosine > max_cosine || n < 0 || min_cosine*max_cosine < 0.0)
	{
		if (std::abs(min_cosine) > GlobalConstants::RELATIVE_TOLERANCE
			&& std::abs(max_cosine) > GlobalConstants::RELATIVE_TOLERANCE) //special case for roundoff error (0 can be negative)
		{
			std::cerr << "Passed a min cosine greater than max cosine to sampling routine, or cell across zero" << std::endl;
			exit(1);
		}
	}
	//sample from pdf f(mu) = norm_const*mu^n, mu in [min_cosine, max_cosine]
	double norm_const = (std::pow(max_cosine, n + 1) - std::pow(min_cosine, (n + 1)));
	if (min_cosine >= 0.0)
	{
		_particle->_mu = _rng->rand_num()*norm_const + pow(min_cosine, n + 1);
		_particle->_mu = pow(_particle->_mu, 1. / (n + 1.));
	}
	else //mu is negative and odd function so special case
	{
		_particle->_mu = _rng->rand_num()*abs(norm_const) + std::abs(pow(max_cosine, n + 1));
		_particle->_mu = -1.*pow(_particle->_mu, 1. / (n + 1.));
	}

}

void Source::sampleAngleCosineLaw(double min_cosine, double max_cosine)
{
	if (min_cosine > max_cosine || min_cosine*max_cosine < 0.0)
	{
		if (std::abs(min_cosine) > GlobalConstants::RELATIVE_TOLERANCE
			&& std::abs(max_cosine) > GlobalConstants::RELATIVE_TOLERANCE) //special case for roundoff error (0 can be negative)
		{
			std::cerr << "Passed a min cosine greater than max cosine to sampling routine, or cell across zero" << std::endl;
			exit(1);
		}
	}
	
	if (min_cosine < -1.*GlobalConstants::RELATIVE_TOLERANCE) //for negative case f(mu) = nrom_const*|mu|, mu in [min_cosine, max_cosine]
	{
		_particle->_mu = -1.*std::sqrt(_rng->rand_num()*(min_cosine*min_cosine-
			max_cosine*max_cosine) + max_cosine*max_cosine);
	}
	else //sample from pdf f(mu) = norm_const*mu, mu in [min_cosine, max_cosine]
	{
		_particle->_mu = std::sqrt(_rng->rand_num()*(max_cosine*max_cosine -
			min_cosine*min_cosine) + min_cosine*min_cosine);
	}
}

double Source::sampleLinDiscFunc1D(std::vector<double> nodal_values, double left_node_coor, double right_node_coor)
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

void Source::mapExtSrcToElement(std::vector<double> & total_src_nodal_values_el, double & tot_src_strength,
	Element* spatial_element, ECMCElement1D* element)
{
	double angular_probability = 0.5*element->getAngularWidth();

	total_src_nodal_values_el = spatial_element->getExtSourceNodalValues(); //initialize to ext source values, this has units of particles/sec-cm
	if (_particle->_method == HoMethods::HOLO_ECMC || _particle->_method == HoMethods::HOLO_STANDARD_MC) //append scattering source
	{
		double sigma_s_el = spatial_element->getMaterial().getSigmaS();
		std::vector<double> scat_src_nodal_values_el = spatial_element->getScalarFluxNodalValues();//This should return 0 if LO system hasnt been solved yet
		if (scat_src_nodal_values_el.size() != total_src_nodal_values_el.size())
		{
			std::cerr << "Scattering source and external source do not have the same number of nodal values" << std::endl;
			exit(0);
		}
		for (int node = 0; node < scat_src_nodal_values_el.size(); node++)
		{
			total_src_nodal_values_el[node] += scat_src_nodal_values_el[node] * sigma_s_el; //phi*_sigma_s, note there is no 1/(4pi) here, because we want (p/sec)
		}
	}
	tot_src_strength = getAreaLinDiscFunction(total_src_nodal_values_el, spatial_element->getElementDimensions()[0])
		* angular_probability; //fraction in this angular element, units of p / (sec-str)
	for (int node = 0; node < total_src_nodal_values_el.size(); ++node) //convert to per steradian
	{
		total_src_nodal_values_el[node] *= 0.5; //azimutthally integrated angular flux values
	}
}

void Source::convertNodalValuesToMoments(std::vector<double> & nodal_values,
	std::vector<double> & ld_moments, bool nodal_values_isotropic)
{
	if (nodal_values.size() != 2)
	{
		std::cerr << "This method is only implemented for 1D, in Source::convertNodalValues\n";
		system("pause");
		exit(1);
	}
	//Set ld_moments to have one more DOF than nodal values
	ld_moments.clear();
	ld_moments.assign(nodal_values.size() + 1, 0.0);

	//compute average and spatial moment
	ld_moments[0] = 0.5*(nodal_values[0] + nodal_values[1]);
	ld_moments[1] = 0.5*(nodal_values[1] - nodal_values[0]);
	if (!nodal_values_isotropic)
	{
		std::cerr << "I have not implemented this method yet" << std::endl;
		system("pause");
		exit(1);
	}
	else //isotropic
	{
		//angular moment is zero because uniform value in angle, the nodal values have taken into account angular fraction already	
	}
}