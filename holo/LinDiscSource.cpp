#include "Source.h"
#include "LinDiscSource.h"
#include "Particle1D.h"
#include "Controller.h"

LinDiscSource::LinDiscSource(Particle1D* particle, string sampling_method) : Source(particle, sampling_method)
{
	//need to initialize alias sampler
	//get the area of the source, and the total external source nodal values

	//local variables
	std::vector<ECMCElement1D*>* elements;
	elements = _particle->_mesh->getElements();
	std::vector<double> source_strength_per_cell;
	std::vector<ECMCElement1D*>::iterator it_el;          //element iterator
	it_el = elements->begin();			            //initialize iterator
	double ext_source_el;		             		//magnitude of external source of curren ECMC element, units of p/sec
	std::vector<double> total_src_nodal_values_el;

	Element* spatial_element; //spatial element corresponding to the current element

	for (; it_el != elements->end(); it_el++)
	{
		try
		{
			//Map external source (and scattering source if ECMC) onto each element
			spatial_element = (*it_el)->getSpatialElement();
			mapExtSrcToElement(total_src_nodal_values_el, ext_source_el, spatial_element, *it_el); //determine external source nodal values over the element

			//Add element values to total array
			_total_src_nodal_values.push_back(total_src_nodal_values_el); 
			source_strength_per_cell.push_back(ext_source_el);
			_vol_src_total += ext_source_el; //units of p / sec
		}
		catch (...)
		{
			std::cerr << "The HoSolver had trouble initialized because the Lo System was not properly initialized and not solved, in Particle1D" << std::endl;
			exit(1);
		}
	}

	if (_sampling_method == HoMethods::STANDARD_SAMPLING) //standard source sampling
	{
		//Create sampler with alias sampling, let it normalize, delete unneccessary data
		_alias_sampler = new AliasSampler(source_strength_per_cell, false);
	}
	else
	{
		std::cerr << "No other samplilng methods are implemented yet" << std::endl;
		exit(1);
	}

	//Initially assume no BC source TODO
	_BC_src_total = 0.0;
}

LinDiscSource::~LinDiscSource()
{
	delete _alias_sampler;
}

void LinDiscSource::sampleSourceParticle()
{
	//Determine if it is volumetric source, or surface source (depending on the mode you are in, may sample scattering source as well)
	//Store the entire source (ext + scattering) into the other one and compute its area.  With the area you can easily determine if sample
	//is from isotropic source or if it is from boundary
	if (_rng->rand_num() < _vol_src_total / (_BC_src_total + _vol_src_total))
	{
		_particle->_current_element_ID = _alias_sampler->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties(); 
		sampleLinDiscSource(_total_src_nodal_values[_particle->_current_element_ID]); //sample position in bin
		
		//determine direction within the element
		double mu_center = _particle->_current_element->getAngularCoordinate();
		double half_angular_width = 0.5*_particle->_current_element->getAngularWidth();

		sampleAngleIsotropic(mu_center - half_angular_width, mu_center + half_angular_width); //sample isotropically within the bin
		
	}
	else //Boundary Source
	{
		std::cerr << "Sampling of BC source is not yet implemented" << std::endl;
		exit(1);
	}

}

void LinDiscSource::sampleLinDiscSource(std::vector<double> nodal_values)
{
	//Sample the position based on the nodal values, should write a function to get the area of the source from the element somehow

	//If this routine is too slow, do a soft check to see if they are different first, then do the check below
	if (abs(nodal_values[0] - nodal_values[1]) / nodal_values[0] < 1.E-10) //then effectively a constant source, sampling is uniform across the cell
	{
		_particle->_position_mfp = _rng->rand_num()*_particle->_element_width_mfp;
	}
	else //need to sample from lin discontinuous source //THIS ROUTINE WORKS
	{
		double left_hat, right_hat; //Normalized nodal values, such that CDF is normalized
		left_hat = 2.0*nodal_values[0] / (nodal_values[1] + nodal_values[0]);
		right_hat = 2.0 - left_hat;
		//use direct inversion of CDF to sample position, based on quadratic formula
		_particle->_position_mfp = -left_hat + sqrt(left_hat*left_hat + 2 * _rng->rand_num()*(right_hat - left_hat));
		_particle->_position_mfp /= (right_hat - left_hat);
		_particle->_position_mfp *= _particle->_element_width_mfp; //convert to mfp
	}
}
