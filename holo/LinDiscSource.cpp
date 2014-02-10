#include "Source.h"
#include "LinDiscSource.h"
#include "Particle1D.h"
#include "Controller.h"

LinDiscSource::LinDiscSource(Particle1D* particle, string sampling_method) : Source(particle, sampling_method)
{
	//need to initialize alias sampler
	//get the area of the source, and the total external source nodal values

	//local variables
	std::vector<Element*>* elements;
	elements = _particle->_mesh->getElements();
	std::vector<double> source_strength_per_cell;
	std::vector<Element*>::iterator it_el;          //element iterator
	it_el = elements->begin();			            //initialize iterator
	double ext_source_el;		             		//magnitude of external source of curr element
	std::vector<double> total_src_nodal_values_el;

	//Local variables to prevent repeated accessing of particle class
	double method = _particle->_method;

	for (; it_el != elements->end(); it_el++)
	{
		//Initialize to the external source strength
		try
		{
			total_src_nodal_values_el = (*it_el)->getExtSourceNodalValues(); //initialize to ext source values
			if (method == HoMethods::HOLO_ECMC || method == HoMethods::HOLO_STANDARD_MC) //append scattering source
			{
				double sigma_s_el = (*it_el)->getMaterial().getSigmaS();
				std::vector<double> scat_src_nodal_values_el = (*it_el)->getScalarFluxNodalValues();//This should return 0 if LO system hasnt been solved yet
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
			ext_source_el = getAreaLinDiscFunction(total_src_nodal_values_el, (*it_el)->getElementDimensions()[0]);
			_total_src_nodal_values.push_back(total_src_nodal_values_el);
			source_strength_per_cell.push_back(ext_source_el);
			_vol_src_total += ext_source_el;
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

void LinDiscSource::sampleSourceParticle()
{
	//Determine if it is volumetric source, or surface source (depending on the mode you are in, may sample scattering source as well)
	//Store the entire source (ext + scattering) into the other one and compute its area.  With the area you can easily determine if sample
	//is from isotropic source or if it is from 
	if (_rng->rand_num() < _vol_src_total / (_BC_src_total + _vol_src_total))
	{
		_particle->_current_element = _alias_sampler->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		sampleLinDiscSource(_total_src_nodal_values[_particle->_current_element]); 
		//sampleLinDiscontSource(_mesh->getElement(_current_element)->getExtSourceNodalValues()); //THIS MIGHT BE USEFUL IN A STANDARD MC CALC, but is essentially equivalent
	}
	else //Volumetric source
	{
		std::cerr << "Sampling of BC source is not yet implemented" << std::endl;
		exit(1);
	}

	//Assume all sources have isotropic angular distribution TODO this will need to be updated for the case of boundary sources, just add over
	//loaded function for sampleAngleIsotropic that takes a start and end range of angle to be allowed in and adjusts weight accordingly
	Source::sampleAngleIsotropic();

	//Update particle properties for the new cell
	_particle->updateElementProperties();
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
