#include "ResidualSource.h"
#include "Particle1D.h"

ResidualSource::ResidualSource(Particle1D* particle, string sampling_method) : Source(particle, sampling_method)
{ 
	//local variables
	std::vector<ECMCElement1D*>* elements;
	elements = _particle->_mesh->getElements();
	std::vector<ECMCElement1D*>::iterator it_el;          //element iterator
	it_el = elements->begin();			            //initialize iterator
	double ext_source_el;		             		//magnitude of external source of curren ECMC element, units of p/sec
	std::vector<double> total_src_nodal_values_el;   //external source nodal values
	std::vector<double> residual_element_magnitudes;  //magnitue of residual source over elements
	std::vector<double> residual_face_magnitudes; //magnitude of the residual for each face

	//Local variables to prevent repeated accessing of particle class
	Element* spatial_element; //spatial element corresponding to the current element
	std::vector<double> residual_element_LD_values; //the dof of the residual for a particular element

	for (; it_el != elements->end(); it_el++)
	{
		//Map external source (and scattering source if ECMC) onto each element
		spatial_element = (*it_el)->getSpatialElement();
		mapExtSrcToElement(total_src_nodal_values_el, ext_source_el, spatial_element, *it_el); //determine external source nodal values over the element
		convertNodalToMoments(total_src_nodal_values_el, residual_element_LD_values); 
		
		








		_vol_src_total += ext_source_el; //units of p / sec

	
	}

	if (_sampling_method == HoMethods::STANDARD_SAMPLING) //standard source sampling
	{
		//Create sampler with alias sampling, let it normalize, delete unneccessary data
		//_alias_sampler = new AliasSampler(source_strength_per_cell, false);
	}
	else
	{
		std::cerr << "No other samplilng methods are implemented yet" << std::endl;
		exit(1);
	}

	//Initially assume no BC source TODO
	_BC_src_total = 0.0; 
}

void ResidualSource::sampleSourceParticle()
{

}