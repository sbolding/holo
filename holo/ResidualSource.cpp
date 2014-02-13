#include "ResidualSource.h"
#include "Particle1D.h"

ResidualSource::ResidualSource(Particle1D* particle, string sampling_method) : Source(particle, sampling_method)
{ 
	//local variables
	std::vector<ECMCElement1D*>* elements;
	elements = _particle->_mesh->getElements();
	std::vector<ECMCElement1D*>::iterator it_el;          //element iterator
	it_el = elements->begin();			            //initialize iterator
	double res_element_mag_el;		         //mag of element residual
	double res_face_mag_el;					//magnitude of residaul on face
	std::vector<double> res_element_mags;  //magnitue of residual source over elements
	std::vector<double> res_face_mags; //magnitude of the residual for each face

	//Local variables to prevent repeated accessing of particle class
	std::vector<double> res_LD_values_el; //the dof of the residual for a particular element


	for (; it_el != elements->end(); it_el++)
	{
		//compute volume integrals //IF CELL IS ACTIVE
		computeElementResidual(*it_el, res_LD_values_el,res_element_mag_el);









		_vol_src_total += res_element_mag_el; //units of p / sec

	
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

void ResidualSource::computeElementResidual(ECMCElement1D* element, 
	std::vector<double> & res_LD_values_el, double & res_mag)
{
	//initialize variables
	res_mag = 0.0;
	res_LD_values_el.assign(3,0.0);
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate();
	double h_x = element->getSpatialWidth();
	std::vector<double> ang_flux_dof = element->getAngularFluxDOF();
	double & psi_avg = ang_flux_dof[0]; //references for better legibility
	double & psi_x = ang_flux_dof[1];
	double & psi_mu = ang_flux_dof[2];

	//Local variables to prevent repeated accessing of particle class
	Element* spatial_element; //spatial element corresponding to the current element
	std::vector<double> ext_source_LD_values; //the external source
	double ext_source_mag;
	std::vector<double> total_src_nodal_values_el;   //external source nodal values

	//Map external source (and scattering source if ECMC) onto each element
	spatial_element = element->getSpatialElement();
	double sigma_tot = spatial_element->getMaterial().getSigmaT(); //total cross section
	mapExtSrcToElement(total_src_nodal_values_el, ext_source_mag, spatial_element, element); //determine external source nodal values over the element
	convertNodalValuesToMoments(total_src_nodal_values_el, ext_source_LD_values, true); //if isotropic special case, Q_mu = 0.0

	//compute residual LD values for the element
	res_LD_values_el[0] = ext_source_LD_values[0] - sigma_tot*psi_avg
		- dir_coor*2. / h_x*psi_x;  //average term
	res_LD_values_el[1] = ext_source_LD_values[1] - sigma_tot*psi_x; //x moment
	res_LD_values_el[2] = ext_source_LD_values[2] - sigma_tot*psi_mu
		- psi_x*h_mu / h_x;  //mu moment

	//temporary computation, totally wrong
	std::cout << "this is so wrong " << std::endl;
	res_mag = h_x*h_mu*res_LD_values_el[0];


}