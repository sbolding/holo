#include "ResidualSource.h"
#include "Particle1D.h"

ResidualSource::ResidualSource(Particle1D* particle, string sampling_method) : Source(particle, sampling_method)
{ 
	//member variables initialize
	_face_src_total = 0.0;

	//local variables
	std::vector<ECMCElement1D*>* elements;
	std::vector<ECMCElement1D*>::iterator it_el;
	elements = _particle->_mesh->getElements();
	double res_element_mag_el;		         //mag of element residual
	double res_face_mag_el;					//magnitude of residaul on face
	std::vector<double> res_element_mags;  //magnitue of residual source over elements
	std::vector<double> res_face_mags; //magnitude of the residual for each face

	//Local variables to prevent repeated accessing of particle class
	std::vector<double> res_LD_values_el; //the dof of the residual for a particular element
	std::vector<double> res_LD_values_face; //the dof of the residual for the delta source

	//for computing face residuals on boundary
	std::vector<int> boundary_cells = _particle->_mesh->findUpwindBoundaryCells(); 

	//initialize res vectors
	int n_elems = _particle->_mesh->getNumElems();
	res_face_mags.resize(n_elems);
	res_element_mags.resize(n_elems);
	_residual_element_LD_values.resize(n_elems);
	_residual_face_LD_values.resize(n_elems);

	it_el = elements->begin();
	for (; it_el != elements->end(); it_el++)
	{
		//reset magnitudes just in case
		res_face_mag_el = 0.0;
		res_element_mag_el = 0.0;

		//compute volume integrals //IF CELL IS ACTIVE
		computeElementResidual(*it_el, res_LD_values_el, res_element_mag_el);
		computeFaceResidual(*it_el, res_LD_values_face, res_face_mag_el, false); //for non boundary cells

		_residual_element_LD_values[(*it_el)->getID()] = res_LD_values_el;
		if ((*it_el)->getDownStreamElement() != NULL) //in this case there is a down wind cell to add value to
		{
			_residual_face_LD_values[(*it_el)->getDownStreamElement()->getID()] = res_LD_values_face;
		}

_vol_src_total += res_element_mag_el; //units of p / sec
_face_src_total += res_face_mag_el; //units of p / sec
	}

	for (int bc_el = 0; bc_el < boundary_cells.size(); bc_el++)
	{
		computeFaceResidual((*elements)[bc_el], res_LD_values_el, res_face_mag_el, true);
		_face_src_total += res_face_mag_el; //units of p / sec
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
	res_LD_values_el.assign(3, 0.0);
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate();
	double h_x = element->getSpatialWidth();
	std::vector<double> ang_flux_dof = element->getAngularFluxDOF();
	double & psi_avg = ang_flux_dof[0]; //references for better legibility
	double & psi_x = ang_flux_dof[1];
	double & psi_mu = ang_flux_dof[2];
	double & res_avg = res_LD_values_el[0]; //aliases for computing residual magnitude
	double & res_x = res_LD_values_el[1];
	double & res_mu = res_LD_values_el[2];

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

	//compute four corner values to figure out which integral to use
	double r_right_plus = res_LD_values_el[0] + res_LD_values_el[1] + res_LD_values_el[2];
	double r_right_minus = res_LD_values_el[0] + res_LD_values_el[1] - res_LD_values_el[2];
	double r_left_plus = res_LD_values_el[0] - res_LD_values_el[1] + res_LD_values_el[2];
	double r_left_minus = res_LD_values_el[0] - res_LD_values_el[1] - res_LD_values_el[2];

	r_right_plus = +2.0;
	r_right_minus = -3.;
	r_left_plus = +7.0;
	r_left_minus = +5.0;


	//this first if block is from jake's code, not sure why it is needed
	if (abs(res_x)/abs(res_avg) < GlobalConstants::RELATIVE_TOLERANCE) //x_slope~0
	{
		std::cerr << "This code has never actually been called, residualSource.cpp" << std::endl;
		system("pause");
		exit(1);
		if (abs(res_mu)/abs(res_avg) < GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0 also
		{
			res_mag = abs(res_avg);
		}
		else //mu slope non-zero
		{
			
			res_mag = (res_mu*res_mu + res_avg*res_avg) / abs(2.*res_mu);
		}
	}
	else if (abs(res_mu)/abs(res_avg)<GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0
	{
		std::cerr << "This code has never actually been called, residualSource.cpp" << std::endl;
		system("pause");
		exit(1);
		res_mag = (res_x*res_x + res_avg*res_avg) / abs(2.*res_x);
	}


	else if (r_left_plus*r_right_plus > 0.0) //no sign change on top
	{
		if (r_left_minus*r_right_minus > 0.0) //no change on bottom
		{
			if (r_left_plus*r_left_minus > 0.0) //no change at all
			{
				res_mag = abs(res_LD_values_el[0]);
			}
			else if (r_left_plus*r_right_minus < 0.0) //change on left and right
			{
				res_mag = (res_avg*res_avg + res_x*res_x / 3. + res_mu*res_mu) / abs(2.*res_mu);
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else if (r_left_minus*r_right_minus < 0.0) //change on bottom (minus)
		{
			if (r_left_plus*r_right_minus < 0.0) // change on bottom  and right
			{
				res_mag = abs(res_avg + pow(r_right_minus, 3) / (12.*res_x*res_mu));
			}
			else if (r_right_plus*r_left_minus < 0.0) //change on bottom and left
			{
				res_mag = abs(res_avg + pow(r_left_minus, 3) / (12.*res_x*res_mu));
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else
		{
			std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
			system("pause");
			exit(1);
		}
	}
	else if (r_left_plus*r_right_plus < 0.0) //change on top (plus)
	{
		if (r_left_minus*r_right_minus>0.0) //no change on bottom
		{
			if (r_right_plus*r_right_minus < 0.0) //change on top and right
			{
				res_mag = abs(res_avg + pow(r_right_plus, 3) / (12.*res_x*res_mu));
			}
			else if (r_left_plus*r_left_minus < 0.0) //change on top and left
			{
				res_mag = abs(res_avg + pow(r_left_plus, 3) / (12.*res_x*res_mu));
			}
			else
			{
				std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
				system("pause");
				exit(1);
			}
		}
		else // change on top and bottom
		{
			res_mag = (res_avg*res_avg + res_x*res_x + res_mu*res_mu / 3.) / (2 * abs(res_x));
		}
	}
	else
	{
		std::cerr << "You should not be here, error in integral logic, ResidualSource.cpp\n";
		system("pause");
		exit(1);
	}

	//multiply by area
	res_mag *= h_x*h_mu;
	std::cout << "this is so wrong " << std::endl;



}

void ResidualSource::computeFaceResidual(ECMCElement1D* element, std::vector<double> & res_LD_values_face, double & res_mag,
	bool on_boundary)
{
	if (element->getDownStreamElement() == NULL) 
	{
		res_LD_values_face.clear();
		res_mag = 0.;
		return;
	}

	//initialize variables
	res_mag = 0.0;
	res_LD_values_face.assign(3, 0.0);
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate();
	double mu_sgn = dir_coor / abs(dir_coor); //negative or positive direction
	double & res_avg = res_LD_values_face[0];
	double & res_mu = res_LD_values_face[1];

	std::vector<double> psi_down = element->getDownStreamElement()->getAngularFluxDOF(); //downstream values
	std::vector<double> psi_up = element->getAngularFluxDOF(); //upstream element values

	if (on_boundary)
	{
		//TODO Here would need to add boundary term
		psi_down.assign(3, 0.);
	}
	//mu_sgn*psi[1] tells you if on left or right face
	res_avg = (psi_up[0] + mu_sgn*psi_up[1]) - (psi_down[0] - mu_sgn*psi_down[1]);
	res_mu = (psi_up[2] - psi_down[2]);

	//compute magnitude of integral
	if (abs(res_avg / res_mu) < 1) //sign change
	{
		double ratio = res_avg / res_mu; //see jakes thesis for terms, this has been checked
		res_mag = h_mu*abs(dir_coor*0.5*res_mu - h_mu*ratio*ratio*res_avg / 12. + dir_coor*0.5*ratio*ratio*res_mu
			+ 0.25*h_mu*res_avg);
		if (abs(ratio) > 1.0e+10)
		{
			std::cerr << "Really high ratio, may cause problems in ResidualSourceFaceCalculation" << std::endl;
			system("pause");
			exit(1);
		}
	}
	else //no sign change
	{
		res_mag = h_mu*abs(dir_coor*res_avg + h_mu*res_mu / 6.);
	}	
}