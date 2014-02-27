#include "ResidualSource.h"
#include "Particle1D.h"
#include <cmath>

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
	std::vector<double> zeros(2, 0.0);

	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren()) //children are in elements vector, so no need to add explicitly here
		{
			//currently zeros are stored for the abstract elements, using a map and 
			_residual_element_LD_values[(*it_el)->getID()] = zeros;
			_residual_face_LD_values[(*it_el)->getID()] = zeros;
			res_face_mags[(*it_el)->getID()] = 0.0;
			res_element_mags[(*it_el)->getID()] = 0.0;
			continue; 
		}

		//reset magnitudes just in case
		res_face_mag_el = 0.0;
		res_element_mag_el = 0.0;

		//compute volume integrals, for cells that have no children
		computeElementResidual(*it_el, res_LD_values_el, res_element_mag_el);
		_residual_element_LD_values[(*it_el)->getID()] = res_LD_values_el;
		res_element_mags[(*it_el)->getID()] = res_element_mag_el;

		if ((*it_el)->getDownStreamElement() != NULL) //in this case there is a down wind cell to add value to the downwind cell
		{
			if ((*it_el)->getDownStreamElement()->hasChildren())
			{
				//need to compute a face residual term for both down stream elements
			}
			else
			{
				computeFaceResidual(*it_el, res_LD_values_face, res_face_mag_el, false); //for non boundary cells
				_residual_face_LD_values[(*it_el)->getDownStreamElement()->getID()] = res_LD_values_face;
				res_face_mags[(*it_el)->getDownStreamElement()->getID()] = res_face_mag_el;
			}
		}

		//add magnitudes to appropriate values,
		_vol_src_total += res_element_mag_el; //units of p / sec
		_face_src_total += res_face_mag_el; //units of p / sec
	}

	//compute the boundary condition LD representation over the face of each boundary element
	computeBCAngularFluxDof();

	int bc_element_ID;
	for (int i = 0; i < boundary_cells.size(); i++)
	{
		//add terms to the elements on the boundary (do not add to downwind element)
		res_face_mag_el = 0.0; //initialize
		bc_element_ID = boundary_cells[i];
		computeFaceResidual((*elements)[bc_element_ID], res_LD_values_face, res_face_mag_el, true);

		//add values to arrays
		_residual_face_LD_values[bc_element_ID] = res_LD_values_face;
		res_face_mags[bc_element_ID] = res_face_mag_el;
		_face_src_total += res_face_mag_el; //units of p / sec
	}

	if (_sampling_method == HoMethods::STANDARD_SAMPLING) //standard source sampling
	{
		//Create sampler with alias sampling, let it normalize, delete unneccessary data
		_element_source = new AliasSampler(res_element_mags, false);
		_face_source = new AliasSampler(res_face_mags, false);
	}
	else
	{
		std::cerr << "No other sampling methods are implemented yet" << std::endl;
		exit(1);
	}

	//No BC source, it is implicitly included in the face residual calculations
}

ResidualSource::~ResidualSource()
{
	delete _element_source;
	delete _face_source;
}

void ResidualSource::sampleSourceParticle()
{
	//Sample if face or element source
	if (_rng->rand_num() < (_vol_src_total / (_vol_src_total + _face_src_total))) //element source
	{
		_particle->_current_element_ID = _element_source->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties();

		//determine location and direction within element using rejection method
		sampleElementSource();
	}
	else //face source
	{
		_particle->_current_element_ID = _face_source->sampleBin(_rng->rand_num(), _rng->rand_num()); //sample bin location
		_particle->_current_element = _particle->_mesh->getElement(_particle->_current_element_ID); //update bin
		_particle->updateElementProperties();

		//determine location and direction within element using rejection method
		sampleFaceSource();
	}
}

void ResidualSource::sampleFaceSource()
{
	ECMCElement1D* element = _particle->_current_element;
	//put particle on the upwind face
	if (element->getAngularCoordinate() > 0.0) //particle moving to the right
	{
		_particle->_position_mfp = 0.0; 
	}
	else //moving to the left, on right face
	{
		_particle->_position_mfp = _particle->_element_width_mfp;
	}

	//Sample angle using rejection method
	double h_mu = element->getAngularWidth();
	double dir_coor = element->getAngularCoordinate();  //mu_m, direction coordinate
	std::vector<double> res_dof = _residual_face_LD_values[element->getID()]; //res_x is 0
	double & res_avg = res_dof[0];
	double & res_mu = res_dof[2];	
	if (res_dof.size() != 3)
	{
		std::cerr << "Evaluating Lin Disc Function that is not 2D, in ResidualSource.cpp\n";
		system("pause");
		exit(1);
	}
	
	//Determine max value of residual pdf function "f"
	double mu_max = 0.5*dir_coor - 0.25*res_avg*h_mu / res_mu;
	double f_minus = std::abs((dir_coor - 0.5*h_mu)*(res_avg - res_mu)); //value on left boundary
	double f_plus = std::abs((dir_coor + 0.5*h_mu)*(res_avg + res_mu)); // value on right boundary
	double f_max = std::fmax(f_minus, f_plus);
	if (std::abs(mu_max - dir_coor) < 0.5*h_mu) //max of quadratic is within range
	{
		f_max = std::fmax(f_max, std::abs( mu_max*(evalLinDiscFunc2D(res_dof,element,0.0,mu_max))));
	}

	//Perform rejection sampling for direction
	double mu_new = 0;
	double f_new;
	while (true)
	{
		mu_new = dir_coor + h_mu*(_rng->rand_num() - 0.5); //pick new direction
		f_new = (mu_new*evalLinDiscFunc2D(res_dof, element, 0.0, mu_new)); //evaluate residual at new mu TODO abs of mu is artifact, not sure why yet, jakes code has this to force to work, probably because deltas are defined going in the wrong direction
		if (_rng->rand_num()*f_max < std::abs(f_new)) //keep mu
		{
			_particle->_mu = mu_new;
			if (f_new < 0.0) //is residual negative here?
			{
				_particle->_weight *= -1.0;
			}
			break;
		}
	}
}

void ResidualSource::sampleElementSource()
{
	ECMCElement1D* element = _particle->_current_element;
	
	//local variables
	double h_mu = element->getAngularWidth();
	double h_x = element->getSpatialWidth();
	double dir_coor = element->getAngularCoordinate();  //mu_m, direction coordinate, cell center mu
	double pos_coor = element->getSpatialCoordinate();  //x_i, cell center
	std::vector<double> res_dof = _residual_element_LD_values[element->getID()]; //res_x is 0
	double & res_avg = res_dof[0];
	double & res_x = res_dof[1];
	double & res_mu = res_dof[2];

	//determine maximum of function over the element (one of the nodes since it is linear)
	double f_max = std::abs(res_avg) + std::abs(res_x) + std::abs(res_mu);

	//perform rejection sampling
	double f_new, mu_new, x_new; //sampled variables
	while (true)
	{
		x_new = pos_coor + h_x*(_rng->rand_num() - 0.5); 
		mu_new = dir_coor + h_mu*(_rng->rand_num() - 0.5); 
		f_new = evalLinDiscFunc2D(res_dof, element, x_new, mu_new);

		if ((_rng->rand_num()*f_max) < abs(f_new)) //keep sample
		{
			_particle->_mu = mu_new;
			_particle->_position_mfp = ((x_new - pos_coor) / h_x + 0.5)
				* _particle->_element_width_mfp; //normalized position in element
			
			if (f_new < 0.0) //residual was negative set weight
			{
				_particle->_weight *= -1.0;
			}
			break;
		}
	}
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


	//this first block is for the special cases to eliminate divide by zero errors
	if (std::abs(res_x/res_avg) < GlobalConstants::RELATIVE_TOLERANCE) //x_slope~0
	{
		if (std::abs(res_mu)/std::abs(res_avg) < GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0 also
		{
			res_mag = std::abs(res_avg);
		}
		else //mu slope non-zero
		{			
			res_mag = (res_mu*res_mu + res_avg*res_avg) / std::abs(2.*res_mu);
		}
	}
	else if (std::abs(res_mu)/std::abs(res_avg)<GlobalConstants::RELATIVE_TOLERANCE) //mu_slope~0
	{
		res_mag = (res_x*res_x + res_avg*res_avg) / std::abs(2.*res_x);
	}
	//---
	else if (r_left_plus*r_right_plus > 0.0) //no sign change on top
	{
		if (r_left_minus*r_right_minus > 0.0) //no change on bottom
		{
			if (r_left_plus*r_left_minus > 0.0) //no change at all
			{
				res_mag = std::abs(res_LD_values_el[0]);
			}
			else if (r_left_plus*r_right_minus < 0.0) //change on left and right
			{
				res_mag = (res_avg*res_avg + res_x*res_x / 3. + res_mu*res_mu) / std::abs(2.*res_mu);
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
				res_mag = std::abs(res_avg + pow(r_right_minus, 3) / (12.*res_x*res_mu));
			}
			else if (r_right_plus*r_left_minus < 0.0) //change on bottom and left
			{
				res_mag = std::abs(res_avg + pow(r_left_minus, 3) / (12.*res_x*res_mu));
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
				res_mag = std::abs(res_avg + pow(r_right_plus, 3) / (12.*res_x*res_mu));
			}
			else if (r_left_plus*r_left_minus < 0.0) //change on top and left
			{
				res_mag = std::abs(res_avg + pow(r_left_plus, 3) / (12.*res_x*res_mu));
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
			res_mag = (res_avg*res_avg + res_x*res_x + res_mu*res_mu / 3.) / (2 * std::abs(res_x));
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
	double mu_sgn = dir_coor / std::abs(dir_coor); //negative or positive direction
	double & res_avg = res_LD_values_face[0];
	double & res_mu = res_LD_values_face[2]; 

	std::vector<double> psi_down; //downstream values
	std::vector<double> psi_up; //upstream element values

	if (on_boundary) //for non isotropic boundary conditions this will not be true
	{
		psi_up = _bc_dof[_bc_element_to_dof_map[element->getID()]];
		psi_down = element->getAngularFluxDOF();
	}
	else
	{	
		ECMCElement1D* ds_element = element->getDownStreamElement();
		if (ds_element->hasChildren()) //downstream element is more refined
		{


		}
		else if (ds_element->getRefinementLevel() == element->getRefinementLevel() - 1) //upstream element is less refined
		{
			std::vector<double> psi_child_minus = element->getChild(1)->getAngularFluxDOF(); //bottom child
			std::vector<double> psi_child_plus = element->getChild(3)->getAngularFluxDOF(); //top child

			//compute upstream values by averaging the upstream terms
			psi_up[0] = 0.5*(psi_child_minus[0] + psi_child_plus[0]);
			psi_up[1] = 0.5*(psi_child_minus[1] + psi_child_plus[1]);
			psi_up[2] = 0.4;

			psi_down = element->getDownStreamElement()->getAngularFluxDOF(); //downstream values
		}
		else if (ds_element->getRefinementLevel() == element->getRefinementLevel()) //same refinement
		{
			psi_up = element->getAngularFluxDOF(); //upstream element values	
			psi_down = element->getDownStreamElement()->getAngularFluxDOF(); //downstream values
		}
		else
		{
			std::cerr << "Refinement is not regularized, one cell is more than 2 levels more refined than another, in ResidualSource::computeFaceResidual\n";
			exit(1);
		}
	}

	//mu_sgn*psi[1] tells you if on left or right face
	//res_X is zero here, note that the residual here is -1.*dpsi/dx, and does not include mu term
	res_avg = mu_sgn*((psi_up[0] + mu_sgn*psi_up[1]) - (psi_down[0] - mu_sgn*psi_down[1]));
	res_mu = mu_sgn*(psi_up[2] - psi_down[2]);
	
	//compute magnitude of integral
	if (std::abs(res_mu / res_avg) < GlobalConstants::RELATIVE_TOLERANCE)
	{
		res_mag = h_mu*std::abs(res_avg*dir_coor);
	}
	if (std::abs(res_avg / res_mu) < 1) //sign change
	{
		double ratio = res_avg / res_mu; //see Jake Peterson A&M thesis for terms, computation has been checked
		res_mag = h_mu*std::abs(dir_coor*0.5*res_mu - h_mu*ratio*ratio*res_avg / 12. + dir_coor*0.5*ratio*ratio*res_mu
			+ 0.25*h_mu*res_avg);
	}
	else //no sign change
	{
		res_mag = h_mu*std::abs(dir_coor*res_avg + h_mu*res_mu / 6.);
	}	
}

void ResidualSource::computeBCAngularFluxDof()
{
	//Calculate the magnitude of the boundary condition source over elements
	std::vector<DirichletBC1D*> bcs = _particle->_mesh->getDirichletBCs();
	std::vector<int> bc_elem_IDs = _particle->_mesh->findUpwindBoundaryCells();
	std::vector<double> bc_dof_el; //LD dof for each bc, for now isotropic so will only be the average term
	ECMCElement1D* element;
	double incident_current, two_pi_incident_flux, direction;
	if (bcs.size() > 2)
	{
		std::cout << "Currently only implemented BCs for 1D problems, in ResidualSource.cpp\n";
	}

	for (int i_bc = 0; i_bc < bcs.size(); i_bc++)
	{
		incident_current = bcs[i_bc]->getValue();	//p/sec entering the spatial cell
		two_pi_incident_flux = 2.*incident_current; //This assumes incident flux is isotropic in halfspace and ECMC elements are azimuthally integrated

		//map flux to boundary elements
		for (int id = 0; id < bc_elem_IDs.size(); id++)
		{
			int el_id = bc_elem_IDs[id]; //id of current boundary element
			element = _particle->_mesh->getElement(el_id);
			if (element->getSpatialElement() != bcs[i_bc]->getElement()) //on the wrong face
			{
				continue;
			}
			else
			{
				bc_dof_el.assign(3, 0.0);
				bc_dof_el[0] = two_pi_incident_flux; //fraction in this bin
				_bc_dof.push_back(bc_dof_el);
				_bc_element_to_dof_map[el_id] = _bc_dof.size() - 1; //index in bc_dof the last element corresponds to
			}
		}
	}
}

double ResidualSource::getTotalSourceStrength()
{
	return _face_src_total + _vol_src_total;
}