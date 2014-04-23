//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : HoSolver.cpp
//  @ Date : 1/27/2014
//  @ Author : 
//
//


#include "HoSolver.h"
#include "GlobalConstants.h"
#include "FEMUtilities.h"
#include <limits>
#include <cmath>

HoSolver::HoSolver()
{
	std::cerr << "You have called HoSolver Default constructor, should never happen" << std::endl;
	exit(1);
}

HoSolver::HoSolver(Mesh* mesh, int n_histories,
	int n_ang_bins_half_range, string method, string sampling_method,
	double exp_convg_constant, int n_batches, int n_batches_to_avg) :
	_rng(), _source(), _psi_minus_dof(), _psi_plus_dof()
{
	_n_batches = n_batches;
	_lo_mesh = mesh;
	_n_histories = n_histories;
	_solver_mode_str = method;
	_solver_mode_int = HoMethods::method_map.at(method);

	//Create high order mesh
	_ho_mesh = new HoMesh(_lo_mesh, n_ang_bins_half_range);

	//initialize particle class with pointer to NULL source
	if (HoController::CONT_WGT_DEPOSITION_PARTICLES)
	{
		_particle = new CWDParticle1D(_ho_mesh, NULL, &_rng, method);
	}
	else
	{
		_particle = new Particle1D(_ho_mesh, NULL, &_rng, method);
	}

	//Initialize the data needed for source sampling
	_sampling_method = sampling_method;
	_sampling_method_index = HoMethods::sampling_map.at(sampling_method); //unsigned int that can be compared against global constants

	//create source based on full particle class, the source constructor will update the particle class pointer
	initializeSamplingSource();

	//create adaptive refinement class
	_mesh_controller = new MeshController(_ho_mesh,exp_convg_constant,n_batches_to_avg);
}

void HoSolver::solveSystem(std::ostream & out)
{
	std::cout << "Solving the HO system..." << std::endl;
	double psi_l2_error; //L2 norm of addiditive error from each batch
	double psi_ref_l2_norm; //L2 norm of psi from the first batch

	//print out initial mesh to mesh file
	if (HoController::WRITE_MESH_EVERY_REFINEMENT)
	{
		//open output file
		ofstream mesh_file("Z:/TAMU_Research/HOLO/results_output_folder/mesh.out", ios::out);
		if (!mesh_file)
		{
			std::cerr << "Can't open mesh output file" << endl;
			exit(1);
		}
		_ho_mesh->printActiveMesh(mesh_file);
		mesh_file.close();
	}

	for (int batch = 0; batch < _n_batches; batch++)
	{
		if (HoController::WRITE_BATCHES_COMPLETE) std::cout << "Running batch " << batch + 1 << " of " << _n_batches << std::endl;	
		for (int hist = 0; hist < _n_histories; hist++) //loop over the number of histories for each batch
		{
			_particle->runHistory();
			if (HoController::WRITE_HISTORIES_COMPLETE)
			{
				if ((hist + 1) % (_n_histories / 10) == 0)
				{
					std::cout << (int)((hist + 1) / (float)_n_histories * 100) << "% of " <<
						_n_histories << " histories complete..." << std::endl;
				}
			}
		}
		//compute the new angular fluxes, determine the L2 norm of the error
		_ho_mesh->computeAngularFluxes(_n_histories, psi_l2_error, _source->getTotalSourceStrength());

		//Convergence criteria for angular flux
		if (batch == 0)
		{
			psi_ref_l2_norm = psi_l2_error;
		}
		else
		{
			//check convergence
			if (psi_l2_error / psi_ref_l2_norm < HoConstants::ECMC_REL_ERR_TOL)
			{
				out << "Error converged on batch " << batch + 1 << " to a rel. tol. of "
					<< HoConstants::ECMC_REL_ERR_TOL << std::endl;
				break;
			}
		}
	
		//debug outputs
		if (HoController::WRITE_RELATIVE_ERROR_NORMS)
		{
			out.setf(ios::scientific);
			out.precision(15);
			out << "Relative Error L2 Norm: " << psi_l2_error / psi_ref_l2_norm << std::endl;

		}
		if (HoController::PARTICLE_BALANCE)  _particle->printParticleBalance(_n_histories);
		if (HoController::WRITE_ALL_ANGULAR_FLUXES) _ho_mesh->printAngularFluxes(out);

		//compute residual or exit
		if (_solver_mode_int == HoMethods::STANDARD_MC) //you are done, exit
		{
			return;
		}
		else //ECMC method
		{
			//compute new residual source
			computeResidualSource();
		}

		//store and print residual norm
		double resid_L1_norm = _source->getTotalSourceStrength();
		_mesh_controller->storeResidualNorm(resid_L1_norm);
		if (HoController::WRITE_RESIDUAL_NORMS && HoController::WRITE_BATCHES_COMPLETE)
		{
			out.setf(ios::scientific);
			out.precision(15);
			out << "Residual L1 Norm: " << resid_L1_norm << std::endl;
		}
		else if (HoController::WRITE_BATCHES_COMPLETE)
		{
			out << std::endl;
		}

		//if necessary refine solution
		if (HoController::ADAPTIVE_REFINEMENT && batch < (_n_batches -1) )
		{
			if (_mesh_controller->meshNeedsRefinement())
			{
				if (HoController::WRITE_BATCHES_COMPLETE) std::cout << "Refining mesh...\n"; //Debug output
				int n_elems_before_refinement = _ho_mesh->getNumActiveElements();
				int n_histories_before_refinement = _n_histories;
				_mesh_controller->refineMesh(); //this will refine if necessary
				_n_histories = (int)((double)_n_histories*((double)_ho_mesh->getNumActiveElements() / (double)n_elems_before_refinement)); //update number of histories before computing new residual, casting is necessary, trust me
				computeResidualSource(); //need to recompute residual for the new cells, if refinement occured

				//Error check for int overflow
				if (_n_histories < n_histories_before_refinement)
				{
					std::cerr << "Int overflow in num histories in HoSolver, need to store n_hist_per_elem instead of n_histories\n";
					exit(1);
				}

				//output the mesh if desired
				if (HoController::WRITE_MESH_EVERY_REFINEMENT)
				{
					std::cout << "New number histories: " << _n_histories << std::endl; //DEBUG	
					std::cout << "New mesh has " << _ho_mesh->getNumActiveElements() << " active elements\n";
					//open output file
					ofstream mesh_file("Z:/TAMU_Research/HOLO/results_output_folder/mesh.out", ios::app);
					if (!mesh_file)
					{
						std::cerr << "Can't open mesh output file" << endl;
						exit(1);
					}
					_ho_mesh->printActiveMesh(mesh_file);
					mesh_file.close();
				}
			}
		}
	}	
}

void HoSolver::getLoData1D(LoData1D & lo_data, int element_id)
{
	//NOTE: all variables in this section are independent of the sources strength in the problem, 
	//not in general true

	//local variables to calculate	
	AveragedCosines surf_cosines;
	AveragedCosines vol_cosines;
	
	//Because HO solution is projected LD solution in half range, the spatial factor is inheriently 2.0 (default LD value)
	lo_data.setSpatialClosureFactor(2.0);

	//Get the angular flux dof for this element
	const std::vector<double> psi_plus_el(_psi_plus_dof[element_id]);
	const std::vector<double> psi_minus_el(_psi_minus_dof[element_id]);

	//-------------------------
	//Calculate face values
	//-------------------------
	double edge_flux_plus = psi_plus_el[0] + psi_plus_el[1]; //int_0,1 psi(x_R,mu) dmu
	double edge_flux_minus = psi_minus_el[0] - psi_minus_el[1]; //int_0,1 psi(x_L,mu) dmu

	if (edge_flux_plus > std::fmax(psi_plus_el[0],psi_plus_el[1]) * GlobalConstants::RELATIVE_TOLERANCE &&
		edge_flux_minus > std::fmax(psi_plus_el[0], psi_plus_el[1])* psi_minus_el[0] * GlobalConstants::RELATIVE_TOLERANCE)
	{
		surf_cosines._mu_right_plus = 0.5*(edge_flux_plus + (psi_plus_el[2]+psi_plus_el[3]) / 3.) / edge_flux_plus; //added cross moment
		surf_cosines._mu_left_minus = -0.5*(edge_flux_minus - (psi_minus_el[2] - psi_minus_el[3]) / 3.) / edge_flux_minus; //negative 0.5 is the mu_center of -1 to 0 angular element
	}
	else
	{
		std::cerr << "Have a zero flux, not sure what to do about it yet, in HoSolver::getLoData1D";
		exit(1);
	}

	//set others to infinity to ensure that they are not being used in algorithm (because of upwinding they shouldnt ever be used)
	surf_cosines._mu_right_minus = std::numeric_limits<double>::infinity();
	surf_cosines._mu_left_plus = std::numeric_limits<double>::infinity();

	//-----------------------
	//Calculate Vol Values
	//-----------------------
		
	//Convert the avg, slope dof to the left and right moments, since psi_mu terms are not effected by basis integrals because of how average is defined
	std::vector<double> basis_moments_plus(2), basis_moments_minus(2);
	std::vector<double> avg_slope_plus(psi_plus_el.begin(), psi_plus_el.begin()+2); //does not include mu and x_mu moments
	std::vector<double> avg_slope_minus(psi_minus_el.begin(), psi_minus_el.begin()+2); //does not include mu and x_mu moments
	FEMUtilities::convertAvgSlopeToBasisMoments1D(avg_slope_plus, basis_moments_plus); //convert to spatial moments
	FEMUtilities::convertAvgSlopeToBasisMoments1D(avg_slope_minus, basis_moments_minus);


	//check that there is non zero moments
	if (basis_moments_plus[0] < GlobalConstants::RELATIVE_TOLERANCE ||
		basis_moments_minus[0] < GlobalConstants::RELATIVE_TOLERANCE )
	{
		std::cerr << "Have a zero flux moment, not sure what to do about it yet, in HoSolver::getLoData1D, probalby just assume diffusions parameters";
		exit(1);
	}

	//Calculate average mu's based on basis moments, very straightforward, use mu moment value from psi_plus and minus
	vol_cosines._mu_left_minus = -0.5*(basis_moments_minus[0] - (psi_minus_el[2] - psi_minus_el[3]) / 3.) / basis_moments_minus[0]; //left moment, minus: <.>L^-
	vol_cosines._mu_right_minus = -0.5*(basis_moments_minus[1] - (psi_minus_el[2] + psi_minus_el[3]) / 3.) / basis_moments_minus[1]; //etc.
	vol_cosines._mu_left_plus = 0.5*(basis_moments_plus[0] + (psi_plus_el[2] - psi_plus_el[3]) / 3.) / basis_moments_plus[0];
	vol_cosines._mu_right_plus = 0.5*(basis_moments_plus[1] + (psi_plus_el[2] + psi_plus_el[3]) / 3.) / basis_moments_plus[1];

	//Set the LoData cosine values
	lo_data.setSurfAveragedCos(surf_cosines);
	lo_data.setVolAveragedCos(vol_cosines);

	//If desired print out the low order half range flux dof's
	if (HoController::WRITE_HALF_RANGE_LO_FLUXES)
	{
		std::cout << "Element ID = " << element_id;
		std::cout << "\n\n Plus DOF: ";
		for (int i = 0; i < psi_plus_el.size(); ++i)
		{
			std::cout << psi_plus_el[i] << "  ";
		}
		std::cout << "\n Minus DOF: ";
		for (int i = 0; i < psi_minus_el.size(); ++i)
		{
			std::cout << psi_minus_el[i] << "  ";
		}
		std::cout << "\n\n";
	}

}

void HoSolver::printAllTallies(std::ostream& out) const
{
	using std::endl;

	out << "\n---------------------------------------------------------\n"
		<< "                   Element Angular Fluxes\n"
		<< "---------------------------------------------------------\n";
	//Print out the homesh tallies
	_ho_mesh->printAngularFluxes(out);
}

void HoSolver::updateSystem()
{
	//calculate projected fluxes
	computeProjectedAngularFlux();
}

void HoSolver::computeResidualSource()
{
	//Free memory for old source (whether old residual or standard MC source)
	delete _source;

	//create new source
	if (_sampling_method_index == HoMethods::STANDARD_SAMPLING)
	{
		_source = new StandardResidualSource(_particle);
	}
	else if (_sampling_method_index == HoMethods::STRATIFIED_SAMPLING)
	{
		//This function will increase the number of histories slightly to make sure that an even number of particles can be sampled from each element
		_source = new StratifiedResidualSource(_particle, _n_histories);
	}
	else
	{
		std::cerr << "Sampling method not implemented in HoSolver::computeResidualSource()\n";
		exit(1);
	}
}

void HoSolver::initializeSamplingSource()
{
	//Initially source is always a standard mc source of some kind
	//_source = new LinDiscSource(_particle); //uses standard sampling, no stratified available for LinDiscSource
	//_source = new StratifiedResidualSource(_particle, _n_histories); //could just use residual source since initially residual is just the ext_source lin_disc, but I use LinDiscSource for debugging and sanity check}
	_source = new StandardResidualSource(_particle);
}

void HoSolver::computeProjectedAngularFlux()
{
	//initialize vectors to proper size
	int n_spatial_elems = _lo_mesh->getNumElems();
	_psi_minus_dof.resize(n_spatial_elems);
	_psi_plus_dof.resize(n_spatial_elems);

	//resize dof vectors, initialize to zeros
	size_t dof_size = _ho_mesh->getElement(0)->getElementDimensions().size()+2; //how many DOF per element, 2 extra for bilinear finite element
	size_t n_dimens = dof_size - 2; //how many dimensions, spatial and angular
	for (int i = 0; i < n_spatial_elems; i++)
	{
		_psi_minus_dof[i].assign(dof_size, 0.0);
		_psi_plus_dof[i].assign(dof_size, 0.0);
	}

	//loop over ECMC elements and add active values to teh appropriate DOF
	std::vector<ECMCElement1D*>::iterator it_el;
	std::vector<ECMCElement1D*>* elements = _ho_mesh->getElements();

	//loop variables
	int spatial_id;
	GaussQuadrature quad(4);
	int n_qps = quad.getNumPoints();
	std::vector<double> x_wgts, mu_wgts, x_pnts, mu_pnts;
	std::vector<double> sp_nodes; //nodeal cordinates of spatial element
	std::vector<double> sp_coors(n_dimens); //the coords of the half range flux element
	std::vector<double> sp_dimens(n_dimens); //the dimensions of the half range flux element
	std::vector<double> el_coors;	//cooridnates of current ECMC element
	std::vector<double> el_dimens;  //dimensions of current ECMC element;
	std::vector<double> el_dof; //angular flux dof for current element
	
	//the half range fluxes always have the same angular width;
	sp_dimens[1] = 1.0;

	//To perform integral over half range of spatial elements, sum integrals over all active elements for a particular spatial element.  
	//This is done by looping over all high order elements and adding to the appropriate spatial elements integrals
	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren())
		{
			continue; //only add values from active elements
		}
		else
		{
			//get the spatial element information
			Element* sp_elem = (*it_el)->getSpatialElement();
			spatial_id = sp_elem->getID();
			sp_nodes = sp_elem->getNodalCoordinates();
			sp_dimens[0] = (sp_nodes[1] - sp_nodes[0]);
			sp_coors[0] = (sp_nodes[1] + sp_nodes[0])*0.5;
			sp_coors[1] = ((*it_el)->getAngularCoordinate() > 0.0) ? 0.5 : -0.5;

			//get ECMC element information
			el_coors = (*it_el)->getElementCoordinates();
			el_dimens = (*it_el)->getElementDimensions();
			el_dof = (*it_el)->getAngularFluxDOF();

			//get quadrature information based on ECMC element, which you are integrating over
			x_pnts = quad.getQuadraturePoints((*it_el)->getSpatialCoordinate(), (*it_el)->getSpatialWidth());
			x_wgts = quad.getQuadratureWeights((*it_el)->getSpatialCoordinate(), (*it_el)->getSpatialWidth());
			mu_pnts = quad.getQuadraturePoints((*it_el)->getAngularCoordinate(), (*it_el)->getAngularWidth());
			mu_wgts = quad.getQuadratureWeights((*it_el)->getAngularCoordinate(), (*it_el)->getAngularWidth());
			
			std::vector<double> sp_sum(dof_size, 0.0); //temp vector for storing calculations
			for (int i_qp=0; i_qp < n_qps; ++i_qp) //x_qps
			{
				for (int j_qp = 0; j_qp < n_qps; ++j_qp) //mu_qps
				{
					//1/(h_x*h_mu)*psi(x,mu)), where hx and hmu are for the half range element
					double psi_avg_qp= x_wgts[i_qp] * mu_wgts[j_qp] / (sp_dimens[0] * sp_dimens[1])*
						FEMUtilities::evalLinDiscFunc2D(el_dof, el_dimens, el_coors, x_pnts[i_qp], mu_pnts[j_qp]); //compute contribution to average flux moment
					
					//compute the integrals for each basis function
					sp_sum[0] += psi_avg_qp;
					sp_sum[1] += 6.0*psi_avg_qp*(x_pnts[i_qp] - sp_coors[0]) / sp_dimens[0]; //these two lines could be replaced by a for loop if x,mu pnts in one vector
					sp_sum[2] += 6.0*psi_avg_qp*(mu_pnts[j_qp] - sp_coors[1]) / sp_dimens[1]; //basis function is based on unrefined half range element, because that is what we are projecting to
					sp_sum[3] += 36.0*(psi_avg_qp*(mu_pnts[j_qp] - sp_coors[1]) / sp_dimens[1]) * ((x_pnts[i_qp] - sp_coors[0]) / sp_dimens[0]);
				}
			}

			//add sp_sum values to the appropriate spatial element
			for (int i_dof=0; i_dof < dof_size; i_dof++)
			{
				if (el_coors[1] > 0.0)
				{
					_psi_plus_dof[spatial_id][i_dof] += sp_sum[i_dof];
				}
				else
				{
					_psi_minus_dof[spatial_id][i_dof] += sp_sum[i_dof];
				}
			}
		}
	}
}

std::vector<double> HoSolver::getScalarFluxDOF(int spatial_elem_id) const
{
	std::vector<double> scalar_flux_dof(_lo_mesh->getSpatialDimension()+1);  
	for (int i = 0; i < scalar_flux_dof.size(); ++i) //avg, spatial moments, then angular moments which are neglected
	{
		scalar_flux_dof[i] = _psi_minus_dof[spatial_elem_id][i] + _psi_plus_dof[spatial_elem_id][i];
	}
	return scalar_flux_dof;
}

void HoSolver::printProjectedScalarFlux(std::ostream &out) const
{
	out << "------------------------------------------------------------\n"
		<< "			Element Scalar Flux Edge Values			\n"
		<< " (Node, flux value)     \n"
		<< "------------------------------------------------------------\n";
	std::vector<double> flux_dof;
	std::vector<double> flux_edge_values;
	std::vector<double> nodal_coords;
	for (int id = 0; id < _lo_mesh->getNumElems(); id++)
	{
		nodal_coords = _lo_mesh->getElement(id)->getNodalCoordinates();
		flux_dof = getScalarFluxDOF(id);	//get moments
		FEMUtilities::convertMomentsToEdgeValues1D(flux_dof, flux_edge_values); //get edge values

		for (int j = 0; j < flux_dof.size(); ++j)
		{
			out.precision(5);
			out << nodal_coords[j] << " " << setw(25) <<
				setprecision(14) << scientific << flux_edge_values[j] << endl;
		}
	}
}