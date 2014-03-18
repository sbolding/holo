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

	for (int batch = 0; batch < _n_batches; batch++)
	{
		std::cout << "\nRunning batch " << batch + 1 << " of " << _n_batches;
		//loop over the number of histories
		for (int hist = 0; hist < _n_histories; hist++)
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
		//compute the new angular fluxes
		_ho_mesh->computeAngularFluxes(_n_histories, _source->getTotalSourceStrength());

		//debug outputs
		if (HoController::PARTICLE_BALANCE)  _particle->printParticleBalance(_n_histories);
		if (HoController::WRITE_ALL_ANGULAR_FLUXES) _ho_mesh->printAngularFluxes(out);

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
		if (HoController::WRITE_RESIDUAL_NORMS)
		{
			out.setf(ios::scientific);
			out.precision(15);
			out << ", Residual L1 Norm: " << resid_L1_norm;
		}

		//if necessary refine solution
		if (HoController::ADAPTIVE_REFINEMENT && batch < (_n_batches -1) )
		{
			if (_mesh_controller->meshNeedsRefinement())
			{
				std::cout << "\nRefining mesh...";
				int n_elems_before_refinement = _ho_mesh->getNumElems();
				_mesh_controller->refineMesh(); //this will refine if necessary
				_n_histories = (int)(_n_histories*_ho_mesh->getNumActiveElements() / (double)n_elems_before_refinement); //update number of histories before computing new residual
				computeResidualSource(); //need to recompute residual for the new cells, if refinement occured
			}
		}
	}	
}


LoData1D HoSolver::getLoData(int element_id)
{
	std::cerr << "This doesnt work whatsoever\n";
	exit(1);
	/*//NOTE: all variables in this section are independent of the sources strength in the problem, 
	//not in general true
	//local variables to calculate	
	LoData1D lo_data;
	double alpha;
	double alpha_plus;
	double alpha_minus;
	AveragedCosines surf_cosines;
	AveragedCosines vol_cosines;

	//Calculate alpha based on phi left and phi right moments, and face flux value
	double phi_left_moment;
	double phi_right_moment;
	double left_face_value;
	double right_face_value;
	int left_face_id = _lo_mesh->getFaceIndex(element_id, 0);
	int right_face_id = _lo_mesh->getFaceIndex(element_id, 1);

	//based on relations between phi_zeta and phi_avg, can be derived from definition moments
	phi_right_moment = 2.*_flux_element_tallies[element_id]->getScoreAngularIntegrated(_n_histories, 1);
	phi_left_moment = 2.*(_flux_element_tallies[element_id]->getScoreAngularIntegrated(_n_histories, 0)) - phi_right_moment;
	left_face_value = _flux_face_tallies[_mesh->getFaceIndex(element_id,0)]->getScoreAngularIntegrated(_n_histories);
	right_face_value = _flux_face_tallies[_mesh->getFaceIndex(element_id, 1)]->getScoreAngularIntegrated(_n_histories);

	//calculate alpha plus and alpha minus
	alpha_plus = (right_face_value - phi_left_moment) / (phi_right_moment - phi_left_moment);
	alpha_minus = (left_face_value - phi_right_moment) / (phi_left_moment - phi_right_moment);

	//Assume alpha is average of these two for now, quite possibly totally wrong though
	alpha = 0.5*(alpha_plus + alpha_minus);
	lo_data.setSpatialClosureFactor(alpha);

	//Calculate the different surface cosines
	surf_cosines._mu_left_minus = _current_face_tallies[left_face_id]->getScore(_n_histories, 0) /
		_flux_face_tallies[left_face_id]->getScore(_n_histories,0);
//	lo_data.setSurfAveragedCos(surf_cosines);
//	lo_data.setVolAveragedCos(vol_cosines); 

	return lo_data; */
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
	//update angular fluxes
	_ho_mesh->computeAngularFluxes(_n_histories);
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
	_source = new LinDiscSource(_particle); //uses standard sampling, no stratified available for LinDiscSource
	//_source = new ResidualSource(_particle); //could just use residual source since initially residual is just the ext_source lin_disc, but I use LinDiscSource for debugging and sanity check
}

void HoSolver::computeProjectedAngularFlux()
{
	//initialize vectors to proper size
	int n_spatial_elems = _lo_mesh->getNumElems();
	_psi_minus_dof.resize(n_spatial_elems);
	_psi_plus_dof.resize(n_spatial_elems);

	//resize dof vectors, initialize to zeros
	size_t dof_size = _ho_mesh->getElement(0)->getElementDimensions().size(); //how many DOF per element
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
	GaussQuadrature quad;
	std::vector<double> wgts = quad.getQuadratureWeights();

	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren())
		{
			continue; //only add values from active elements
		}
		else
		{
			//add values to correct spatial element
			spatial_id = (*it_el)->getSpatialElement()->getID();
			//compute the integrals for each basis function
			//loop over gauss points, add term for each basis function
			

		}
	}
	
}