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

HoSolver::HoSolver()
{
	std::cerr << "You have called HoSolver Default constructor, should never happen" << std::endl;
	exit(1);
}

HoSolver::HoSolver(Mesh* mesh, int n_histories,
	int n_ang_bins_half_range, string method, double exp_convg_constant, int n_batches, int n_batches_to_avg) :
	_rng()
{
	_n_batches = n_batches;
	_lo_mesh = mesh;
	_n_histories = n_histories;
	_current_face_tallies.resize(mesh->getNumEdges());		//initialize tallies to appropriate size
	_current_element_tallies.resize(mesh->getNumElems());  
	_flux_face_tallies.resize(mesh->getNumEdges());		//initialize tallies to appropriate size
	_flux_element_tallies.resize(mesh->getNumElems());
	_solver_mode_str = method;
	_solver_mode_int = HoMethods::method_map.at(method);

	//initialize all tallies
	for (size_t face = 0; face < _current_face_tallies.size(); ++face)
	{
		_flux_face_tallies[face] = new FluxFaceTally(2);
		_current_face_tallies[face] = new CurrentFaceTally(2);
	}
	for (size_t elem = 0; elem < _current_element_tallies.size(); ++elem)
	{
		_flux_element_tallies[elem] = new FluxElementTally(2,2);
		_current_element_tallies[elem] = new CurrentElementTally(2,2);
	}

	//Create high order mesh
	_ho_mesh = new HoMesh(_lo_mesh, n_ang_bins_half_range);

	//initialize particle class
	_particle = new Particle1D(_ho_mesh, &_rng, method, _current_face_tallies,
		_current_element_tallies, _flux_face_tallies, _flux_element_tallies);

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
		_ho_mesh->computeAngularFluxes(_n_histories, _particle->getTotalSourceStrength());

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
			_particle->computeResidualSource();
		}

		//store and print residual norm
		double resid_L1_norm = _particle->getTotalSourceStrength();
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
				_particle->computeResidualSource(); //need to recompute residual for the new cells, if refinement occured
				_n_histories = (int)(_n_histories*_ho_mesh->getNumActiveElements() / (double)n_elems_before_refinement);
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
		<< "                   Element tallies\n"
		<< "(ID) (Type) (- dir. spat. moments) (+ dir. spat. moments)\n"
		<< "---------------------------------------------------------\n";
	for (int elem = 0; elem < _lo_mesh->getNumElems(); elem++)
	{
		printElementTally(elem, out);
	}
	out << "\n---------------------------------------------------------\n"
		<< "             Face tallies\n"
		<< "(ID) (Type) (- dir. tally) (+ dir. tally)\n"
		<< "---------------------------------------------------------\n";
	for (size_t face = 0; face < _current_face_tallies.size(); face++)
	{
		printFaceTally(face, out);
	}

	//Print out the homesh tallies
	_ho_mesh->printAngularFluxes(out);
}

void HoSolver::printElementTally(int element_id, std::ostream &out) const
{
	std::vector<std::vector<double>> values;

	out << "ID = " ;
	out.width(5);
	out << element_id << " Current = ";
	out.setf(ios::scientific);
	out.precision(8);
	out.width(12);
	
	//get the values then print them
	values = _current_element_tallies[element_id]->getScores(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";
	out << "  Flux = ";

	//get the values then print them
	values = _flux_element_tallies[element_id]->getScores(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";

	out << endl;

	//print out the errors
	out << "Std. Dev.:             ";

	//get the values then print them
	values = _current_element_tallies[element_id]->getStdDevs(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";
	out << "         ";

	//get the values then print them
	values = _flux_element_tallies[element_id]->getStdDevs(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";

	out << endl;

}

void HoSolver::printFaceTally(int face_id, std::ostream &out) const
{
	std::vector<std::vector<double>> values;

	out << "ID = ";
	out.width(5);
	out << face_id << " Current = ";
	out.setf(ios::scientific);
	out.precision(8);
	out.width(12);

	//get the values then print them
	values = _current_face_tallies[face_id]->getScores(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";
	out << "  Flux = ";

	//get the values then print them
	values = _flux_face_tallies[face_id]->getScores(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";

	out << endl;

	//Now print out the errors
	out << "Std. Dev.:             ";

	//get the values then print them
	values = _current_face_tallies[face_id]->getStdDevs(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";
	out << "         ";

	//get the values then print them
	values = _flux_face_tallies[face_id]->getStdDevs(_n_histories);

	for (size_t i = 0; i < values.size(); ++i)
		for (size_t j = 0; j < values[i].size(); ++j) out << values[i][j] << " ";

	out << endl;

}

void HoSolver::updateSystem()
{
	//update angular fluxes
	_ho_mesh->computeAngularFluxes(_n_histories);
}