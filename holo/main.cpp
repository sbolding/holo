#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>
#include "RNG.h"
#include "DataTransfer.h"
#include <cmath>

int main()
{
	using std::cout;
	using std::endl;

	cout << "We're rolling..." << endl;
	LoSolver* lo_solver;
	HoSolver* ho_solver;

	//open output file
	ofstream out_file("Z:/TAMU_Research/HOLO/results_output_folder/results.txt", std::ofstream::out);
	if (!out_file)
	{
		std::cerr << "Can't open file" << endl;
		exit(1);
	}

	//Temporarily hard coded dimensions until there is stuff for reading from input file
	int dimension = 1;
	double width = 3; //cm
	double ext_source = 1.0; //(p/(sec cm^3)), do not use non-zero values << 1, or some logic may be wrong currently
	int num_elems = 5;
	int n_ang_elements = 2; //number angles in half ranges
	//Temporarily hard coded monte carlo parameters
	int n_histories = num_elems*2*n_ang_elements*120; //50000000
	int n_batches = 100;
	double exp_convg_rate = 0.05;
	double convergence_tolerance = 1.E-3;
	string solver_mode = "holo-ecmc"; //"standard-mc", "holo-ecmc", "holo-standard-mc"
	string sampling_method = "stratified";
					  // ID, sig_a, sig_s
	MaterialConstant mat(10, 0.25, 0.0);

	//Create the mesh and elements;
	Mesh mesh_1D(dimension, num_elems, width, &mat);
	mesh_1D.setExternalSource(ext_source);
	mesh_1D.print(cout);

	size_t n_holo_solves = 100;

	//Assemble the LoSystem
	lo_solver = new LoSolver1D(&mesh_1D); //uses default some estimated lo order parameters and LD

	size_t i_holo_solves = 0; //counter

	//Variables for checking convergence
	std::vector<double> old_flux_vector(num_elems*(dimension+1),0.0);
	std::vector<double> new_flux_vector(num_elems*(dimension+1));

	while (true)
	{
		//solve lo order system
/*		lo_solver->solveSystem();
		lo_solver->updateSystem(); //Update lo order system scalar flux values to current solution
		mesh_1D.getDiscScalarFluxVector(new_flux_vector); */
		
		//Print LO scalar flux estimate
		mesh_1D.printLDScalarFluxValues(cout);

		if (i_holo_solves == n_holo_solves) //do one extra LO solve, because it is the only solution really being calculated
		{
			mesh_1D.printLDScalarFluxValues(out_file);
			mesh_1D.printLDScalarFluxValues(std::cout);
			break;
		}

		//Solve high order system
		ho_solver = new HoSolver(&mesh_1D, n_histories, n_ang_elements, solver_mode, sampling_method, exp_convg_rate, n_batches,2);
		ho_solver->solveSystem();
		ho_solver->updateSystem();
		ho_solver->printProjectedScalarFlux(std::cout);

		//Transfer HO estimated parameters to the LO system
		DataTransfer data_transfer(ho_solver, &mesh_1D);
		std::cout << "Change in LoData (before and after): \n";
		data_transfer.printLoDataEl(0, std::cout);
		data_transfer.updateLoSystem();
		data_transfer.printLoDataEl(0, std::cout);

		//update counter
		i_holo_solves++;

		//Check convergence of solution
		double diff_sum_sq = 0.;
		double old_sum_sq = 0.;
		for (int i_vec = 0; i_vec < new_flux_vector.size(); i_vec++)
		{
			double diff = new_flux_vector[i_vec] - old_flux_vector[i_vec];
			diff_sum_sq += diff*diff;
			old_sum_sq += old_flux_vector[i_vec] * old_flux_vector[i_vec];
		}
		std::cout << "Convergence Norm of Flux Vector: " << std::sqrt(diff_sum_sq) / std::sqrt(old_sum_sq);
		if (std::sqrt(diff_sum_sq) < convergence_tolerance*std::sqrt(old_sum_sq))
		{
			std::cout << "\nConverged on iteration " << i_holo_solves << " to a relative precision"
				<< " of " << convergence_tolerance << std::endl;
			lo_solver->solveSystem();
			lo_solver->updateSystem();
			mesh_1D.printLDScalarFluxValues(out_file);
			mesh_1D.printLDScalarFluxValues(std::cout);
		//	ho_solver->printProjectedScalarFlux(out_file);
			break;
		}
		else
		{
			old_flux_vector = new_flux_vector;
		}
	}




	return 0;
}