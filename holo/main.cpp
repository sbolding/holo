#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>
#include "RNG.h"
#include "DataTransfer.h"
#include <cmath>
#include "ConstFixedSource.h"
#include "MMSFixedSource.h"
#include "numMatrixBanded.h"
#include "numVector.h" 

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
	double width = 3.0; //cm
	double sigma_a = 1.0;
	double sigma_s = 1.0;
	double ext_source = 1.0; //(p/(sec cm^3)), do not use non-zero values << 1, or some logic may be wrong currently
	double bc_left = 0.0;
	double bc_right = 0.0;
	int num_elems = 20;
	int n_ang_elements = 5; //number angles in half ranges
	//Temporarily hard coded monte carlo parameters
	int n_histories = num_elems*2*n_ang_elements*100; //50000000
	int n_batches = 100;
	double exp_convg_rate = 0.05;
	double convergence_tolerance = 10.e-4;
	string solver_mode = "holo-ecmc"; //"standard-mc", "holo-ecmc", "holo-standard-mc"
	string sampling_method = "stratified";
					  // ID, sig_a, sig_s
	MaterialConstant mat(10, sigma_a, sigma_s);

	//MMS factors
	double a = 4.0;
	double b = 2.0;
	double c = 3.0;

	//for isotropic case
//	bc_left = a / 2.;
//	bc_right = a / 2. + b*width / 2.;

	//Array for simple isotropic boundary conditions
	double* bc_values = new double[2];
	bc_values[0] = bc_left;
	bc_values[1] = bc_right;

	//Vector of incident flux LD moments for more complicated anisotropic boundary conditions
	std::vector<std::vector<double>> bc_moments;
	bc_moments = { {a+0.5*c, c/2.}, {a+b*width-0.5*c, c/2.} };

	//Create a constant external source
	MMSFixedSource q(sigma_a*a,b*sigma_a,(b+c*(sigma_a + sigma_s))); //bilinear function, average, x coeff, mu coeff, that gives matching bc's
//	ConstFixedSource q(ext_source/2.0); //constant source (p/(sec-cm^3-str))

	//Create the mesh and elements;
	Mesh mesh_1D(dimension, num_elems, width, &mat, bc_values);
	mesh_1D.setExternalSource(q);
	mesh_1D.setBoundaryConditions(bc_moments);
	mesh_1D.print(cout);

	size_t n_holo_solves = 100;

	//Assemble the LoSystem
	lo_solver = new LoSolver1D(&mesh_1D); //uses default some estimated lo order parameters and LD

	size_t i_holo_solves = 0; //counter

	//Variables for checking convergence
	std::vector<double> old_flux_vector(num_elems*(dimension+1),0.0);
	std::vector<double> new_flux_vector(num_elems*(dimension+1));
	double old_delta_phi_norm = 0.01;
	double spectral_radius; //estimate of spectral radius estimated by change in scalar flux values
	ho_solver = new HoSolver(&mesh_1D, n_histories, n_ang_elements, solver_mode, sampling_method, exp_convg_rate, n_batches, 3); //just for debugging purposes

	while (true)
	{
		//solve lo order system
		lo_solver->solveSystem();
		lo_solver->updateSystem(); //Update lo order system scalar flux values to current solution
		mesh_1D.getDiscScalarFluxVector(new_flux_vector);
		
		if (i_holo_solves == 0)
		{
			mesh_1D.printLDScalarFluxValues(out_file); //TEMPORARY DEBUG print out Mark Diffusion Solution
		}

		//Check convergence of solution 
		double diff_sum_sq = 0.;
		double old_sum_sq = 0.;
		for (int i_vec = 0; i_vec < new_flux_vector.size(); i_vec++)
		{
			double diff = new_flux_vector[i_vec] - old_flux_vector[i_vec];
			diff_sum_sq += diff*diff;
			old_sum_sq += old_flux_vector[i_vec] * old_flux_vector[i_vec];
		}
		double diff_norm = std::sqrt(diff_sum_sq);
		double relative_diff_norm = diff_norm / std::sqrt(old_sum_sq);
		spectral_radius = std::sqrt(diff_sum_sq) / old_delta_phi_norm;
		double error = std::abs(relative_diff_norm / (1. - spectral_radius));
		std::cout << "Error estimate based on spectral radius: " << error << std::endl;
		std::cout << "Spectral Radius = " << spectral_radius << " Relative Difference = " << 
			relative_diff_norm << " Abs Diff Phi = " << diff_norm << std::endl;
		if (error < convergence_tolerance)
		{
			std::cout << "\nConverged on iteration " << i_holo_solves << " to a relative precision"
				<< " of " << convergence_tolerance << std::endl;
			//		lo_solver->solveSystem();
			//		lo_solver->updateSystem();
			mesh_1D.printLDScalarFluxValues(out_file);
			mesh_1D.printLDScalarFluxValues(std::cout);
			ho_solver->printProjectedScalarFlux(out_file);
			break;
		}
		else
		{
			old_delta_phi_norm = diff_norm;
			old_flux_vector = new_flux_vector;
		}

		if (i_holo_solves == n_holo_solves) //do one extra LO solve, because it is the only solution really being calculated
		{
			mesh_1D.printLDScalarFluxValues(out_file);
			mesh_1D.printLDScalarFluxValues(std::cout);
			break;
		}

		//Solve high order system
		ho_solver = new HoSolver(&mesh_1D, n_histories, n_ang_elements, solver_mode, sampling_method, exp_convg_rate, n_batches, 5);
		ho_solver->solveSystem();
		ho_solver->updateSystem();
		ho_solver->printProjectedScalarFlux(std::cout);
		
		//Transfer HO estimated parameters to the LO system
		DataTransfer data_transfer(ho_solver, &mesh_1D);
		data_transfer.updateLoSystem();
		if (HoController::WRITE_ALL_LO_DATA)
		{
			data_transfer.printAllLoData(std::cout);
		}

		//update counter
		i_holo_solves++;

	}




	return 0;
}

