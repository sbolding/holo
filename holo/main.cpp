#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>
#include "RNG.h"
#include "DataTransfer.h"

int main()
{
	using std::cout;
	using std::endl;

	cout << "We're rolling..." << endl;
	LoSolver* lo_solver;
	HoSolver* ho_solver;

	//open output file
	ofstream lo_out_file("./lo_results.txt", std::ofstream::out);
	if (!lo_out_file)
	{
		std::cerr << "Can't open file" << endl;
		exit(1);
	}

	//Temporarily hard coded dimensions until there is stuff for reading from input file
	int dimension = 1;
	double width = 2.0; //cm
	double ext_source = 2.0; //(p/(sec cm^3)), do not use non-zero values << 1, or some logic may be wrong currently
	int num_elems =  2;
	int n_ang_elements = 2; //number angles in half ranges
	//Temporarily hard coded monte carlo parameters
	int n_histories = 20000; //50000000
	int n_batches = 20;
	double exp_convg_rate = 0.20;
	string solver_mode = "holo-ecmc"; //"standard-mc", "holo-ecmc", "holo-standard-mc"
	string sampling_method = "stratified";
					  // ID, sig_a, sig_s
	MaterialConstant mat(10, 0.4, 0.50);

	//Create the mesh and elements;
	Mesh mesh_1D(dimension, num_elems, width, &mat);
	mesh_1D.setExternalSource(ext_source);
	mesh_1D.print(cout);

	size_t n_holo_solves = 25;

	//Assemble the LoSystem
	lo_solver = new LoSolver1D(&mesh_1D); //uses default some estimated lo order parameters and LD

	size_t i_holo_solves = 0; //counter
	while (true)
	{
		//solve lo order system
		lo_solver->solveSystem();

		//Update lo order system to current solution
		lo_solver->updateSystem();

		//Print LO scalar flux estimate
		mesh_1D.printLDScalarFluxValues(cout);
		mesh_1D.printLDScalarFluxValues(lo_out_file);

		if (i_holo_solves == n_holo_solves) //do one extra LO solve, because it is the only solution really being calculated
		{
			break;
		}

		ho_solver = new HoSolver(&mesh_1D, n_histories, n_ang_elements, solver_mode, sampling_method, exp_convg_rate, n_batches);
		ho_solver->solveSystem();
		ho_solver->updateSystem();

		//Transfer HO data to the LO system
		DataTransfer data_transfer(ho_solver, &mesh_1D);
		data_transfer.updateLoSystem();

		//update counter
		i_holo_solves++;
	}
	
	return 0;
}