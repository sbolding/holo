#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>
#include "RNG.h"
#include "DataTransfer.h"

using namespace std;

int main()
{
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
	double width = 2.0;
	double ext_source = 1.;
	int num_elems = 10;
	string solver_mode = "standard-mc"; //"standard-mc", "holo-ecmc", "holo-standard-mc"
					  // ID, sig_a, sig_s
	MaterialConstant mat(10, 0.5, 0.0);

	//Create the mesh and elements;
	Mesh mesh_1D(dimension, num_elems, width, &mat);
	mesh_1D.setExternalSource(ext_source);
	mesh_1D.print(cout);

	//Assemble and solve the Lo system
	/*lo_solver = new LoSolver1D(&mesh_1D);
	lo_solver->solveSystem();

	//Update lo order system to current solution
	lo_solver->updateSystem();*/

	//Temporarily hard coded monte carlo parameters
	int n_histories = 500000; //50000000

	//Solve the low order system
	ho_solver = new HoSolver(&mesh_1D, n_histories, ext_source, solver_mode);
/*	ho_solver->solveSystem();
	ho_solver->printAllTallies(cout); */

	//Transfer HO data to the LO system
	/*DataTransfer data_transfer(ho_solver, &mesh_1D);
	data_transfer.updateLoSystem();*/
	
	//temporary return
	system("pause");

	return 0;

	//Print scalar flux
	mesh_1D.printLDScalarFluxValues(cout);
	mesh_1D.printLDScalarFluxValues(lo_out_file);

	//Produce scalar flux on nodes.
	system("pause");

	return 0;
}