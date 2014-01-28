#include "Mesh.h"
#include "MaterialConstant.h"
#include <iostream>
#include "LoSolver1D.h"
#include "HoSolver.h"
#include <fstream>

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
	double width = 5;
	double ext_source = 2;
	int num_elems = 35;
					  // ID, sig_a, sig_s
	MaterialConstant mat(10, 0.5, 0.6);

	//Create the mesh and elements;
	Mesh mesh_1D(dimension, num_elems, width, &mat);
	mesh_1D.setExternalSource(ext_source);
	mesh_1D.print(cout);

	//Temporarily hard coded monte carlo parameters
	int n_histories = 10;

	//Solve the low order system
	ho_solver = new HoSolver(&mesh_1D, n_histories, ext_source);
	ho_solver->solveSystem();

	//temporary return
	system("pause");

	return 0;

	//Assemble and solve the system
	lo_solver = new LoSolver1D(&mesh_1D);
	lo_solver->solveSystem();

	//Update lo order system to current solution
	lo_solver->updateSystem();

	//Print scalar flux
	mesh_1D.printLDScalarFluxValues(cout);
	mesh_1D.printLDScalarFluxValues(lo_out_file);

	//Produce scalar flux on nodes.
	system("pause");

	return 0;
}