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

HoSolver::HoSolver()
{
	std::cerr << "You have called HoSolver Default constructor, should never happen" << std::endl;
	exit(1);
}

HoSolver::HoSolver(Mesh* mesh, int n_histories, double ext_source) :
	_rng()
{
	_mesh = mesh;
	_n_histories = n_histories;
	_face_tallies.resize(mesh->getNumEdges());		//initialize tallies to appropriate size
	_element_tallies.resize(mesh->getNumElems());  
	_particle = new Particle1D(mesh, &_rng);

	

}

void HoSolver::solveSystem()
{
	std::cout << "Solving the HO system..." << std::endl;

	//loop over the number of histories
	for (int hist=1; hist <= _n_histories; hist++) 
	{
		_external_source->sampleSource(_particle, _mesh);
		//sample the source distribution
		//Stream teh particle across a cell
		//Determine interactions etcs
		//Once particle leaves the cell, tally the events appropriately
	}

	

	
}
