// This object assembles the matrix and solves system by looping over all the elements in a mesh
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : LoSolver.cpp
//  @ Date : 11/1/2013
//  @ Author : 
//
//


#include "LoSolver.h"
#include "numMatrixFull.h"
#include "Controller.h"

using namespace std;

void LoSolver::solveLinearSystem()
{
	//Currently just a really slow gaussian elimination solver
	_system_matrix->solve(_system_vec, _sol_vec); 
}

LoSolver::LoSolver(int dim, Mesh* mesh) :
_mesh(mesh),
_dim(dim)
{
	if (dim != 1)
	{
		cerr << "Only 1D for now" << endl;
		exit(2);
	}
	_sol_vec = NULL;
	_system_matrix = NULL;
	_system_vec = NULL;		
}

LoSolver::LoSolver()
{
	//Default constructor
	_sol_vec = NULL;
	_system_matrix = NULL;
	_system_vec = NULL;
}

void LoSolver::deleteMatrixVector()
{
	delete _system_matrix;
	delete _system_vec;
}

void LoSolver::printSystem(std::ostream & out)
{
	using namespace LoController; //Get all bools and debug output from here
	if (WRITE_MATRIX)
	{
		_system_matrix->printCompact(out);
	}
	if (WRITE_LOAD_VECTOR)
	{
		_system_vec->print(out);
	}
}

void LoSolver::updateSystem(void) const
{
	//Map DOF from solution back on to element DOF values
	//TODO this could be really slow
	std::vector<Element *>* elements;
	std::vector<int> eqns;
	std::vector<double> el_dof_values;
	elements = _mesh->getElements();
	int el_n_dof = elements->operator[](0)->getNumDof(); //TODO, this could be different for each element...maybe?

	std::vector<Element *>::iterator it_el;
	it_el = elements->begin();
	for (; it_el != elements->end(); ++it_el)
	{
		el_dof_values.clear();
		eqns = (*it_el)->getEqnNumbers(); //Get equations for element
		for (int i = 0; i < el_n_dof; ++i) 	//Get global solution DOF values
			el_dof_values.push_back(_sol_vec->getCoeff(eqns[i]));

		(*it_el)->setElementDof(el_dof_values);
	}
}
