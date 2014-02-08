//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : DataTransfer.cpp
//  @ Date : 2/7/2014
//  @ Author : 
//
//


#include "DataTransfer.h"

DataTransfer::DataTransfer(HoSolver* ho_solver, Mesh* mesh)
{
	_ho_solver = ho_solver;
	_mesh = mesh;
}

void DataTransfer::updateLoSystem()
{
	//Loop over elements in the mesh
	Element* current_element;
	LoData1D* lo_data;
	
	for (int el = 0; el < _mesh->getNumElems(); ++el)
	{
		current_element = _mesh->getElement(el);
		current_element->setLoData(_ho_solver->getLoData(el));
	}

}

