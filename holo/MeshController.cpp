//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : MeshController.cpp
//  @ Date : 2/23/2014
//  @ Author : SRB
//
//


#include "MeshController.h"
#include <cmath>

MeshController::MeshController(HoMesh* mesh, double exp_conv_constant, int n_batches_to_check):
	_required_conv_rate(exp_conv_constant),
	_mesh(mesh), 
	_n_batches_to_check(n_batches_to_check),
	_batch_residual_norms()
{
	createConnectivityArray();
}

void MeshController::computeJumpError()
{

}

void MeshController::createConnectivityArray()
{
	//loop over elements
	std::vector<ECMCElement1D*>* elements = _mesh->getElements();
	std::vector<ECMCElement1D*>::iterator it_el;
	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		//TODO create the connectivity array
	}
}

void MeshController::storeResidualNorm(double residual_norm)
{
	if (_batch_residual_norms.size() == (_n_batches_to_check+1)) //one extra since you need ratios to check refinement
	{
		_batch_residual_norms.erase(_batch_residual_norms.begin()); //remove the oldest element
	}
	_batch_residual_norms.push_back(residual_norm);
}

void MeshController::refineMesh()
{
	//Refine mesh
	std::vector<ECMCElement1D*> new_elements;
	_mesh->getElement(1)->refine(_mesh->_n_elems-1); //pass the id of last element made
	//add the new elements to the list and update number of elements
	new_elements = _mesh->getElement(1)->getChildren();
	_mesh->_elements.insert(_mesh->_elements.end(), new_elements.begin(), new_elements.end());
	_mesh->_n_elems += new_elements.size();
	_batch_residual_norms.erase(_batch_residual_norms.begin(), _batch_residual_norms.end() - 1); //clear all but the last one

}

bool MeshController::meshNeedsRefinement()
{
	if (_batch_residual_norms.size() != (_n_batches_to_check + 1))
	{
		return false; //not enough batches to check convergence, because of noise
	}
	else
	{
		double alpha_avg = 0.0;
		for (int i = 0; i < _n_batches_to_check; ++i) //loop over batches
		{
			alpha_avg += std::log(_batch_residual_norms[i] / _batch_residual_norms[i+1]);
		}
		alpha_avg /= (float)_n_batches_to_check;

		//check convergence
		if (alpha_avg < _required_conv_rate) //if error increases, alpha will be negative and go back up
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}