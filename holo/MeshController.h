//
// Batch errors is reset each time the mesh is adapted, thus the length of it can be checked quickly to determine if
// refinement is necessary
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : MeshController.h
//  @ Date : 2/23/2014
//  @ Author : SRB
//
//


#if !defined(_MESHCONTROLLER_H)
#define _MESHCONTROLLER_H

#include "HoMesh.h"
#include <vector>
#include <map>
#include "ECMCElement1D.h"
#include "Controller.h"
#include "GlobalConstants.h"
#include <cmath>
#include <algorithm>

struct ElementNeighbors
{
	std::vector<ECMCElement1D*> _element_pntrs;// (4, NULL);
	ECMCElement1D* & _left; //element to the left of the current element (in negative x)
	ECMCElement1D* & _minus; //elemetn below the current element (in negative mu) 
	ECMCElement1D* & _right; //element to the right of current element (in positive x)
	ECMCElement1D* & _plus; //element above the current element

	//can make more general by having the constructor take an argument for 
	//the size, and based on size only initialize appropriate elements, set
	//rest to nulls.  Also, could use element id's instead of pointers, would
	//need to change the way that findNeighbors is implemented to handle ID's
	//and case of NULL ds_elements

	ElementNeighbors() :
		_element_pntrs(4,NULL), //currently 4, so only really works for 1D
		_left(_element_pntrs[0]),
		_minus(_element_pntrs[1]),
		_right(_element_pntrs[2]),
		_plus(_element_pntrs[3])
	{}

	ElementNeighbors& operator = (const ElementNeighbors& copy)
	{
		_element_pntrs = copy._element_pntrs;
		return *this;
	}
};

class MeshController
{
protected:

	std::vector<double> _batch_residual_norms; //errors calculated for previous batches, number to keep set by user
	double _required_conv_rate;	//Required exponential convergence rate alpha, i.e., e^-alpha*batch_number
	int _n_batches_to_check; //how many batches to keep and average for checking convergence, default of 3
	std::map<int, ElementNeighbors> _connectivity_array; //for each cell that does NOT have children, which cells are on its face (for 1D 4 cells), indexed by ID so that this can be used for LoElements more easily
	std::vector<int> _newly_refined_elements; //a list of elements that were refined this iteration, needed for updating connectivity array
	HoMesh* _mesh; //ho Mesh
	MeshController(); //don't use

	//protected functions
	void createConnectivityArray(); //go through the mesh and update the connectivity array
	double computeElementJumpError(int elem_id); //compute jump errors for a single, active, element
	double computeJumpIntegral(std::vector<double> avg_slope_1, std::vector<double> avg_slope_2, double width); //DOF are avg, then slope, along face
	void refineElement(int elem_id);
	void updateConnectivityArray(int refined_element_id); //will find new neighbors list for refined element, as well as all elements it points to 
	ElementNeighbors findNeighbors(int elem_id); //find the neighbors of an element, sets to NULL if they don't exist
	
public:

	MeshController(HoMesh* mesh, double exp_convergence_rate, int n_batches_to_check); //n_batches_to_check is how many batches to check convergence on 
	void refineMesh();
	bool meshNeedsRefinement(); //Check is refinement needed?
	void storeResidualNorm(double L1_norm_of_residual); //add the L1_norm of the residual to the list
};

#endif  //_MESHCONTROLLER_H
