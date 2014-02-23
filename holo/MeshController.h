//
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

class MeshController
{
protected:

	std::vector<double> _batch_errors; //errors calculated for previous batches, number to keep set by user
	double _required_conv_rate;	//Required exponential convergence rate alpha, i.e., e^-alpha*batch_number
	std::map<int, std::vector<int>> _connectivity_array; //for each cell that does NOT have children, which cells are on its face (for 1D 4 cells), indexed by ID so that this can be used for LoElements more easily
	HoMesh* _mesh; //ho Mesh
	void computeJumpError();

public:

	MeshController(HoMesh* mesh);
	void refineMesh();
};

#endif  //_MESHCONTROLLER_H
