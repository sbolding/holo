// 
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : HoMesh.h
//  @ Date : 2/10/2014
//  @ Author : srb
//
//
//  IF YOU WANT TO USE COORDINATES: This mesh always assumes that the x coordinate starts at 0.0 and then goes out from there.  
//  This is fine for 1D, but in multi-D you would need to have a lo_mesh.getNodes() function 
//  to determine the center of each cell, etc.
//


#if !defined(_HOMESH_H)
#define _HOMESH_H

#include "ECMCElement.h"
#include "ECMCElement1D.h"
#include "Mesh.h"
#include "DirichletBC1D.h"
#include <ostream>


class HoMesh
{
protected:
	std::vector<ECMCElement1D*> _elements; //I think this can be generalized to ECMCElements by using a virtual destructor and dynamic_cast in particle1D class
	int _n_elems;
	Mesh* _lo_mesh;
	HoMesh(); //don't use
public:
	HoMesh(Mesh* mesh, int n_ang_cells_half_range); //how much to divide the angular cells into
	ECMCElement1D* getElement(int element_id) const; 
	int getNumElems() const;
	std::vector<ECMCElement1D* >* getElements(void);

	//computing angular flux
	void computeAngularFluxes(int n_histories, double total_src_strength = 1.0); //compute moments of angular flux based on n_histories and total source strength (assumed 1 if normalizing to per source particle)
	void printAngularFluxes(std::ostream &out);
	std::vector<int> findUpwindBoundaryCells() const; //return a list of cells who have boundary on their upwind face
	std::vector<DirichletBC1D*> getDirichletBCs() const;
};

#endif  //_HOMESH_H
