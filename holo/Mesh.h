// Contains the code for a mesh.  Currently only a 1D mesh based on the number of elems
// and the size of the domain is allowed. It also creates the Element objects based on the given mesh
//
// NOTE: the read access functions for this class return pointers, so data can be changed externally,
// however this is much faster than passing the element list around.
//
// TODO: need a function that returns individual elements by index.  Probably have a function that returns an
// iterator you can use to access the elements individually.  This will encapsulate data much more effectively
//
//  @ Project : Untitled
//  @ File Name : Mesh.h
//  @ Date : 11/1/2013
//  @ Author : SRB
//
//


#if !defined(_MESH_H)
#define _MESH_H

#include <fstream>
#include <cstdlib>
#include <iostream>
#include "Element.h"
#include "Node.h"
#include <vector>
#include "MaterialConstant.h"
#include "DirichletBC1D.h"

class Mesh
{
protected:

	//Variables
    int _n_nodes;			//Total number of nodes in a mesh
    int _n_elems;			//Number of elements in a mesh
    int _n_edges;			//number of edges, for 1D this is just number of nodes
	int _n_dirichlet_bc;	//Number of essential BC
    //std::ifstream & _input_file; //reference to input file TODO, currently no file needed
    std::vector<Element *> _elements; //Pointers to element objects
	std::vector<Node *> _nodes; //Pointers to node objects
	std::vector<DirichletBC1D *> _dirichlet_bcs; //Pointers to bc's;
    int _dim;					      //dimension of the problem
    std::vector<std::vector<int>> _connectivity_array; //first index is element, second is face of element
	std::vector<MaterialConstant *> _mat_list;

public:

    Mesh();                                          //Default constructor
    Mesh(int dim, int number_elements, double width, MaterialConstant* material);   // pass in the number of elements as well as dimension
    Mesh(int dim, std::ifstream  input_file); //Constructor with a input file TODO not yet implemented

	//Printing munctions
	void print(std::ostream &out) const;      //echo mesh parameters
	void printLDScalarFluxValues(std::ostream &out) const; //print results of lo order flux LD values (phi_r and phi_l)

	//Reading functions
	int getNumElems() const;
	int getNumNodes() const;
	int getNumEdges() const;
	int getNumDirichletBC() const;
	std::vector<Element* >* getElements(void);		//WARNING returns a pointer to elements, you can change it outside of class
	std::vector < DirichletBC1D *> getDirichletBCs(void);
	Element* getElement(int element_id) const;
	int getFaceIndex(int element_id, int face_id) const; //get a face index using the connectivity array for tallying, face_id is 0 for left 1 for right, etc.

	//Setting functions
	void setExternalSource(double) const;
	//void setExternalSource(void(*func_handle)); not implemented yet
};

#endif  //_MESH_H
