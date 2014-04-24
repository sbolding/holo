// Contains class for 1D elements.  Consist of only 2 nodes.  Material is assumed to be uniform
// across the element for now (cross sections).  
//
// The 1D element has a container for the 
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : Element1D.h
//  @ Date : 11/1/2013
//  @ Author : 
//
//


#if !defined(_ELEMENT1D_H)
#define _ELEMENT1D_H

#include "Dof.h"
#include "Element.h"
#include "MaterialConstant.h"
#include <iostream>
#include "numMatrixFull.h"
#include "numVector.h"
#include "GlobalConstants.h"

using EquationMaps1D::DOF_MAP;  //

class Element1D : public Element
{
private:

	double _h;	//Width of element
	Element1D(); //Default constructor, never should be called
	void addElementUpwindedTerms(numMatrix &M);

	// --- Alias Accessors ----
	//These are here just to make ordering equations less confusing, and to alias DOF vector.  
	//In 1D when constructing matrices access these, otherwords, access the _elem_dof vector.
	Dof & _phi_left_plus; //for setElementDof, this is DOF_MAP[0]th member
	Dof & _phi_right_plus;  //for setElementDof, this is DOF_MAP[1]th member
	Dof & _phi_left_minus; //etc.
	Dof & _phi_right_minus; //

	//Methods
	virtual void getScalarFluxValues(double alpha, std::vector<double> &scalar_flux_values,
		std::vector<double> &locations) const;  //Return the flux values on the faces (nodes in 1D)

public:

	Element1D(int id, MaterialConstant* mat, std::vector<Node *> nodes);
	virtual void print(std::ostream &out) const;
	virtual void setElementDof(std::vector<double> elem_dofs);  //Set element DOF's based on solution vector TODO not verified yet
	virtual std::vector<int> getEqnNumbers(void) const;		//Get the equation numbers corresponding to each DOF
	virtual void getElementMomentMatrix(numMatrix* M, numVector* b, 
		std::vector<int> &eqns, FixedSourceFunctor* q=NULL) const;	//Returns the moment equations for an element, including upwinded terms.  
	virtual void addDirichletBC(numVector* b, std::vector<int> &eqns,
		double value, Node* node) const;  //Will return values to modify global matrix and vector correctly
	virtual void getPosUpwinding(std::vector<double> &pos_upwind_values, int &eqn, std::vector<int> &cols) const;		//Will return upwinding terms in + dir, which row, and which cols
	virtual void getNegUpwinding(std::vector<double> &neg_upwind_values, int &eqn, std::vector<int> &cols) const;		//Will return upwinding terms in - dir, which row, and which cols
	virtual void getScalarFluxLinDisc(std::vector<double> &scalar_flux_values,
		std::vector<double> &locations) const; //Returns scalar flux on faces using Linear Discontinuous approximation
	virtual void getScalarFluxHOClosure(std::vector<double> &scalar_flux_values,
		std::vector<double> &locations) const; //This is for scalar flux based on alpha closure, shouldnt be used except verification
	virtual std::vector<double> getScalarFluxNodalValues(void) const; //returns a vector of the scalar flux nodal values needed for HO solver
	virtual void printLDScalarFluxValues(std::ostream &out) const; //Print out the LD scalar flux and nodes
	virtual std::vector<double> getElementDimensions() const; //For 1D, return the width of the element
	virtual std::vector<double> getNodalCoordinates(void) const; //returns a vector of the nodal coordinates of the element
	virtual std::vector<double> getSpatialCoordinates(void) const; //returns a vector of the coordinates corresponding to the center of spatial element
};

#endif  //_ELEMENT1D_H
