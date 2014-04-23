/*
	Class contains all of the information needed to handle a dirichlet BC.

	The boundary condition knows which element and node it belongs to.  The element
	has to figure out how to handle the boundary condition when it constructs the 
	local matrices.  This may be done seperately from the orginial matrix construction
	and just overwrite the values it has set.

*/

#ifndef _DIRICHLETBC1D_H
#define _DIRICHLETBC1D_H

#include "Element.h"
#include "Node.h"
#include "FixedSourceFunctor.h"

class DirichletBC1D
{
private:

	//TODO you dont really need to know the node for 1D, but it will help in 2D
	Node *  _node;			//The node for which the BC is specified on
	Element * _element;		//The element ''						  "
	double _value_current;	//Value of boundary condition based on incident current
	double _inc_flux_avg, _inc_flux_mu; //angular flux avg and moment, constructors must compute these
	int _id;			    //BC id

	//Never used constructor and copiers
	DirichletBC1D();
	DirichletBC1D operator=(const DirichletBC1D &);
	DirichletBC1D(const DirichletBC1D &);

public:

	DirichletBC1D(int id, Element* element, Node* node, double incid_current); //Construction based on specified incoming current
	DirichletBC1D(int id, Element* element, Node* node, const FixedSourceFunctor & q); //construct based on evaluation of q at boundaries.  Useful for MMS solutions

	//Public access functions
	int getID() const;
	double getCurrent() const;
	int getElementID() const;
	Element* getElement() const;
	Node* getNode() const;
	int getNodeID() const;

};

#endif