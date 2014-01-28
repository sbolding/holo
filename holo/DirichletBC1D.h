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

class DirichletBC1D
{
private:

	//TODO you dont really need to know the node for 1D, but it will help in 2D
	Node *  _node;			//The node for which the BC is specified on
	Element * _element;		//The element ''						  "
	double _value;			//Value of boundary condition
	int _id;			    //BC id

	//Never used constructor
	DirichletBC1D();

public:

	DirichletBC1D(int id, Element* element, Node* node, double val); //Constructor

	//Public access functions
	int getID() const;
	double getValue() const;
	int getElementID() const;
	Element* getElement() const;
	Node* getNode() const;
	int getNodeID() const;

};

#endif