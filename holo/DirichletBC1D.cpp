#include "DirichletBC1D.h"

//Constructor
DirichletBC1D::DirichletBC1D(int id, Element* element, Node* node, double val)
{
	_id = id;
	_element = element;
	_node = node;
	_value_current = val;

	//Need to set the other stuff
}

DirichletBC1D::DirichletBC1D(int id, Element* element, Node* node, std::vector<double> inc_flux_moments) :
_id(id),
_element(element),
_node(node)
{
	//Calculate the half range current using q
	
	q.getLoNodalValues()
	

}

int DirichletBC1D::getElementID() const
{
	return _element->getID();
}

Element* DirichletBC1D::getElement() const
{
	return _element;
}

Node* DirichletBC1D::getNode() const
{
	return _node;
}

double DirichletBC1D::getCurrent() const 
{
	return _value_current;
}

int DirichletBC1D::getID() const
{
	return _id;
}

int DirichletBC1D::getNodeID() const
{
	return _node->getID();
}


