#include "DirichletBC1D.h"

//Constructor
DirichletBC1D::DirichletBC1D(int id, Element* element, Node* node, double val)
{
	_id = id;
	_element = element;
	_node = node;
	_value = val;
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

double DirichletBC1D::getValue() const 
{
	return _value;
}

int DirichletBC1D::getID() const
{
	return _id;
}

int DirichletBC1D::getNodeID() const
{
	return _node->getID();
}
