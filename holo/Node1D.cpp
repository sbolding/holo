//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : 1DNode.cpp
//  @ Date : 11/1/2013
//  @ Author : 
//
//


#include "Node1D.h"

//Intialize based on x pnt and id
Node1D::Node1D(int id, double pnt) : Node(id)
{
	_x = pnt;
}

void Node1D::print(std::ostream &out)
{
	out << " ID = ";
	out.width(4);
	out << _id;
	out << ", X = ";
	out.width(4);
	out.left;
	out << _x;
	out << std::endl;
}
