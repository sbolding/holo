//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : ECMCElement1D.cpp
//  @ Date : 2/10/2014
//  @ Author : 
//
//


#include "ECMCElement1D.h"

ECMCElement1D::ECMCElement1D(Element* element, ECMCElement1D* down_stream_element, std::vector<double> dimensions, std::vector<double> coordinates):
ECMCElement(element,dimensions,coordinates), _psi_average(_ang_flux_dof[0]), _psi_x(_ang_flux_dof[1]), _psi_mu(_ang_flux_dof[2]),
_width_spatial(_dimensions[0]), _width_angle(_dimensions[1]),
_position_center(_coordinates[0]), _mu_center(_coordinates[1])
{
	_down_stream_element = down_stream_element;
	_tally = new ECMCTally(); //intializes the correct size
}

ECMCElement1D* ECMCElement1D::getDownStreamElement() const
{
	return _down_stream_element;
}