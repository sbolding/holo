//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : ECMCElement1D.h
//  @ Date : 2/10/2014
//  @ Author : 
//
//


#if !defined(_ECMCELEMENT1D_H)
#define _ECMCELEMENT1D_H

#include "ECMCElement.h"

class ECMCElement1D : public ECMCElement
{
protected:

	//references for easy use in 1D calculations
	double & _psi_average;	    //Average flux
	double & _psi_x;			//Moment in x direction
	double & _psi_mu;			//Moment in y direction
	double & _width_spatial;	//delta_x
	double & _width_angle;		//delta_mu
	double & _position_center;	//x_i (center of element in x)
	double & _mu_center;		//mu_i (center of element in mu)
	ECMCElement1D();			//don't call me

	//Pointer to the downstream element (different for positive and negative flow)
	ECMCElement1D* _down_stream_element;

public:
	ECMCElement1D(Element* element, ECMCElement1D* down_stream_element, std::vector<double> dimensions, std::vector<double> coordinates);
	virtual ECMCElement1D* getDownStreamElement() const;
};

#endif  //_ECMCELEMENT1D_H
