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
	double & _psi_average;	    //Average flux, ALL FLUXES CURRENTLY HAVE UNITS OF PER SOURCE PARTICLE
	double & _psi_x;			//Moment in x direction
	double & _psi_mu;			//Moment in mu direction
	double & _width_spatial;	//delta_x
	double & _width_angle;		//delta_mu
	double & _position_center;	//x_i (center of element in x)
	double & _mu_center;		//mu_i (center of element in mu)
	ECMCElement1D();			//don't call me

	//Pointer to the downstream element (different for positive and negative flow)
	ECMCElement1D* _down_stream_element;

public:
	ECMCElement1D(int id, Element* element, ECMCElement1D* down_stream_element, std::vector<double> dimensions, std::vector<double> coordinates);

	//methods only used by Particle1D
	virtual ECMCElement1D* getDownStreamElement() const;
	double getAngularWidth() const;	//width of element in angle
	double getSpatialWidth() const; //width of element in x
	double getAngularCoordinate() const; //center of element in angle
	double getSpatialCoordinate() const; //center of element in x
	void incrementTallyScores(double weight, double path_length_cm, double normalized_dir_cosine,
		double volume, double normalized_position);

	//virtual methods
	virtual void printAngularFluxDOF(std::ostream & out) const;
	virtual void computeAngularFLuxDOF(int n_histories, double ext_src_str = 1.0);

};
#endif  //_ECMCELEMENT1D_H
