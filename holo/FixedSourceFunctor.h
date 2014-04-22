// Class to handle arbitrary external SourceFunctors by overriding 
// the evaluate function.  Uses quadrature to compute 
// the moments. This is a virtual class, that needs to have
// the "getValue" function overridden by derived classes.  The
// getValue function is the external SourceFunctor evaluated at 
// some set of cooridnates, passed in as a vector.  To
// use a function read from a file, such as using boost, then
// just need to override the constructor for the spatial moments
// class.
//
// This class can create Dirichlet boundary condition classes
// based on function evaluations, if desired
//
//
//  @ Project : Untitled
//  @ File Name : FixedSourceFunctorFunctor.h
//  @ Date : 4/22/2014
//  @ Author : srb
//
//

#ifndef _FIXEDSOURCEFUNCTOR_H
#define _FIXEDSOURCEFUNCTOR_H

#include <vector>
#include "FEMUtilities.h"
#include "DirichletBC1D.h"

class FixedSourceFunctor
{
protected:
	//copy constructor
	FixedSourceFunctor(const FixedSourceFunctor&);
	FixedSourceFunctor& operator=(const FixedSourceFunctor&);

public:

	//Constructors
	FixedSourceFunctor();
	virtual ~FixedSourceFunctor();

	//Primary Functions
	std::vector<double> getLoNodalValues(const std::vector<double> & coors, const std::vector<double> & dimens) const; //only need spatial dimensions, must correspond to evals dimensions
	std::vector<double> getHoMoments(const std::vector<double> & coors, const std::vector<double> & dimens) const;

	//virtual function evalation
	virtual double getValue(const std::vector<double> & coordinates) const = 0; //evaluate function at coordinates, e.g., x,mu, in domain
};






#endif //_FIXEDSourceFunctor_H