// Contains all the data an element would need for constructing the low order system.
// 
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : LoData.h
//  @ Date : 11/1/2013
//  @ Author : SRB
//
// TODO: Really need a different class for 1D and 2D, this is hardcoded for 1D currently


#if !defined(_LODATA1D_H)
#define _LODATA1D_H
// const double TWO_THIRDS = 0.666666666666667;
const double TWO_THIRDS = 0.57735026919;

#include <vector>
#include "AverageCosineData.h"

class LoData1D
{
protected:

    double _alpha;				//For 1D the spatial closure constant psi_(i+1/2) = \alpha*<psi>_R + (1 - \alpha)*<psi>_L etc.
	AveragedCosines _surf_avg_cosines; 
	AveragedCosines  _vol_avg_cosines;

public:

	//Constructors
	LoData1D();	//Default constructor.  Here alpha will be set initially to linear Discontinuous and cosines to that of diffusion theory

	//Read and set functions
    double getSpatialClosureFactor() const;		//Return the value of alpha
    AveragedCosines getVolAveragedCos() const; //This is inefficent, but will prevent accidentally changing data
	AveragedCosines getSurfAveragedCos() const; //Return the Surface Averaged Cosine struct
	void setVolAveragedCos(AveragedCosines);	//Must be set initially because changes with every HO solve
	void setSurfAveragedCos(AveragedCosines);	//Must be set initially because changes with every HO solve
	void setSpatialClosureFactor(double);		//Initially is set to 2, which is equivalent to Linear Discon.
	std::vector<double> getExtSourceNodalValues(); //Returns the LD nodal values of the external source

private:
	void init(); //For use by constructors because they intialize most of the data the same way
};

#endif  //_LODATA1D_H
