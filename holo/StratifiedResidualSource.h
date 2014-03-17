//  There is a special case for when the downstream cell is less refined than the upstream cell that is handled in
// both the samping routine and in the constructor.  It effectively splits the less refined cell into two ghost cells and then samples
// the angle based on the LD values of the two ghost cells along the upwind face.  The two ghost cells LD values are stored in
// the vector from bottom to top.
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : ResidualSource.h
//  @ Date : 2/8/2014
//  @ Author : 
//
//

#if !defined(_STRATIFIEDRESIDUALSOURCE_H)
#define _STRATIFIEDRESIDUALSOURCE_H



#include "Source.h"
#include "ResidualSource.h"
#include <iostream>
#include <vector>

class StratifiedResidualSource : public ResidualSource
{
protected:

	unsigned int _n_samples_per_element;//how many samples to take from each element
	double _stratified_probability; //the pdf for the current element of the stratified distribution, for now it is assumed uniform number created in each cell

	//info for current element you are sampling
	int _current_element_id;  //which element you are currently sampling from
	double _curr_el_face_probability; //total probability of a particel being born on face of element, not relative probablity
	double _curr_el_element_probability; //total probability of a particel being born on volumetric portion of residual, not relative probablity
	unsigned int _n_sampled_from_current_element;  //how many samples taken from the current element

public:

	~StratifiedResidualSource();
	StratifiedResidualSource(Particle1D* particle, int & n_histories); //need to know n_histories to determine how many to put in each cell, will increase number histories to make even match
	virtual void sampleSourceParticle();
};




#endif // _STRATIFIEDRESIDUALSOURCE_H