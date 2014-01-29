//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : Particle1D.cpp
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
//


#include "Particle1D.h"
#include <iostream>

Particle1D::Particle1D(Mesh* mesh, RNG* rng)
{
	_rng = rng;
	_mesh = mesh;
	_position = -999999999999.; //initialize outside teh domain to check that source samplign works correctly
}

//Sample a path length in cm
double Particle1D::samplePathLength()
{
	return -1.*log(_rng->rand_num())*_mfp_tot;
}

//Sample a path length in units of number of MFP, useful for streaming through many cells
double Particle1D::samplePathLengthMFP()
{
	return -1.*log(_rng->rand_num());
}

//Update the particle position based on a path length that has been sampled, make sure hasnt streamed out of cell
void Particle1D::updatePosition(double path_length)
{
	std::cerr << "Need to check and make sure mu has not already been resampled" << std::endl;
	exit(1);
	
	

	_position += path_length*_mu; 
}

void Particle1D::streamThroughCell()
{

}

//Determine if a scatter or an absorption, and then do teh appropriate behavior after that, depending on the mode
void Particle1D::sampleCollision()
{
	switch (_method)
	{
	case HoMethods::HOLO_ECMC:
		//then a pure absorber problem
		//tally the scores, and end the history
	case HoMethods::HOLO_STANDARD_MC:
		//usual MC sample scattering event, angle, and new direction
	case HoMethods::STANDARD_MC:
		//usual MC sample scattering even, angle, and new direction
		_mu += _rng->rand_num() * 2;
	default:
		std::cerr << "Input an incorrect mode of operation\n";
		exit(1);
	}

}

void Particle1D::runHistory()
{
	sampleSourceParticle();

	//Stream teh particle across a cell
	//Determine interactions etcs
	//Once particle leaves the cell, tally the events appropriately
}

//return random numbers for use by source or whoever needs one
double Particle1D::getRandNum() const
{
	return _rng->rand_num();
}

