//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : HoSolver.h
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
//


#if !defined(_HOSOLVER_H)
#define _HOSOLVER_H

#include "Mesh.h"
#include <vector>
#include "FaceTally.h"
#include "ElementTally.h"
#include "Particle1D.h"
#include "Source.h"
#include "SourceConstant.h"
#include "RNG.h"

class HoSolver
{
protected:

	Mesh* _mesh;	//pointer to the mesh to be use
	std::vector<FaceTally> _face_tallies;  //vector of all the face tallies, indexed using connectivity array
	int _n_histories;	//number of histories
	std::vector<ElementTally> _element_tallies; //vector of all the volume tallies, indexed using connectivity array
	std::string _solver_mode;  //TODO either "holo-scat", "mc", or "holo-ecmc"
	Particle1D* _particle;	//One particle that has all the methods to stream, cross interfaces, etc.
	Source* _external_source; //external source
	Source* _scattering_source; //the in scattering source
	RNG _rng; //random number generator
	
public:

	HoSolver(); //Default constructor, should probably never be called
	HoSolver(Mesh* _mesh, int n_histories, double ext_source); //For a constant external source
	void solveSystem(); //run the high order problem

};

#endif  //_HOSOLVER_H
