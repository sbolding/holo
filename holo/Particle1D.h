//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  High Order Solver
//  @ File Name : Particle1D.h
//  @ Date : 1/27/2014
//  @ Author : Simon R Bolding
//
//


#if !defined(_PARTICLE1D_H)
#define _PARTICLE1D_H

#include <string>
#include "Mesh.h"
#include "RNG.h"
#include "FaceTally.h"
#include "ElementTally.h"
#include "GlobalConstants.h"

class Particle1D
{
protected:

	Mesh* _mesh;

	//general particle properties
	double _position_mfp;
	double _mu;
	double _weight;

	//debugging properties, generally won't be used
	int _n_scat;
	int _n_abs;
	int _n_leak;
	int _n_terminations; //how many times was terminateHistory() called?

	//material properties for the current element
	double _sigma_tot;
	double _mfp_tot;
	double _scat_ratio;
	double _sigma_abs;
	double _sigma_scat;
	unsigned int _method;
	double _element_width_mfp; //to ray trace easily

	//data for tracking across elements
	int _current_element;	//which element are you in
	size_t _n_elements;		//for sampling which element a particle is born in
	int _element_mat_ID;  //the material ID of the current element
	bool _is_dead;		 //for terminating particle history
	RNG* _rng;			

	//tally arrays
	std::vector<FaceTally>* _face_tallies;  //vector of all the face tallies, indexed using connectivity array
	std::vector<ElementTally>* _element_tallies; //vector of all the volume tallies, indexed using connectivity array

	//protected methods
	//---------------------------------------------
	//Streaming and collision methods
	double samplePathLength();
	double samplePathLengthMFP();
	double sampleAngleIsotropic();	//returns a cosine sampled from uniform distribution
	void sampleCollision();
	void streamAcrossGeometry();  //Most tallies get called in this routine
	void updateElementProperties();
	void leaveElement();	//Called when leaving an element and moving into next geometrical region
	void terminateHistory(); //kill particle, do other appropriate things

	//tallies
	void scoreFaceTally();
	void scoreElementTally();
	void scoreTallies();
	//Sampling the source methods
	void sampleSourceParticle();
	void sampleLinDiscontSource(std::vector<double>);
	void sampleConstExtSource();
	void sampleDirichletBCSource();
	
public:

	//constructors
	Particle1D(Mesh* mesh, RNG* rng, string method_str); //Standard constructor, pass a pointer for rng to make sure you dont resample random numbers, method_str is which method to use from HoSolver
	
	//public functions
	double getRandNum();
	void runHistory();
	void printParticleBalance(int n_histories);


};

#endif  //_PARTICLE1D_H
