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
#include "HoMesh.h"
#include "RNG.h"
#include "GlobalConstants.h"
#include "AliasSampler.h"

//forward declarations for friend classes
class Source; 
class LinDiscSource;
class ResidualSource;
class StandardResidualSource;
class StratifiedResidualSource;

class Particle1D
{
protected:

	HoMesh* _mesh;
	unsigned int _method; //type of MC
	RNG* _rng;

	//Sampling properties, use friend classes that handle the source sampling
	friend class Source;
	friend class LinDiscSource;
	friend class ResidualSource;
	friend class StandardResidualSource;
	friend class StratifiedResidualSource;
	Source* _source;
	string _sampling_method; 
	unsigned int _sampling_method_index;

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
	double _element_width_mfp; //in MFP to ray trace easily
	double _element_angular_width; //for tallying

	//data for tracking across elements
	int _current_element_ID;	//which element are you in
	int _current_mat_ID;  //Material of current element
	ECMCElement1D* _current_element; //Pointer to the current ECMC element
	Element* _spatial_element; //The spatial element you are currently in
	size_t _n_elements;		//for sampling which element a particle is born in
	bool _is_dead;		 //for terminating particle history

	//protected methods
	//---------------------------------------------
	//Streaming and collision methods
	double samplePathLength();      //Sample a path length in cm
	double samplePathLengthMFP();   //Sample a path length in units of number of MFP, useful for streaming through many cells
	double sampleAngleIsotropic();	//returns a cosine sampled from uniform distribution
	void sampleCollision();  //Determine if a scatter or an absorption, and then do teh appropriate behavior after that, depending on solver method
	void updateElementProperties();
	void leaveElement();	//Called when leaving an element and moving into next geometrical region
	void terminateHistory(); //kill particle, do other appropriate things

	//tallies
	void scoreFaceTally(); //this doesnt make sense in ECMC context, just for verifying particles are being tracked properly
	void scoreElementTally(double path_start_mfp, double path_end_mfp); //where the track begin and ended, in terms of x-coordinate

	//Sampling the source methods
	void sampleSourceParticle();	//base method, called to create a source particle
	void initializeSamplingSource(); //initialize method that determines where to put the particle

	inline void initializeHistory(); //This is in the particle class to ensure it is called everytime
	
public:

	//constructors
	Particle1D(HoMesh* mesh, RNG* rng, string method_str, string sampling_method); //Standard constructor, pass a pointer for rng to make sure 
	//you dont resample random numbers, method_str is which method
	//to use from HoSolver, sampling method for stratified or not
	
	//public functions
	double getRandNum();
	void runHistory();
	void printParticleBalance(int n_histories);
	void computeResidualSource();
	double getTotalSourceStrength(); //get the total source strength from the source


};

#endif  //_PARTICLE1D_H
