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
#include <iostream>
#include "Particle1D.h"
#include "Source.h"
#include "RNG.h"
#include "Controller.h"
#include "AverageCosineData.h"
#include "MeshController.h"
#include "StratifiedResidualSource.h"
#include "StandardResidualSource.h"
#include "LinDiscSource.h"
#include "CWDParticle1D.h"
#include "FixedSourceFunctor.h"
#include <vector>

class HoSolver
{
protected:

	Mesh* _lo_mesh;	//pointer to the mesh to be used
	HoMesh* _ho_mesh; //pointer to the ho mesh to be created;
	MeshController* _mesh_controller; //For adapting the high order mesh, dynamic because it is based on ho_mesh object
	Source* _source; //pointer to the source to be used for the problem, the source will get recomputed often

	//other problem parameters
	int _n_histories;	//number of histories
	int _n_batches; //number of batches for ECMC simulations
	std::string _solver_mode_str;  ////either "holo-ecmc", 'holo-standard-mc', or 'standard-mc'
	int _solver_mode_int; //integer solver mode, use HoSolver map from GlobalConstant.h to map the string to int
	Particle1D* _particle;	//One particle that has all the methods to stream, cross interfaces, etc.
	RNG _rng; //random number generator

	//source parameters
	string _sampling_method;
	unsigned int _sampling_method_index;  //corresponding value in HoMethods related to the sampling method string

	//residual source handling functions
	void computeResidualSource();
	void initializeSamplingSource(FixedSourceFunctor & q);

	//For computing projected half range fluxes on coarsest mesh, useful for data transfer to low order date
	//Currently they are projected onto a BILINEAR finite element to ensure the LO data is computed exactly
	//All other calculations are identical
	std::vector<std::vector<double>> _psi_plus_dof; //row: spatial element, column: angular flux dof
	std::vector<std::vector<double>> _psi_minus_dof; //row: spatial element, column: angular flux dof
	
	//For data transfer and printing angular flux
	void computeProjectedAngularFlux(); //computes the angular flux dof  element in half range
	std::vector<double> getScalarFluxDOF(int spatial_elem_id) const; //get the scalar flux LD DOF for a particular spatial element
		
public:

	HoSolver(); //Default constructor, should probably never be called
	~HoSolver();
	HoSolver(Mesh* _mesh, int n_histories, int n_bins_half_range, string solver_mode, string sampling_method,
		double required_exp_convg_constant, int n_batches = 1, int n_batches_to_avg = 3); //How many angular cells to split each spatial mesh cell into initially, n_batches_to_keep is for convergence check
	void solveSystem(std::ostream & out = std::cout); //run the high order problem, default output to screen
	void updateSystem(); //compute the angular fluxes, currently unused

	//reader, printer, and interface functions
	void printAllTallies(std::ostream &out) const;
	void printProjectedScalarFlux(std::ostream &out) const;

	//data transfer functions
	virtual void getLoData1D(LoData1D & lo_data, int spatial_element_id);		//calculate the LoData parameters based on projected half range fluxes;
};

#endif  //_HOSOLVER_H
