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


#if !defined(_RESIDUALSOURCE_H)
#define _RESIDUALSOURCE_H

#include "Source.h"
#include <vector>
#include <map>

class ResidualSource : public Source
{
protected:

	double _face_src_total; //magnitude of the face source total, vol src total is stored in base class

	//The boundary conditions are only needed for computing the residual, the result is stored in the face source sampler
	std::vector<std::vector<double>> _bc_dof; //dof for boundary conditions for elements, access using map, this could be done dynamically and deleted each time
	std::map<int, int> _bc_element_to_dof_map;  //key: element, value: index in bc_dof array

	//Elementwise information
	std::vector<std::vector<double>> _residual_element_LD_values; //each member contains average, slope x, slope mu;
	std::vector<std::vector<double>> _residual_face_LD_values; //each memeber contains average, slope x (always zero), and slope in mu for sampling
	std::vector<double> _res_element_mags;  //magnitue of residual source over elements, these are only stored temporarily for derived classes to use, must carefully free the memory
	std::vector<double> _res_face_mags; //magnitude of the residual for each face, only stored temporarily for derived classes to use, must carefully free the memory

	//Functions for building sampler
	void computeElementResidual(ECMCElement1D* element, std::vector<double> & residual_LD_values_el, double & residual_element_magnitude); 
	void computeFaceResidual(ECMCElement1D* element, ECMCElement1D* down_str_elem, std::vector<double> & res_LD_values_face, 
		double & residual_face_magnitude, bool on_boundary = false); //special case for boundary cells, the up_stream_element is NULL
	void computeBCAngularFluxDof(); //determine LD values for equivalent incident flux of boundary conditions and corresponding elements
	double evalDeltaResidualIntegral(std::vector<double> res_LD_face_values, double ang_width, double ang_coordinate); //compute magnitude of delta integrals

	//Functions for sampling position and direction
	void sampleElementSource();
	void sampleFaceSource();

public:

	~ResidualSource();
	ResidualSource(Particle1D* particle);
	virtual void sampleSourceParticle() = 0;
	virtual double getTotalSourceStrength();

};

#endif  //_RESIDUALSOURCE_H
