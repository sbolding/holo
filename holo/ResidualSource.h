//
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
#include "AliasSampler.h"
#include <vector>

class ResidualSource : public Source
{
protected:
	AliasSampler* _element_source;  //Sampler to determine which source you are in
	AliasSampler* _face_source;	    //Sampler to determine which face you are on

	double _face_src_total; //magnitude of the face source total

	std::vector<std::vector<double>> _residual_element_LD_values; //each member contains average, slope x, slope mu;
	std::vector<std::vector<double>> _residual_face_LD_values; //each memeber contains average and slope in mu for sampling

	//Functions for building sampler
	void computeElementResidual(ECMCElement1D* element, std::vector<double> & residual_LD_values_el, double & residual_element_magnitude); 
	void computeFaceResidual(ECMCElement1D* element, std::vector<double> & res_LD_values_face, double & residual_face_magnitude,
		bool on_boundary = false); //special case for boundary cells

	//Functions for sampling position and direction
	void sampleElementSource();
	void sampleFaceSource();
	double evalLinDiscFunc2D(std::vector<double> dof, ECMCElement1D* element, double x, double mu); //evaluate based on dof and dimensions of particular element

public:
	~ResidualSource();
	ResidualSource(Particle1D* particle, string sampling_method);
	virtual void sampleSourceParticle();

};

#endif  //_RESIDUALSOURCE_H
