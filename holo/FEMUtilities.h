//	File contains quadrature and other useful functions for computations on
//  2D LD elements.
//
//  @ File Name : FEMUtilies.h
//  @ Date : 3/18/2013
//  @ Author : Simon R Bolding
//
//

#ifndef _FEMUTILITIES_H
#define _FEMUTILITIES_H

#include <vector>
#include <iostream>

class GaussQuadrature
{
private:

	//Static variables, values set in constructor, currently only 2D
	static std::vector<double> _weights; //weights corresponding to each gauss point
	static std::vector<double> _gauss_points; //1D gauss points on reference element [-1,1]
	static int _n_points;

public:

	GaussQuadrature(); //default 2 gauss points
	GaussQuadrature(int n_gauss_points); //currently only works for 2
	std::vector<double> getQuadraturePoints(double a, double b); //returns gauss points over the interval a to b
	std::vector<double> getQuadratureWeights(double a, double b); //note, for a higherDimension domain weights have to be multiplied together
	int getNumPoints();

};


namespace FEMUtilties
{
}








#endif // _FEMUTILITIES_H