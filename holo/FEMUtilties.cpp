#include "FEMUtilities.h"



GaussQuadrature::GaussQuadrature()
{
	GaussQuadrature(2);
}

GaussQuadrature::GaussQuadrature(int n_gauss_points)
{
	if (n_gauss_points != 2)
	{
		std::cerr << "Only 2 point gauss quadrature available\n";
		exit(1);
	}
	_n_points = 2;
	_gauss_points = { -0.5773502691896257, 0.5773502691896257 }; //only works in c++11
	_weights = { 1., 1. };
}

int GaussQuadrature::getNumPoints()
{
	return _n_points;
}

std::vector<double> GaussQuadrature::getQuadraturePoints(double a, double b)
{
	std::vector<double> points(_n_points);
	if (a > b)
	{
		std::cerr << "Passed in gauss points that were not in increasing order, in FEMUtilites.cpp\n";
		exit(1);
	}
	double h = b - a;
	for (int i = 0; i < _n_points; i++)
	{
		points[i] = a + 0.5*h + _gauss_points[i] * 0.5*h;
	}
	return points;
}

std::vector<double> GaussQuadrature::getQuadratureWeights(double a, double b)
{
	std::vector<double> wgts(_n_points);
	for (int i = 0; i < _n_points; ++i)
	{
		wgts[i] = _weights[i] * std::abs(b - a);
	}
	return wgts;
}