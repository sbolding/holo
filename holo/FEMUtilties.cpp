#include "FEMUtilities.h"

std::vector<double> GaussQuadrature::_weights = { 1., 1. };
std::vector<double> GaussQuadrature::_gauss_points = { -0.5773502691896257, 0.5773502691896257 };

GaussQuadrature::GaussQuadrature()
{
	_n_points = 2;
}

GaussQuadrature::GaussQuadrature(int n_gauss_points)
{
	if (n_gauss_points != 2)
	{
		std::cerr << "Only 2 point gauss quadrature available\n";
		exit(1);
	}
	_n_points = 2;
}

int GaussQuadrature::getNumPoints()
{
	return _n_points;
}

std::vector<double> GaussQuadrature::getQuadraturePoints(double center, double width)
{
	std::vector<double> points(_n_points);
	if (width < 0.0)
	{
		std::cerr << "Passed in a negative width to quadrature, in FEMutitlites.cpp\n";
		exit(1);
	}
	double a = center - 0.5*width;
	double b = center + 0.5*width;
		
	for (int i = 0; i < _n_points; i++)
	{
		points[i] = a + 0.5*width + _gauss_points[i] * 0.5*width;
	}
	return points;
}

std::vector<double> GaussQuadrature::getQuadratureWeights(double center, double width)
{
	double a = center - 0.5*width;
	double b = center + 0.5*width;

	std::vector<double> wgts(_n_points);
	for (int i = 0; i < _n_points; ++i)
	{
		wgts[i] = _weights[i] * 0.5*width;
	}
	return wgts;
}

