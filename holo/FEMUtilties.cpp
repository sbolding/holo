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
		wgts[i] = _weights[i] * 0.5*width; //multiply by 0.5 because gauss quadrature is over intervale [-1,1]
	}
	return wgts;
}

void FEMUtilities::convertMomentsToEdgeValues1D(std::vector<double> moment_dof, std::vector<double> & nodal_values) //convert moment dof to edge values
{
	if (moment_dof.size() != 2)
	{
		std::cerr << "Passed in incorrect length of moment vector to FEMUtilties::convertMomentsToEdgeValues1D\n";
		exit(1);
	}
	nodal_values.resize(2);
	nodal_values[0] = moment_dof[0] - moment_dof[1];
	nodal_values[1] = moment_dof[0] + moment_dof[1];
}

void FEMUtilities::convertAvgSlopeToBasisMoments1D(std::vector<double> const & moment_dof, std::vector<double> & left_right_moments)
{
	if (moment_dof.size() != 2)
	{
		std::cerr << "Passed in incorrect length of moment vector to FEMUtilties::convertAvgSlopeToBasis\n";
		std::cerr << "Length = " << moment_dof.size() << std::endl;
		exit(1);
	}
	left_right_moments.resize(2);
	left_right_moments[0] = moment_dof[0] - moment_dof[1] / 3.; //left basis moment
	left_right_moments[1] = moment_dof[0] + moment_dof[1] / 3.; //right basis moment
}