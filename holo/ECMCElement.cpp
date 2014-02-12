#include "ECMCElement.h"
#include <iostream>

ECMCElement::ECMCElement(Element* spatial_element, std::vector<double> dimensions, std::vector<double> coordinates)
{
	_spatial_element = spatial_element;
	_dimensions = dimensions;
	_coordinates = coordinates;
	_ang_flux_dof.assign(_dimensions.size() + 1, 0.0); //1 extra for average flux, set all to 0
	_tally = NULL; //make derived classes create this

	if (_dimensions.size() != _coordinates.size())
	{
		std::cerr << "Vector sizes do not agree in ECMCElement constructor" << std::endl;
		exit(1);
	}
}

ECMCElement::ECMCElement(Element* spatial_element, int n_spatial_dimensions, int n_angular_dimensions)
{
	_spatial_element = spatial_element;
	_tally = NULL;

	//Initialize arrays for use with derived classes
	int n_dimensions = n_angular_dimensions + n_spatial_dimensions;
	_dimensions.resize(n_dimensions);
	_coordinates.resize(n_dimensions);
	_ang_flux_dof.resize(n_dimensions + 1); //1 extra for average flux
}

std::vector<double> ECMCElement::getElementDimensions() const
{
	return _dimensions;
}
