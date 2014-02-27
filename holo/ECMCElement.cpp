#include "ECMCElement.h"
#include <iostream>

ECMCElement::ECMCElement(int id, Element* spatial_element, std::vector<double> dimensions, std::vector<double> coordinates)
{
	_id = id;
	_spatial_element = spatial_element;
	_dimensions = dimensions;
	_coordinates = coordinates;
	_ang_flux_dof.assign(_dimensions.size() + 1, 0.0); //1 extra for average flux, set all to 0
	_tally = NULL; //make derived classes create this

	//mesh refinement info, leave vector uninitialized
	_has_children = false;
	_refinement_level = 0;

	if (_dimensions.size() != _coordinates.size())
	{
		std::cerr << "Vector sizes do not agree in ECMCElement constructor" << std::endl;
		exit(1);
	}
}

ECMCElement::ECMCElement(int id, Element* spatial_element, int n_spatial_dimensions, int n_angular_dimensions)
{
	_id = id;
	_spatial_element = spatial_element;
	_tally = NULL;

	//Initialize arrays for use with derived classes
	int n_dimensions = n_angular_dimensions + n_spatial_dimensions;
	_dimensions.resize(n_dimensions);
	_coordinates.resize(n_dimensions);
	_ang_flux_dof.resize(n_dimensions + 1); //1 extra for average flux

	//mesh refinement info, leave vector uninitialized
	_has_children = false;
	_refinement_level = 0;
}

std::vector<double> ECMCElement::getElementDimensions() const
{
	return _dimensions;
}

std::vector<double> ECMCElement::getAngularFluxDOF() const
{
	return _ang_flux_dof;
}

int ECMCElement::getID() const
{
	return _id;
}

unsigned int ECMCElement::getRefinementLevel() const
{
	return _refinement_level;
}