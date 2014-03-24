//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : Element1D.cpp
//  @ Date : 11/1/2013
//  @ Author : 
//
//

#include <iostream>
#include <iomanip>
#include "Element1D.h"
#include "AverageCosineData.h"
#include "GlobalConstants.h"
#include "DirichletBC1D.h"

using namespace std;

Element1D::Element1D() : Element(), _phi_left_plus(_elem_dofs[DOF_MAP[0]]), _phi_right_plus(_elem_dofs[DOF_MAP[1]]),
_phi_left_minus(_elem_dofs[DOF_MAP[2]]), _phi_right_minus(_elem_dofs[DOF_MAP[3]])
{
	using std::cerr;
	//Default constructor, this will never actually be called
	cerr << "This should not be called" << std::endl;
	exit(1);

}

//Initialization list is done purely for references to alias base class dof members
Element1D::Element1D(int id, MaterialConstant* mat, std::vector<Node* >  nodes) : Element(1, id, mat, nodes),
_phi_left_plus(_elem_dofs[DOF_MAP[0]]), _phi_right_plus(_elem_dofs[DOF_MAP[1]]),
_phi_left_minus(_elem_dofs[DOF_MAP[2]]), _phi_right_minus(_elem_dofs[DOF_MAP[3]])
{
	using EquationMaps1D::DOF_MAP; //This is for mapping degrees of freedom to their aliases

	//Make sure element nodes are ordered left to right and compute width of element
	_h = _nodes[1]->getX() - _nodes[0]->getX();
	if (_h < 0.0)
	{
		std::cerr << "Coordinates are not properly ordered" << std::endl;
		exit(0);
	}
	_lo_data = new LoData1D(); //default for constructor is diffusion theory values
	_mat = mat;
	_n_elem_dof = 4; //This is hard coded for 1D system
	_ext_source_nodal_values.resize(2);

	//Initialize the global degree numbers based on element ID.
	for (int i = 0; i < _n_elem_dof; ++i)
		_elem_dofs[i].setEqnNum(_id * _n_elem_dof + i);
}


void Element1D::print(std::ostream &out) const
{
	using std::endl;
	out << " ID = ";
	out.width(4);
	out << _id;
	out << ", Width = ";
	std::setprecision(4);
	out << _h;
	out << ", Nodes: ";
	for (int i = 0; i < _n_nodes; i++){
		out.width(4);
		out << _nodes[i]->getID();
	}
	out << endl;
}

void Element1D::setElementDof(std::vector<double> elem_dofs) 
{
	if (elem_dofs.size() != _elem_dofs.size())
	{
		cout << "This function needs to be checked";
		exit(1);
	}
	for (int i = 0; i < _n_elem_dof; ++i)
		_elem_dofs[i].setValue(elem_dofs[i]);

}

std::vector<int> Element1D::getEqnNumbers(void) const
{
	std::vector<int> eqn_numbers;
	for (int i = 0; i < _n_elem_dof; ++i)
		eqn_numbers.push_back(_elem_dofs[i].getEqnNum());
	return eqn_numbers;
}

void Element1D::getElementMomentMatrix(numMatrix* M, numVector* b, std::vector<int> &eqns) const
{
	using EquationMaps1D::DOF_MAP;
	using EquationMaps1D::EQN_MAP;

	//Get out constants
	using GlobalConstants::FOUR_PI;

	//Local variables
	double alpha = _lo_data->getSpatialClosureFactor();
	AveragedCosines vol_mu, surf_mu; 
	double value, sigma_t, sigma_s;
	//local row and column for each term are based on DOF map

	//Set local variables
	vol_mu = _lo_data->getVolAveragedCos();
	surf_mu = _lo_data->getSurfAveragedCos();
	sigma_t = _mat->getSigmaT();
	sigma_s = _mat->getSigmaS();

	//DOF_MAP is used to add coefficients to appropriate equations
	//DOF_MAP[0] refers to < . >_L^+ term
	//DOF_MAP[1] refers to < . >_R^+ term
	//DOF_MAP[2] refers to < . >_L^- term
	//DOF_MAP[3] refers to < . >_R^- term
	//
	//EQN_MAP is used to add terms to appropriate row, same ordering as DOF_MAP
	
	//Same source term for + and - equations for each spatial moment eqn (_L and _R)
	double source_value_left = _h/FOUR_PI * (2. / 3.*_ext_source_nodal_values[0] + 1. / 3.*_ext_source_nodal_values[1]);
	double source_value_right = _h/FOUR_PI * (1. / 3.*_ext_source_nodal_values[0] + 2. / 3.*_ext_source_nodal_values[1]);

	//Add in terms for < . >_L^+ equation
	// The order is changed here, because we want this equation to always be first to limit band size of matrix
	//----------------------------------------------------------------------------------------------------
	value = vol_mu._mu_left_plus + (_h * (sigma_t - sigma_s / (FOUR_PI)));
	M->addCoeff(EQN_MAP[0], DOF_MAP[0], value);
	value = vol_mu._mu_right_plus;
	M->addCoeff(EQN_MAP[0], DOF_MAP[1], value);
	value = -_h* (sigma_s / FOUR_PI);
	M->addCoeff(EQN_MAP[0], DOF_MAP[2], value);
	M->addCoeff(EQN_MAP[0], DOF_MAP[3], 0.0);
	b->addCoeff(EQN_MAP[0], source_value_left);  

	//Add in terms for < . >_R^+ equation
	//----------------------------------------------------------------------------------------------------
	value = 2.0*  surf_mu._mu_right_plus *(1 - alpha) - vol_mu._mu_left_plus;
	M->addCoeff(EQN_MAP[1], DOF_MAP[0], value);
	value = (2. * surf_mu._mu_right_plus * alpha) - vol_mu._mu_right_plus + sigma_t*_h - (sigma_s*_h / FOUR_PI);
	M->addCoeff(EQN_MAP[1], DOF_MAP[1], value);
	M->addCoeff(EQN_MAP[1], DOF_MAP[2], 0.0);
	value = -_h* (sigma_s / FOUR_PI);
	M->addCoeff(EQN_MAP[1], DOF_MAP[3], value);
	b->addCoeff(EQN_MAP[1], source_value_right); 

	//Add in terms for < . >_L^- equation
	//----------------------------------------------------------------------------------------------------
	value = -_h* (sigma_s / FOUR_PI);
	M->addCoeff(EQN_MAP[2], DOF_MAP[0], value);
	M->addCoeff(EQN_MAP[2], DOF_MAP[1], 0.0);
	value = vol_mu._mu_left_minus + sigma_t*_h - (sigma_s*_h / FOUR_PI) - 2.*surf_mu._mu_left_minus*alpha;
	M->addCoeff(EQN_MAP[2], DOF_MAP[2], value);
	value = -2.*surf_mu._mu_left_minus*(1. - alpha) + vol_mu._mu_right_minus;
	M->addCoeff(EQN_MAP[2], DOF_MAP[3], value);
	b->addCoeff(EQN_MAP[2], source_value_left);  

	//Add in terms for < . >_R^- equation
	// The order is changed here, because we want this equation to always be last to limit band size of matrix
	//----------------------------------------------------------------------------------------------------
	M->addCoeff(EQN_MAP[3], DOF_MAP[0], 0.0);
	value = -_h* (sigma_s / FOUR_PI);
	M->addCoeff(EQN_MAP[3], DOF_MAP[1], value);
	value = -1.*vol_mu._mu_left_minus;
	M->addCoeff(EQN_MAP[3], DOF_MAP[2], value);
	value = sigma_t*_h - (sigma_s*_h / FOUR_PI) - vol_mu._mu_right_minus;
	M->addCoeff(EQN_MAP[3], DOF_MAP[3], value);
	b->addCoeff(EQN_MAP[3], source_value_right);  //TODO need to make this not a constant

	eqns = getEqnNumbers(); //set the equation numbers correctly
	M->printCompact(std::cout);
}

void Element1D::getPosUpwinding(std::vector<double> &pos_upwind_values, int &eqn, std::vector<int> &cols) const
{
	using EquationMaps1D::EQN_MAP;

	//Local variables
	double mu_times_2 = 2.0*_lo_data->getSurfAveragedCos()._mu_right_plus;
	double alpha = _lo_data->getSpatialClosureFactor();

	//Positive upwinding terms go to < . >_L^+ equation, for next cell
	//TODO probably could clean this up
	eqn = EQN_MAP[0] + _n_elem_dof*(_id+1);

	//Resize passed in vectors
	pos_upwind_values.resize(2);
	cols.resize(2);

	//Add the phi_R^+ term, and append appropriate column number
	pos_upwind_values[0] = -1.*mu_times_2*alpha;
	cols[0] = _phi_right_plus.getEqnNum();

	//Add the phi_L+ term, and append appropriate column number
	pos_upwind_values[1] = -1.*mu_times_2*(1. - alpha);
	cols[1] = _phi_left_plus.getEqnNum();
}

void Element1D::getNegUpwinding(std::vector<double> &neg_upwind_values, int &eqn, std::vector<int> &cols) const
{
	using EquationMaps1D::EQN_MAP;

	//Local variables
	double mu_times_2 = 2.0*_lo_data->getSurfAveragedCos()._mu_left_minus;
	double alpha = _lo_data->getSpatialClosureFactor();

	//Negative upwinding terms go to < . >_R^- equation, for previous cell
	eqn = EQN_MAP[3] + _n_elem_dof*(_id - 1);

	//Resize passed in vectors
	neg_upwind_values.resize(2);
	cols.resize(2);

	//Add the phi_L^- term, and append appropriate column number
	neg_upwind_values[0] = mu_times_2*alpha;
	cols[0] = _phi_left_minus.getEqnNum();

	//Add the phi_R^- term, and append appropriate column number
	neg_upwind_values[1] = mu_times_2*(1. - alpha);
	cols[1] = _phi_right_minus.getEqnNum();
}

void Element1D::addDirichletBC(numVector* b, std::vector<int> &eqns, double dirichlet_value, Node* node) const
{
	eqns = getEqnNumbers(); //TODO this will need to change if not passing in whole element b
	using EquationMaps1D::EQN_MAP;

	//Here, it is performed for a Marshak Boundary Condition, where "dirichlet_value" is incoming current
	//and does not to be weighted by mu

	if (node->getID() == _nodes[0]->getID())  //Is it a left BC?
	{
		//Term gets added to the phi_L+ equation
		b->zero();
		b->addCoeff(EQN_MAP[0], 2.*dirichlet_value);
	}
	else if (node->getID() == _nodes[1]->getID()) //Is it a right BC?
	{
		//Term gets added to the phi_R- equations 
		b->zero();
		b->addCoeff(EQN_MAP[3], 2.*abs(dirichlet_value));
	}
	else
	{
		std::cerr << "Your boundary condition does not correspond to this element" << std::endl;
	}
}

//Locations is the x coordinate of the nodes
void Element1D::getScalarFluxValues(double alpha, 
	std::vector<double> &scalar_flux_values, std::vector<double> &locations) const   //Return the flux and location of the nodes.  In 1D locations is just x coordinates
{

	if (scalar_flux_values.size() != 2) //Resize
	{
		scalar_flux_values.resize(2);
		locations.resize(2);
	}

	//Compute left flux and right scalar fluxes
	double _phi_right_minus_value = _phi_right_minus.getValue();
	double _phi_left_minus_value = _phi_left_minus.getValue();
	double _phi_right_plus_value = _phi_right_plus.getValue();
	double _phi_left_plus_value = _phi_left_plus.getValue();

	scalar_flux_values[0] = alpha*(_phi_left_minus_value + _phi_left_plus_value)
		+ ( (1. - alpha) * (_phi_right_minus_value + _phi_right_plus_value) ) ;
	scalar_flux_values[1] = alpha*(_phi_right_minus_value + _phi_right_plus_value)
		+ ((1. - alpha) * (_phi_left_minus_value + _phi_left_plus_value));

	//Store left and right locations of nodes, the nodes are guaranteed to be in teh right order when they were created, but double check
	locations[0] = _nodes[0]->getX();
	locations[1] = _nodes[1]->getX();

	//TODO, if code is slow, you can remove this check
	if (locations[0] > locations[1])
	{
		std::cerr << "Your nodal mesh points are out of order, error in Element1D.cpp Line 271" << std::endl;
		exit(1);
	}
}

void Element1D::getScalarFluxLinDisc(std::vector<double> &scalar_flux_values, std::vector<double> &locations) const //Returns scalar flux on faces using Linear Discontinuous approximation
{
	getScalarFluxValues(2.0, scalar_flux_values, locations);  //For lin dis. alpha = 2
}

std::vector<double> Element1D::getElementDimensions() const
{
	std::vector<double> dimensions;
	dimensions.push_back(_h);
	return dimensions;
}

void Element1D::getScalarFluxHOClosure(std::vector<double> &scalar_flux_values,
	std::vector<double> &locations) const //This is for scalar flux based on alpha closure, shouldnt be used except verification
{
	getScalarFluxValues(_lo_data->getSpatialClosureFactor(), scalar_flux_values, locations);
}

//This function returns the source based on the LD closure
std::vector<double> Element1D::getScalarFluxNodalValues() const
{
	std::vector<double> flux_values;
	std::vector<double> dummy_vector; //locations is value of nodes, not needed here
	getScalarFluxLinDisc(flux_values, dummy_vector);
	return flux_values;
}

void Element1D::printLDScalarFluxValues(std::ostream &out) const
{
	std::vector<double> flux_values;
	std::vector<double> flux_coordinates;

	getScalarFluxLinDisc(flux_values, flux_coordinates);
	
	for (unsigned int i = 0; i < flux_coordinates.size(); ++i)
	{
		out.precision(3);
		out <<  flux_coordinates[i] << " " << setw(14) <<
			setprecision(15) << scientific << flux_values[i] << endl;
	}
}

std::vector<double> Element1D::getNodalCoordinates() const
{
	std::vector<double> locations;
	locations.resize(2);
	locations[0] = _nodes[0]->getX();
	locations[1] = _nodes[1]->getX();
	return locations;
}