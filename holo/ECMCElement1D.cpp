//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : ECMCElement1D.cpp
//  @ Date : 2/10/2014
//  @ Author : 
//
//

#include "GlobalConstants.h"
#include "ECMCElement1D.h"
#include "FEMUtilities.h"

ECMCElement1D::ECMCElement1D(int id, Element* element, ECMCElement1D* down_stream_element, 
		std::vector<double> dimensions, std::vector<double> coordinates, int refinement_level):
	ECMCElement(id, element,dimensions,coordinates),
	_psi_average(_ang_flux_dof[0]), _psi_x(_ang_flux_dof[1]), _psi_mu(_ang_flux_dof[2]),
	_width_spatial(_dimensions[0]), _width_angle(_dimensions[1]),
	_position_center(_coordinates[0]), _mu_center(_coordinates[1])
{
	_refinement_level = refinement_level; //overwrite refinement level in the base class
	_down_stream_element = down_stream_element;
	_tally = new ECMCTally(); //intializes the correct size
}

ECMCElement1D* ECMCElement1D::getDownStreamElement() const
{
	return _down_stream_element;
}

void ECMCElement1D::setDownStreamElement(ECMCElement1D* new_ds_elem)
{
	_down_stream_element = new_ds_elem;
}

double ECMCElement1D::getSpatialWidth() const
{
	return _width_spatial;
}

double ECMCElement1D::getAngularWidth() const
{
	return _width_angle;
}

double ECMCElement1D::getAngularCoordinate() const
{
	return _mu_center;
}

double ECMCElement1D::getSpatialCoordinate() const
{
	return _position_center;
}

void ECMCElement1D::incrementTallyScores(double weight, double path_length_cm, double dir_cosine,
	double normalized_position) //this method will normalize the direction cosine
{
	//Score Element tally, need to convert path_length and volume
	//to cm, rather than mfp
	double volume_cm_str = _width_spatial*_width_angle;//*1.0cm*1.0cm*delta_mu = h_xh_mu(cm^3-str)
	double normalized_direction = ((dir_cosine - _mu_center) / _width_angle) + 0.5; //normalized angle with in the element
	normalized_position *= 6.;

	_tally->incrementScores(weight, path_length_cm, normalized_direction, volume_cm_str, normalized_position);
}

void ECMCElement1D::computeAngularFLuxDOF(int n_histories, double & l2_error_element_sq, double total_src_strength) 
{
	if (_has_children) //in the future may want to map from fine elements on to parent elements
	{
		l2_error_element_sq = 0.0;
		return;
	}

	std::vector<double> spatial_moments = _tally->getSpatialMoments(n_histories); //0th and 1st spatial moment
	double angular_moment = _tally->getAngularMoment(n_histories); //angular moment

	//References for the additive error DOF
	std::vector<double> err_dof(3);
	double & psi_avg_err = err_dof[0];
	double & psi_x_err = err_dof[1];
	double & psi_mu_err = err_dof[2];

	//calculate moments based on LD closure, adding the moments from the tallies
	//if standard MC, tallies are the the angular flux, else tallies are the additive error
	psi_avg_err = spatial_moments[0] * total_src_strength;
	psi_x_err = 6. * (spatial_moments[1] - 0.5*spatial_moments[0])*total_src_strength;
	psi_mu_err = 6. * (angular_moment - 0.5*spatial_moments[0])*total_src_strength;

	//DEBUG STATEMENTS, THESE ARE NOT NECESSARILY PROBLEMS AS LONG AS THE RESIDUAL IS ZERO IN THIS ELEMENT
	/*
	if (std::abs(psi_avg_err) < GlobalConstants::RELATIVE_TOLERANCE*_psi_average)
	{
		std::cout << "Average flux is too small, Element ID = "<< _id << " " << _psi_average << "\n";
		printData(std::cout);
	}
	if (std::abs(psi_x_err) < GlobalConstants::RELATIVE_TOLERANCE *std::abs(_psi_x))
	{
		std::cout << "X flux is too small, Element ID = " << _id << " " << psi_x_err << "PSIX_IS " << _psi_x << " AVERAGE_IS:  " << _psi_average << std::endl;
		printData(std::cout);
	}
	if (std::abs(psi_mu_err) < GlobalConstants::RELATIVE_TOLERANCE *std::abs(_psi_mu))
	{
		std::cout << "Mu flux is too small, Element ID = " << _id << " " << psi_mu_err << " " << _psi_mu << std::endl;
	}*/

	//update angular flux
	_psi_average += psi_avg_err;
	_psi_x += psi_x_err;
	_psi_mu += psi_mu_err;

	//reset tallies to zero
	_tally->reset();

	//compute l2 relative error over element with quadrature
	GaussQuadrature quad;
	int n_qps = quad.getNumPoints();
	std::vector<double> x_wgts(quad.getQuadratureWeights(_position_center, _width_spatial));
	std::vector<double> mu_wgts(quad.getQuadratureWeights(_mu_center, _width_angle));
	std::vector<double> x_pnts(quad.getQuadraturePoints(_position_center, _width_spatial));
	std::vector<double> mu_pnts(quad.getQuadraturePoints(_mu_center, _width_angle));

	l2_error_element_sq = 0.0; //zero out sum
	double err_sq_qp; //error in angular flux evaluated at quadrature point
	for (int i_qp = 0; i_qp < n_qps; ++i_qp) //x_qps
	{
		for (int j_qp = 0; j_qp < n_qps; ++j_qp) //mu_qps
		{
			err_sq_qp = FEMUtilities::evalLinDiscFunc2D(err_dof, _dimensions, _coordinates, x_pnts[i_qp], mu_pnts[j_qp]);
			err_sq_qp *= err_sq_qp; //square it
			l2_error_element_sq += x_wgts[i_qp] * mu_wgts[j_qp]*err_sq_qp; //w_x_i w_mu_i * psi^2(x_i,mu_i)
		}
	}
}

void ECMCElement1D::printAngularFluxDOF(std::ostream &out) const
{
	using std::ios;

	if (_has_children)
	{
		out << " Parent Element\n";
		return;
	}

	out << " psi avg. = ";
	out.setf(ios::scientific);
	out.precision(15);
	out.width(12);
	out << _psi_average << " psi x. = ";
	out.width(12);
	out << _psi_x << " psi mu. =";
	out.width(12);
	out << _psi_mu << std::endl;
}

//-----------------------------------------------------------------------------------
//Adaptive refinement functions
//-----------------------------------------------------------------------------------

void ECMCElement1D::refine(int last_element_id)
{
	if (_has_children)
	{
		std::cout << "Cannot refine elements that already have children\n";
		exit(1);
	}
	_children.clear(); //clear out the children list

	//local variables
	double sgn_mu = abs(_mu_center) / _mu_center; //+1 or -1
	ECMCElement1D* child;
	ECMCElement1D* child_ds_elem;
	std::vector<double> child_dimens(2, NULL);
	std::vector<double> child_coor(2, NULL);
	double & child_h_mu = child_dimens[1];
	double & child_h_x = child_dimens[0];
	double & child_x_coor = child_coor[0];
	double & child_mu_coor = child_coor[1];

	//initialize variables before loops
	child_h_mu = _width_angle*0.5;
	child_h_x = _width_spatial*0.5;
	child_mu_coor = _mu_center - child_h_mu*0.5; //intialize to first row of element
	int child_id = last_element_id + 1; //initialized to last element+1

	//update refinement 
	_has_children = true;

	//cells are created from downwind to upwind, then from minus mu to plus mu
	for (int i = 0; i < 2; ++i) //start with downwind cells
	{
		//initialize to downstream element as less refined or NULL
		child_ds_elem = _down_stream_element;
		if (_down_stream_element != NULL) 
		{
			if (_down_stream_element->hasChildren()) //downwind cell is refined, so get the up wind element of the row
			{
				std::vector<ECMCElement1D*> down_str_children = _down_stream_element->getChildren();
				child_ds_elem = down_str_children[2 * i + 1];
			}
		}

		child_x_coor = _position_center + sgn_mu*child_h_x*0.5;
		child = new ECMCElement1D(child_id, _spatial_element, child_ds_elem,
			child_dimens,child_coor,_refinement_level+1);  		//create first element in row
		_children.push_back(child);
		child_id++; //increment id
		
		//create second element in row
		child_x_coor -= sgn_mu*child_h_x;  //update x_coor
		child_ds_elem = child; //the first child is downstream of the next one
		child = new ECMCElement1D(child_id, _spatial_element, child_ds_elem,
			child_dimens, child_coor, _refinement_level + 1);  		
		_children.push_back(child);
		child_id++;
		
		child_mu_coor += child_h_mu; //update mu for next row
	}

	//map current angular flux estimates to the new cells
	mapAngFluxToChildren();

	//delete info that is no longer important for the parent class
	delete _tally; //will also ensure no errors in tracking routines
}

std::vector<ECMCElement1D*> ECMCElement1D::getChildren() const
{
	if (_has_children)
	{
		return _children;
	}
	else
	{
		std::cerr << "Tried to return children, but this element is not refined, in ECMCElement1D.cpp\n";
		exit(1);
	}
}

void ECMCElement1D::mapAngFluxToChildren() 
{
	ECMCElement1D* child;
	double sgn_x, sgn_mu; 
	for (int child_id = 0; child_id < _children.size(); child_id++)
	{
		child = _children[child_id];
		child->_psi_mu = _psi_mu*0.5; //same for all elements
		child->_psi_x = _psi_x*0.5;  //same for all elements

		//determien which quadrant cell is in
		sgn_x = (child->_position_center > _position_center ? 1. : -1.);
		sgn_mu = (child->_mu_center > _mu_center ? 1. : -1.);

		//calculate psi_average based on quadrant and child's slopes
		child->_psi_average = _psi_average + sgn_x*child->_psi_x + sgn_mu*child->_psi_mu;
	}
}

ECMCElement1D*  ECMCElement1D::findChildEntered(double mu) const
{
	//since particle is moving downstream, must be either 1 or 3 in child vector, since 
	//the elements are built upstream
	return (mu < _mu_center) ? _children[1] : _children[3];	
}

ECMCElement1D* ECMCElement1D::getChild(double x, double mu)
{
	//find child based on an x and mu coordinate
	if ((std::abs(mu - _mu_center) > _width_angle*0.5) ||
		(std::abs(x - _position_center) > _width_spatial*0.5))
	{
		std::cerr << "Passed in x and mu that are not within this element, in ECMCElement1D::getChild(x,mu)\n";
		exit(1);
	}

	int row, column;
	if (mu < _mu_center) //On bottom
	{
		row = 0;
	}
	else
	{
		row = 1;
	}
	if (x <= _position_center) //on left side
	{
		if (_mu_center > 0.0)
		{
			column = 1;
		}
		else
		{
			column = 0;
		}
	}
	else //on right side
	{
		if (_mu_center > 0.0)
		{
			column = 0;
		}
		else
		{
			column = 1;
		}
	}
	return _children[2 * row + column];
}

ECMCElement1D* ECMCElement1D::getChild(int index) const
{
	if (index > _children.size() - 1)
	{
		std::cerr << "Tried to access a child outside of possible indices, in ECMCElement1D.cpp\n";
		exit(1);
	}
	return _children[index];
}