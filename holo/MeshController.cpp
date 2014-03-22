//
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : MeshController.cpp
//  @ Date : 2/23/2014
//  @ Author : SRB
//
//


#include "MeshController.h"

//Functor for comparing pairs, needed for sorting jump errors
class CompareBySecond
{
public:
	bool operator () (const std::pair<int, double> & left, const std::pair<int, double> & right)
	{
		return (left.second > right.second);
	}
};

MeshController::MeshController(HoMesh* mesh, double exp_conv_constant, int n_batches_to_check):
	_required_conv_rate(std::fmax(0.0,exp_conv_constant)), //error must at least stay constant
	_mesh(mesh), 
	_n_batches_to_check(n_batches_to_check),
	_batch_residual_norms()
{
	createConnectivityArray();
}

double MeshController::computeJumpIntegral(std::vector<double> avg_slope_1, std::vector<double> avg_slope_2, double width)
{
	if (avg_slope_1.size() != 2 && avg_slope_2.size() != 2)
	{
		std::cerr << "MeshController::computeJumpIntegral only implemented for 1D\n";
		exit(12);
	}

	//define coefficients for integral
	double diff_avgs = std::abs(avg_slope_2[0] - avg_slope_1[0]);
	double diff_slopes = std::abs(avg_slope_2[1] - avg_slope_1[1]);

	//check if there is a sign change
	if (diff_avgs > diff_slopes || diff_slopes < GlobalConstants::RELATIVE_TOLERANCE)  //no slope, or no jump
	{
		return width*diff_avgs;
	}
	else //sign change occured, integral becomes:
	{
		return 0.5*width*(diff_avgs*diff_avgs/diff_slopes + diff_slopes);
	}
}

double MeshController::computeElementJumpError(int element_id)
{
	//Computed values
	std::vector<double> element_jump_errors;

	//Get element information
	ECMCElement1D* element = _mesh->getElement(element_id);
	std::vector<double> dimens = element->getElementDimensions();
	std::vector<double> coords = element->getElementCoordinates();
	double & mu = coords[1];
	double & x = coords[0];
	double & h_x = dimens[0];
	double & h_mu = dimens[1];
	unsigned int refinement_lev = element->getRefinementLevel();
	std::vector<double> el_dof = element->getAngularFluxDOF();
	std::vector<double> el_face_dof(2, 0.0);

	//Initialize adjacent element information
	ECMCElement1D* adj_element;
	std::vector<double> adj_dimens(2,0.);
	std::vector<double> adj_coords(2,0.);
	std::vector<double> adj_dof;
	std::vector<double> adj_face_dof(2,0.0);
	
	//Values of -1,1,or 0 for mapping to correct faces in a general way
	double avg_mu_val, avg_x_val;
	double slope_mu_val, slope_x_val;
	double edge_width;

	//loop over neighbors
	ElementNeighbors neighbors_struct = _connectivity_array[element_id];
	std::vector<ECMCElement1D*> neighbors = _connectivity_array[element_id]._element_pntrs;
	for (int i_nbr = 0; i_nbr < neighbors.size(); i_nbr++)
	{
		adj_element = neighbors[i_nbr];
		if (adj_element == NULL)
		{
			continue;
		}
		
		//Set values to map fluxes to the correct faces
		if (adj_element == neighbors_struct._left)
		{
			//For example
			avg_mu_val = 0.0; //there is no term for mu added to average moment on face
			avg_x_val = -1.; //Average is psi_a - psi_x on left face
			slope_mu_val = 1.; //In the slope moment you are integrating the jump along mu, thus slope moment is psi_mu
			slope_x_val = 0.0; //not integrating along x face
			edge_width = h_mu; //integrating along x, so the width of the integration domain along the face is in mu
		}
		else if (adj_element == neighbors_struct._right)
		{
			avg_mu_val = 0.0;
			avg_x_val = 1.;
			slope_mu_val = 1.;
			slope_x_val = 0.0;
			edge_width = h_mu;
		}
		else if (adj_element == neighbors_struct._minus)
		{
			avg_mu_val = -1.0;
			avg_x_val = 0.0;
			slope_mu_val = 0.0;
			slope_x_val = 1.0;
			edge_width = h_x;
		}
		else if (adj_element == neighbors_struct._plus)
		{
			avg_mu_val = 1.0;
			avg_x_val = 0.0;
			slope_mu_val = 0.0;
			slope_x_val = 1.0;
			edge_width = h_x;
		}
		else //neighbors list has extra term
		{
			std::cerr << "Logic Error in compute jump error in MeshController";
		}

		//Compute integrals based on relative refinement levels
		if (adj_element->hasChildren())
		{
			//adj_element is more refined, compute jump error for both children on the face
			double jump_error = 0.0;
			edge_width *= 0.5;
			std::vector<double> el_child_dof(el_dof); //copy to manipulate for each of the children

			//map el_dof slopes to smaller element children
			el_child_dof[1] *= 0.5;
			el_child_dof[2] *= 0.5;

			//Initialize to the bottom element, nearest to appropriate face, for horizontal edges will be "left" child
			el_child_dof[0] += avg_x_val*el_child_dof[1] + avg_mu_val*el_child_dof[2] - //get to face edge this line, next line get to top or bottom cell 
				(slope_mu_val*el_child_dof[2] + slope_x_val*el_child_dof[1]); //map to smaller element, i know this works, but very convoluted

			//Initialize the location of the child element
			adj_coords = adj_element->getElementCoordinates();
			adj_dimens = adj_element->getElementDimensions();
			double & adj_mu = adj_coords[1];
			double & adj_x = adj_coords[0];
			double & adj_h_x = adj_dimens[0];
			double & adj_h_mu = adj_dimens[1];
			double adj_child_x = adj_x - 0.25*adj_h_x*(avg_x_val + std::abs(avg_mu_val)); //(neighbor:child) right:left, left:right, plus:left, minus:left
			double adj_child_mu = adj_mu - 0.25*adj_h_mu*(avg_mu_val+std::abs(avg_x_val));	//(neighbor:child) plus:minus, , plus:left, minus:left
			
			for (int i_child = 0; i_child < 2; i_child++)
			{
				//get the child element values
				ECMCElement1D* child(adj_element->getChild(adj_child_x, adj_child_mu)); //copy constructor
				adj_dof = child->getAngularFluxDOF();

				//set average values on face
				el_face_dof[0] = el_child_dof[0] + avg_x_val*el_child_dof[1] + avg_mu_val*el_child_dof[2];
				adj_face_dof[0] = adj_dof[0] - avg_x_val*adj_dof[1] - avg_mu_val*adj_dof[2];

				//set slope values on face;
				el_face_dof[1] = slope_x_val*el_child_dof[1] + slope_mu_val*el_child_dof[2];
				adj_face_dof[1] = slope_x_val*adj_dof[1] + slope_mu_val*adj_dof[2];

				jump_error += computeJumpIntegral(el_face_dof, adj_face_dof, edge_width);

				//move el_dof average to the second child element (either top or right element)
				el_child_dof[0] += 2*(el_child_dof[2] * slope_mu_val + el_child_dof[1] * slope_x_val); //move to top if integrating along mu, right if integrating along x
				adj_child_x += slope_x_val*0.5*adj_h_x; //if integrating along x, move to next position
				adj_child_mu += slope_mu_val*0.5*adj_h_mu; //if integrating along mu, move to next position
			}
			element_jump_errors.push_back(jump_error);
		}
		else if (adj_element->getRefinementLevel() == refinement_lev)
		{
			//same refinement level
			adj_dof = adj_element->getAngularFluxDOF(); 

			//set average values on face
			el_face_dof[0] = el_dof[0] + avg_x_val*el_dof[1] + avg_mu_val*el_dof[2];
			adj_face_dof[0] = adj_dof[0] - avg_x_val*adj_dof[1] - avg_mu_val*adj_dof[2];
			
			//set slope values on face;
			el_face_dof[1] = slope_x_val*el_dof[1] + slope_mu_val*el_dof[2];
			adj_face_dof[1] = slope_x_val*adj_dof[1] + slope_mu_val*adj_dof[2];

			element_jump_errors.push_back(computeJumpIntegral(el_face_dof, adj_face_dof, edge_width));
		}
		else if (adj_element->getRefinementLevel() == refinement_lev - 1) //adjacent element is less refined
		{
			//map adj_element dof to smaller element
			adj_dof = adj_element->getAngularFluxDOF();
			adj_coords = adj_element->getElementCoordinates();
			double & adj_mu = adj_coords[1];
			double & adj_x = adj_coords[0];
			double & adj_h_x = adj_dimens[0];
			double & adj_h_mu = adj_dimens[1];
			adj_dof[1] *= 0.5; //slope is half on smaller element
			adj_dof[2] *= 0.5;

			//determine if the "top" or "bottom" element, for horizontal edges will be "right" or "left"
			double top_or_bottom = (mu*slope_mu_val + x*slope_x_val > adj_mu*slope_mu_val + adj_x*slope_x_val ? 1. : -1.);
			adj_dof[0] = adj_dof[0] - avg_mu_val*adj_dof[2] - avg_x_val*adj_dof[1] + top_or_bottom*(slope_mu_val*adj_dof[2] + slope_x_val*adj_dof[1]); //map to smaller element, i know this works, but very convoluted

			//set average values on face
			el_face_dof[0] = el_dof[0] + avg_x_val*el_dof[1] + avg_mu_val*el_dof[2];
			adj_face_dof[0] = adj_dof[0] - avg_x_val*adj_dof[1] - avg_mu_val*adj_dof[2];

			//set slope values on face;
			el_face_dof[1] = slope_x_val*el_dof[1] + slope_mu_val*el_dof[2];
			adj_face_dof[1] = slope_x_val*adj_dof[1] + slope_mu_val*adj_dof[2];

			element_jump_errors.push_back(computeJumpIntegral(el_face_dof, adj_face_dof, edge_width));
		}
		else
		{
			std::cerr << "Refinement logic error in MesHController.cpp::computeElementJumpError()\n";
			exit(1);
		}
	}
	
	if (HoController::USE_MAX_JUMP_ERROR)
	{
		return *(std::max_element(element_jump_errors.begin(), element_jump_errors.end()));
	}
	else
	{
		std::cerr << "Only max jump error indicator implemented, in MeshController.cpp\n";
		exit(1);
	}
}

void MeshController::createConnectivityArray()
{
	//loop parameters
	std::vector<ECMCElement1D*>* elements = _mesh->getElements();
	std::vector<ECMCElement1D*>::iterator it_el;

	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		if ((*it_el)->hasChildren())
		{
			continue;
		}	
		_connectivity_array[(*it_el)->getID()] = findNeighbors((*it_el)->getID());
	}
}

void MeshController::storeResidualNorm(double residual_norm)
{
	if (_batch_residual_norms.size() == (_n_batches_to_check+1)) //one extra since you need ratios to check refinement
	{
		_batch_residual_norms.erase(_batch_residual_norms.begin()); //remove the oldest element (first in)
	}
	_batch_residual_norms.push_back(residual_norm);
}

void MeshController::refineMesh()
{ 
	//compute the jump error for each active element
	std::vector<std::pair<int,double>> jump_errors; //a vector of <cell_id, jump_error_value>
	jump_errors.resize(_connectivity_array.size()); //initialize to the number of active elements

	std::vector<ECMCElement1D*>::iterator it_el;
	size_t index = 0;
	for (it_el = _mesh->_elements.begin(); it_el != _mesh->_elements.end(); it_el++)
	{
		if ((*it_el)->hasChildren()) //element is inactive
		{
			continue;
		}
		else
		{
			int element_id = (*it_el)->getID();
			jump_errors[index].first = element_id;
			jump_errors[index].second = computeElementJumpError(element_id);
			index++;
		}
	}

	//sort jump errors, biggest to smallest
	std::sort(jump_errors.begin(), jump_errors.end(), CompareBySecond());

	//Determine number of new elements to create
	int n_refinements = std::ceil(HoController::FRACTION_CELLS_TO_REFINE*(double)_mesh->_n_elems); //round up so always at least one
	
	//Refine elements, and neighboring elements as needed
	index = 0;  //which element to refine

	
	while (true)
	{
		refineElement(jump_errors[index].first);
		if (_newly_refined_elements.size() >= n_refinements) //limit on number of refinements, including extra refinements to keep mesh regular
		{
			break;
		}
		else
		{
			index++;
		}
	}

	//update the boundary cells if needed
	_mesh->findUpwindBoundaryCells();

	//update the connectivity array for all the new elements and their neighbors
	for (int el_id = 0; el_id < _newly_refined_elements.size(); el_id++)
	{
		updateConnectivityArray(_newly_refined_elements[el_id]);
	}
	_newly_refined_elements.clear(); //no need to keep this
	
	//reset convergence rate criteria
	_batch_residual_norms.erase(_batch_residual_norms.begin(), _batch_residual_norms.end() - 1); //clear all but the last one
}

bool MeshController::meshNeedsRefinement()
{
	if (_batch_residual_norms.size() != (_n_batches_to_check + 1))
	{
		return false; //not enough batches to check convergence, because of noise
	}
	else
	{
		double alpha_avg = 0.0;
		for (int i = 0; i < _n_batches_to_check; ++i) //loop over batches
		{
			alpha_avg += std::log(_batch_residual_norms[i] / _batch_residual_norms[i+1]);
		}
		alpha_avg /= (float)_n_batches_to_check;

		//check convergence
		if (alpha_avg < _required_conv_rate) //if error increases, alpha will be negative, and thus less than required rate
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

void MeshController::refineElement(int element_id)
{

	//Check that the element to be refined has not already been refined to maintain mesh regularity, in this case another call to this function has already refined the cell
	std::vector<int>::iterator it_refined = std::find(_newly_refined_elements.begin(),
		_newly_refined_elements.end(), element_id);
	if (it_refined != _newly_refined_elements.end())
	{
		return;
	}

	//add refined element to the array
	_newly_refined_elements.push_back(element_id);

	//probably better to have this function return the new elements made
	ECMCElement1D* unrefined_element = _mesh->_elements[element_id];
	unsigned int refinement_level = unrefined_element->getRefinementLevel();
	std::vector<ECMCElement1D*> new_elements;

	//Check that all neighbors are at least same level of refinement, if not, refine those elements recursively
	std::vector<ECMCElement1D*> element_neighbors(_connectivity_array[element_id]._element_pntrs);
	std::vector<ECMCElement1D*>::iterator it_nbr;
	for (it_nbr = element_neighbors.begin(); it_nbr != element_neighbors.end(); it_nbr++)
	{
		if (*it_nbr != NULL)
		{
			if ((*it_nbr)->getRefinementLevel() < refinement_level)
			{
				if (!(*it_nbr)->hasChildren()) //this will handle case that this elemetn has already been refined but neighbor list hasnt been updated yet
				{
					refineElement((*it_nbr)->getID());
					updateConnectivityArray((*it_nbr)->getID());
				}
			}
		}
	}

	//refine the element
	unrefined_element->refine(_mesh->_n_elems - 1); //pass the id of last element made

	//add the new elements to the list and update number of elements
	new_elements = _mesh->getElement(element_id)->getChildren();
	_mesh->_elements.insert(_mesh->_elements.end(), new_elements.begin(), new_elements.end());
	_mesh->_n_elems += new_elements.size();

	//update the number of active elements
	_mesh->_n_active_elements += new_elements.size() - 1; //subtract one for the element you just refined

	//determine if the refined cell was on a boundary
	std::vector<int>::iterator it_bound;
	it_bound = std::find(_mesh->_boundary_cells.begin(), _mesh->_boundary_cells.end(), element_id);
	if (it_bound != _mesh->_boundary_cells.end()) 
	{
		_mesh->_boundary_cells_need_update = true; //created cell was on a boundary, update boundary cells
		return; //no upwind boundary poitners to update
	}
	
	//Update downstream element pointers to the refined elements if necessary
	ECMCElement1D* upstream_el = _mesh->findJustUpwindElement(element_id);
	if (upstream_el->getRefinementLevel() == new_elements[0]->getRefinementLevel()) 
	{
		std::cerr << "This is old code when findJustUpwindElement used to return the refined cells, it should never be called, except for potentially on a boundary cell, in MeshController::refineElement()\n";
		exit(1);
		//the upstraem element is of same refinement, need to set the other daughter too
		if (upstream_el->getAngularCoordinate() < unrefined_element->getAngularCoordinate()) //found bottom (minus mu) refined cell
		{
			upstream_el->setDownStreamElement(new_elements[1]);
			upstream_el = _mesh->findJustUpwindElement(new_elements[3]->getID()); //get the top cell
			upstream_el->setDownStreamElement(new_elements[3]); //set the top cells ds element
		}
		else //you found the top (plus mu) refined cell
		{
			//This case is not normally called, but it works and may happen in the event of round off
			upstream_el->setDownStreamElement(new_elements[3]);
			upstream_el = _mesh->findJustUpwindElement(new_elements[1]->getID());
			upstream_el->setDownStreamElement(new_elements[1]);
		}
	}
	else if (upstream_el->getRefinementLevel() == unrefined_element->getRefinementLevel()) 
	{
		if (upstream_el->hasChildren())
		{
			//This is called because the boundary cells are only being updated once per total refinement,
			//but it is not a problem.  The findJustUpwind function works in this case, it just returns
			//the parent element rather than the children because children are not listed as boundary cell yet
			upstream_el->getChild(0)->setDownStreamElement(new_elements[1]);
			upstream_el->getChild(2)->setDownStreamElement(new_elements[3]);
		}
		else
		{
			//no ds elem pnters to update
		}
	}
}

ElementNeighbors MeshController::findNeighbors(int element_id)
{
	ElementNeighbors neighbors;

	//loop parameters
	std::vector<ECMCElement1D*>* elements = _mesh->getElements();
	std::vector<ECMCElement1D*>::iterator it_el;
	
	//reference element parameters
	ECMCElement1D* element = (*elements)[element_id];
	int el_refinement = (int)element->getRefinementLevel();
	double el_mu_coor = element->getAngularCoordinate();
	double el_mu_minus = el_mu_coor - 0.5*element->getAngularWidth();
	double el_mu_plus = el_mu_coor + 0.5*element->getAngularWidth();
	double el_x_coor = element->getSpatialCoordinate();
	double el_x_left = el_x_coor - 0.5*element->getSpatialWidth();
	double el_x_right = el_x_coor + 0.5*element->getSpatialWidth();


	//Set the downstream and upstream element values, based on flow direction
	_mesh->findUpwindBoundaryCells(); //must have boundary cells for this function to work
	ECMCElement1D* up_str_element = _mesh->findJustUpwindElement(element_id);
	ECMCElement1D*  ds_element = element->getDownStreamElement();
	if (el_mu_coor > 0)
	{
		neighbors._right = ds_element;
		neighbors._left = up_str_element;
	}
	else
	{
		neighbors._right = up_str_element;
		neighbors._left = ds_element;
	}

	//find the top (plus mu) and bottom (minus mu) elements
	bool top_found = false;
	bool bottom_found = false;
	for (it_el = elements->begin(); it_el != elements->end(); it_el++)
	{
		//check that still have elements to find
		if (bottom_found) 
		{
			if (top_found)
			{
				break;
			}
		}

		//don't compare to itself
		if (*it_el == element) continue;

		//make sure element is a possible match, based on refinement
		int it_refinement_level = (int)(*it_el)->getRefinementLevel();
		if (std::abs(it_refinement_level - el_refinement) > 1)
		{
			continue;
		}

		double it_x_coor = (*it_el)->getSpatialCoordinate();
		double it_h_x = (*it_el)->getSpatialWidth();

		if (std::abs(el_x_coor - it_x_coor) < 0.5*it_h_x) //element matches vertically
		{
			double it_mu_coor = (*it_el)->getAngularCoordinate();
			double it_h_mu = (*it_el)->getAngularWidth();
			double it_mu_minus = it_mu_coor - 0.5*it_h_mu;
			double it_mu_plus = it_mu_coor + 0.5*it_h_mu;

			if (std::abs(el_mu_plus - it_mu_minus) < GlobalConstants::RELATIVE_TOLERANCE)
			{
				//check that refinement level isnt greater, if that is the case then the parent cell is top cell and will be found elsewhere
				if (it_refinement_level > el_refinement)
				{
					continue;
				}
				else //found top cell
				{
					neighbors._plus = (*it_el);
					top_found = true;
				}
			}
			else if (std::abs(el_mu_minus - it_mu_plus) < GlobalConstants::RELATIVE_TOLERANCE)
			{
				//check that refinement level isnt greater, if that is the case then the parent cell of it_el is bot cell and will be found elsewhere in loop
				if (it_refinement_level > el_refinement)
				{
					continue;
				}
				else //found top cell
				{
					neighbors._minus = (*it_el);
					bottom_found = true;
				}
			}
			else
			{
				continue;
			}
		}
	}

	if (! HoController::REFINE_ACROSS_MU_ZERO)
	{
		//check that element does not have a mu edge on zero
		if (std::abs(el_mu_minus) < GlobalConstants::RELATIVE_TOLERANCE)
		{
			neighbors._minus = NULL;
		}
		else if (std::abs(el_mu_plus) < GlobalConstants::RELATIVE_TOLERANCE)
		{
			neighbors._plus = NULL;
		}
	}
	
	return neighbors;
}

void MeshController::updateConnectivityArray(int refined_element_id)
{
	//Check all neighbors to see if they need updating before you update the actual elements
	ElementNeighbors neighbors = _connectivity_array[refined_element_id];
	for (int i = 0; i<neighbors._element_pntrs.size(); ++i)
	{
		ECMCElement1D* neighbor = neighbors._element_pntrs[i];
		if (neighbor != NULL)
		{
			if (neighbor->hasChildren()) //need to update pntrs for its children, else pointers still point to the parent of refined_element_id which is ok
			{
				std::vector<ECMCElement1D*> children(neighbor->getChildren());
				std::vector<ECMCElement1D*>::iterator it_child;
				for (it_child = children.begin(); it_child != children.end(); it_child++)
				{
					if ((*it_child)->hasChildren()) //since the neighbors lists can be outdated, it is possible some refined elements are still in list
					{
						_connectivity_array.erase((*it_child)->getID());
					}
					else
					{
						int child_id = (*it_child)->getID();
						_connectivity_array[child_id] = findNeighbors(child_id);
					}
				}
			}
		}
	}

	//Now update the connectivity array for the children of this element
	std::vector<ECMCElement1D*> children = _mesh->_elements[refined_element_id]->getChildren();
	std::vector<ECMCElement1D*>::iterator it_child;
	for (it_child = children.begin(); it_child != children.end(); it_child++)
	{
		int child_id = (*it_child)->getID();
		_connectivity_array[child_id] = findNeighbors(child_id);
	}

	//Delete the parent cells member in map
	_connectivity_array.erase(refined_element_id);
}