// 
//
//  Generated by StarUML(tm) C++ Add-In
//
//  @ Project : Untitled
//  @ File Name : HoMesh.cpp
//  @ Date : 2/10/2014
//  @ Author : 
//
//


#include "HoMesh.h"
#include "GlobalConstants.h"

HoMesh::HoMesh(Mesh* lo_mesh, int n_ang_cells_half_range) :
	_boundary_cells(),
	_boundary_cells_need_update(true)
{
	_lo_mesh = lo_mesh; //stored for finding boundary cells

	if (n_ang_cells_half_range < 1)
	{
		n_ang_cells_half_range = 1; //must have at least 1 cell in each direction
		std::cout << "Creating Ho mesh with 2 elements..." << std::endl;
	}

	std::vector<Element*>* spatial_elements;
	spatial_elements = lo_mesh->getElements();
	Element* spat_elem;        //pointer to current spatial element

	//TODO this is an angular connectivity array, currently not used, but correctly constructed
	std::vector<std::vector<int>> connectivity_array;
	connectivity_array.resize(lo_mesh->getNumElems()); //each member is a list of ECMC elements that correspond to each spatial cell


	//create the number of angular elements needed for both half ranges
	if (lo_mesh->getSpatialDimension() == 1)
	{
		//Local dimensions to pass into element constructors
		std::vector<double> dimensions, coordinates; //width_X,width_mu, X_i, mu_i
		dimensions.resize(2);
		coordinates.resize(2);
		double mu_center=0.0; //the center of each element in mu
		double x_center;  // the center of each element in x
		ECMCElement1D* down_stream_element = NULL;
		ECMCElement1D* last_element = NULL;
		int element_id=-1;
		
		//all angular widths are the same before refinement
		dimensions[1] = 1. / (double)(n_ang_cells_half_range);
		mu_center = -1 - 0.5*dimensions[1]; //starting value for mu_center (ghost cell)

		//create the elments for negative flow direction first, currently makes in reverse order because it is easiest
		//to create the elements upstream, and the order in the array doesnt actually matter
		for (int i_ang = 0; i_ang < n_ang_cells_half_range; i_ang++)
		{	
			mu_center += dimensions[1];
			coordinates[1] = mu_center;
			
			for (int it_el = 0; it_el<spatial_elements->size(); it_el++)
			{	
				element_id++;
				spat_elem = (*spatial_elements)[it_el];
				if (it_el == 0) //Edge elements have a null pointer
				{
					x_center = -0.5*spat_elem->getElementDimensions()[0];
					down_stream_element = NULL;
				}
				else
				{
					down_stream_element = last_element;
				}
				
				//determine dimensions for the element
				dimensions[0] = spat_elem->getElementDimensions()[0]; //width_x
				x_center += dimensions[0]; //move the center forward each time
				coordinates[0] = x_center;

				//create the element
				last_element = new ECMCElement1D(element_id, spat_elem, down_stream_element, dimensions, coordinates);
				_elements.push_back(last_element);
				connectivity_array[it_el].push_back(_elements.size() - 1); //add last element to correct connectivity array
			}
		}

		//store x_center at the end for easy access next time
		double x_center_start = x_center;

		//Construct elements for positive flow direction, traverse loop up stream
		for (int i_ang = 0; i_ang < n_ang_cells_half_range; i_ang++)
		{
			mu_center += dimensions[1];
			coordinates[1] = mu_center;
			for (int it_el = spatial_elements->size()-1; it_el >= 0; it_el--)
			{
				element_id++;
				spat_elem = (*spatial_elements)[it_el];
				if (it_el == spatial_elements->size()-1) //Edge elements have a null pointer
				{
					x_center = x_center_start+spat_elem->getElementDimensions()[0]; //start at ghost cell;
					down_stream_element = NULL;
				}
				else
				{
					down_stream_element = last_element;
				}

				//determine dimensions for the element
				dimensions[0] = spat_elem->getElementDimensions()[0]; //width_x
				x_center -= dimensions[0]; //move the center backward each time, moving to the left
				coordinates[0] = x_center;

				//create the element
				last_element = new ECMCElement1D(element_id, spat_elem, down_stream_element, dimensions, coordinates);
				_elements.push_back(last_element);
				connectivity_array[it_el].push_back(_elements.size() - 1); //add last element to correct connectivity array
			}
		}
	}
	else
	{
		std::cerr << "Multi dimension not implemented for HoMesh constructor\n";
		system("pause");
		exit(1);
	}



	_n_elems = _elements.size();
}

ECMCElement1D* HoMesh::getElement(int element_id) const
{
	return _elements[element_id];
}

int HoMesh::getNumElems() const
{
	if (_n_elems != _elements.size())
	{
		std::cerr << "Error in updating element sizes in HoMesh.cpp\n";
		exit(1);
	}
	return _n_elems;
}

std::vector<ECMCElement1D* >* HoMesh::getElements(void)
{
	return &_elements; //pointer, be careful not to change these values
}

void HoMesh::computeAngularFluxes(int n_histories, double total_src_strength)
{
	for (int i = 0; i < _n_elems; ++i)
	{
		_elements[i]->computeAngularFLuxDOF(n_histories, total_src_strength);
	}
}

void HoMesh::printAngularFluxes(std::ostream &out)
{
	using std::endl;

	out << "\n---------------------------------------------------------\n"
		<< "                   Angular Flux Moments\n"
		<< "---------------------------------------------------------\n";

	for (int i = 0; i < _n_elems; ++i)
	{
		out << "ID = " << i;
		_elements[i]->printAngularFluxDOF(out);
	}
}

std::vector<int> HoMesh::findUpwindBoundaryCells() 
{
	if (_lo_mesh->getSpatialDimension() != 1)
	{
		std::cerr << "HoMesh::findBoundaryCells only works on 1D meshes currently" << std::endl;
		system("pause");
		exit(1);
	}
	if (!_boundary_cells_need_update)
	{
		return _boundary_cells;
	}

	_boundary_cells.clear();
	int spatial_ID;
	int last_element_ID = _lo_mesh->getNumElems() - 1;

	for (int i = 0; i < _n_elems; ++i)
	{
		if (_elements[i]->hasChildren()) //skip parent cells
		{
			continue;
		}
		spatial_ID = _elements[i]->getSpatialElement()->getID();
		if ( spatial_ID == 0 || spatial_ID == last_element_ID) //spatial cell on a boundary, need to make sure boundary is upwind
		{
			if (_elements[i]->getDownStreamElement() != NULL || last_element_ID == 0) // potential boundary element, id=0 case is for 1 cell geometry
			{
				//check if element is on edge of boundary
				std::vector<double> x_edges;
				x_edges = _elements[i]->getSpatialElement()->getNodalCoordinates();
				double x_coor = _elements[i]->getSpatialCoordinate();
				double h_x = _elements[i]->getSpatialWidth();
				if (_elements[i]->getAngularCoordinate() > 0.0)
				{
					double edge_left = x_coor - 0.5*h_x;
					if (std::abs((edge_left - x_edges[0]) / h_x) < GlobalConstants::RELATIVE_TOLERANCE) 
					{ 
						_boundary_cells.push_back(i); //on left boundary face
					}
				}
				else //negative flow
				{
					double edge_right = x_coor + 0.5*h_x;
					if (std::abs((edge_right - x_edges[1]) / h_x) < GlobalConstants::RELATIVE_TOLERANCE)
					{
						_boundary_cells.push_back(i); //on right boundary face
					}
				}
			}
		} 
	}

	//reset the boundary_cells flag and return new list of boundary cells
	_boundary_cells_need_update = false;
	return _boundary_cells;
}

std::vector<DirichletBC1D*> HoMesh::getDirichletBCs() const
{
	return _lo_mesh->getDirichletBCs();
}

ECMCElement1D* HoMesh::findJustUpwindElement(int down_str_element_id)
{
	//find the boundary cell that is on the same mu level as current element
	std::vector<int>::iterator it_bc_id;
	findUpwindBoundaryCells(); //update just in case

	double bc_mu_center;
	double ds_mu_center = _elements[down_str_element_id]->getAngularCoordinate();
	double ds_h_mu = _elements[down_str_element_id]->getAngularWidth();
	double ds_mu_minus = ds_mu_center - 0.5*ds_h_mu;
	double ds_mu_plus = ds_mu_center + 0.5*ds_h_mu;
	unsigned int ds_refinement_level = _elements[down_str_element_id]->getRefinementLevel();

	for (it_bc_id = _boundary_cells.begin(); it_bc_id != _boundary_cells.end(); it_bc_id++)
	{
		bc_mu_center = _elements[*it_bc_id]->getAngularCoordinate();
		if (bc_mu_center > ds_mu_minus) //os abp
		{
			if (bc_mu_center < ds_mu_plus)
			{
				break;
			}
		}	
	}
	if (it_bc_id == _boundary_cells.end())
	{
		std::cerr << "Boundary element not found in findJustUpwindElement, HoMesh.cpp\n";
		exit(1);
	}

	//make sure that the boundary element not a child of the ds element (if it has any), since the
	//boundary cells have been updated before this function is called, children are active
	//and thus in boundary cell list
	if (_elements[down_str_element_id]->hasChildren())
	{
		std::vector<ECMCElement1D*> ds_element_children = _elements[down_str_element_id]->getChildren();
		ECMCElement1D* up_str_element = _elements[*it_bc_id];
		std::vector<ECMCElement1D*>::iterator it_children;
		it_children = std::find(ds_element_children.begin(), ds_element_children.end(), up_str_element);
		if (it_children != ds_element_children.end())
		{
			return NULL; //the element has children on the boundary
		}
	}
	else
	{
		if (*it_bc_id == down_str_element_id)
		{
			return NULL; //the element is on a boundary
		}
	}
	
	//Else, there is an interior cell upwind of this element
	//start on boundary and search to find the upwind element, note 
	//that this function will work correctly even if the upwind boundary elements
	//have yet to be updated, this saves computing the new boundary elements every
	//time you add a new element
	ECMCElement1D* new_down_str_elem = NULL;
	ECMCElement1D* up_str_elem = _elements[*it_bc_id]; //start on boundary
	while (up_str_elem != NULL)
	{
		new_down_str_elem = up_str_elem->getDownStreamElement();
		if (new_down_str_elem->getID() == down_str_element_id)
		{
			break;
		}
		else if ( ds_refinement_level >new_down_str_elem->getRefinementLevel())
		{
			//need to move back down a refinement level when possible (at some point it will be possible)
			if (new_down_str_elem->hasChildren())
			{
				new_down_str_elem = new_down_str_elem->findChildEntered(ds_mu_center + GlobalConstants::RELATIVE_TOLERANCE*ds_h_mu); //add tiny delta h to prevent roundoff errors in this function, shoudl be find though

				//check that refined cell is not the droid you were looking for
				if (new_down_str_elem->getID() == down_str_element_id)
				{
					break;
				}
				else
				{
					up_str_elem = new_down_str_elem;
				}
			}
			else //a later cell will be able to go back to the correct refinement level
			{
				up_str_elem = new_down_str_elem;
			}
		}
		else
		{
			up_str_elem = new_down_str_elem; //check next element
		}
	}
	if (up_str_elem == NULL)
	{
		std::cerr << "Could not find downstream element, in HoMesh.cpp::finJustUpwind\n";
		exit(1);
	}
	return up_str_elem;
}