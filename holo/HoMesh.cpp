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

HoMesh::HoMesh(Mesh* lo_mesh, int n_ang_cells_half_range)
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
					x_center = x_center_start;
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

std::vector<int> HoMesh::findUpwindBoundaryCells() const
{
	if (_lo_mesh->getSpatialDimension() != 1)
	{
		std::cerr << "HoMesh::findBoundaryCells only works on 1D meshes currently" << std::endl;
		system("pause");
		exit(1);
	}

	std::vector<int> boundary_cells;
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
				boundary_cells.push_back(i);
			}
		}
	}
	
	//exclude boundary cells that have an upwind cell, i.e., cells that are down wind of another cell. This is for the refined case
	//THIS CODE WORKS BUT REALLY NEEDS TO BE DONE IN A CLEANER MANNER
	int down_wind_id = -9999;
	ECMCElement1D* down_wind_el;
	std::vector<int>::iterator exclusion_id;
	std::vector<int> final_boundary_cells;
	final_boundary_cells = boundary_cells; //copy over
	for (int bc_id=0; bc_id < boundary_cells.size(); bc_id++)
	{
		down_wind_el = _elements[boundary_cells[bc_id]]->getDownStreamElement();
		if (down_wind_el != NULL) //special case for when there is only one spatial element in LO mesh
		{
			down_wind_id = down_wind_el->getID();
			if (down_wind_el->hasChildren()) //need to remove all children from options since they are not on boundary
			{
				std::vector<ECMCElement1D*> children = down_wind_el->getChildren();
				for (int child = 0; child < children.size(); child++)
				{
					exclusion_id = std::find(final_boundary_cells.begin(), final_boundary_cells.end(), children[child]->getID()); //find if the downwind cell in boundary list
					if (exclusion_id != final_boundary_cells.end())
					{
						final_boundary_cells.erase(exclusion_id); //found id, erase that memeber
					}
				}
			}
		} 
		exclusion_id = std::find(final_boundary_cells.begin(), final_boundary_cells.end(), down_wind_id); //is the downwind cell in boundary list?
		if (exclusion_id != final_boundary_cells.end()) 
		{
			final_boundary_cells.erase(exclusion_id); //found id, erase that memeber
		}
	}
	return final_boundary_cells;
}

std::vector<DirichletBC1D*> HoMesh::getDirichletBCs() const
{
	return _lo_mesh->getDirichletBCs();
}